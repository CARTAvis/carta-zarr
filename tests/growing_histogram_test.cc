/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The one-pass histogram's arithmetic, against chosen values and exact counts.
//
// Every assertion here is an exact count over a few pixels, not a tolerance against a two-pass
// oracle: how a bin straddling a target edge is split, and where the range seeds itself, are each a
// defect this can state in its own terms.
//
// This target links nothing. A growing histogram takes floats and gives back counts.

#include "reduce/growing_histogram.h"

#include "support/check.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace {

using carta::zarr::internal::GrowingHistogram;

using carta::zarr::testing::Require;

std::string Show(const std::vector<std::uint64_t>& counts) {
    std::string text = "{";
    for (std::size_t i = 0; i < counts.size(); ++i) {
        text += (i == 0 ? "" : ", ") + std::to_string(counts.at(i));
    }
    return text + "}";
}

void RequireCounts(const std::vector<std::uint64_t>& got, const std::vector<std::uint64_t>& want,
                   const std::string& what) {
    Require(got == want, what + ": expected " + Show(want) + ", got " + Show(got));
}

// The range seeds itself around the first pixel and then only ever grows, so the two directions are
// different operations: growing upward leaves the lower bound where it is, growing downward moves
// it by a whole range. Both merge in pairs, and neither may drop a count.
//
// The grid here is chosen so that aggregating onto it is the identity -- the target bins are the
// provisional ones -- which is what makes these exact rather than approximate. The first pixel is
// 2^24, where a float's spacing is 2, so the seed is four bins of 2 centred on it.
void TestTheRangeSeedsAndGrowsBothWays() {
    GrowingHistogram histogram(4);
    constexpr double kAt = 16777216.0;

    // Seeded around 2^24: [2^24 - 4, 2^24 + 4) in four bins of 2, so this lands in [2^24, 2^24 + 2).
    histogram.Add(static_cast<float>(kAt));
    RequireCounts(histogram.Aggregate(4, kAt - 4, kAt + 4), {0, 0, 1, 0}, "the pixel the range was seeded around");

    // Below the range: it doubles downward to [2^24 - 12, 2^24 + 4) in bins of 4, and the pixel that
    // was in bin 2 of four is now in bin 3 of four. 2^24 - 5 lands in [2^24 - 8, 2^24 - 4).
    histogram.Add(static_cast<float>(kAt - 5));
    RequireCounts(histogram.Aggregate(4, kAt - 12, kAt + 4), {0, 1, 0, 1}, "after growing downward");

    // Above it: the lower bound stays and the range doubles to [2^24 - 12, 2^24 + 20) in bins of 8,
    // merging the two halves of the old range into the first two bins. 2^24 + 10 lands in
    // [2^24 + 4, 2^24 + 12).
    histogram.Add(static_cast<float>(kAt + 10));
    RequireCounts(histogram.Aggregate(4, kAt - 12, kAt + 20), {1, 1, 1, 0}, "after growing upward");
}

// Pixels close together far from zero keep the resolution between them. The first range is as
// narrow as the floats around the first pixel are apart, and grows only as far as the pixels need;
// seeded as wide as the pixel's own magnitude, it was 2000 wide around 1000 and never shrank, so a
// pixel a sixteenth above thirty-nine others shared their provisional bin and was spread evenly over
// all four target bins.
void TestANarrowDistributionFarFromZeroKeepsItsResolution() {
    GrowingHistogram histogram(64);
    for (int i = 0; i < 39; ++i) {
        histogram.Add(1000.0F);
    }
    histogram.Add(1000.0625F);
    RequireCounts(histogram.Aggregate(4, 1000.0, 1000.0625), {39, 0, 0, 1},
                  "a narrow distribution around a nonzero baseline");
}

// At the largest float there is no float above it to take the spacing from -- the next one up is
// infinity, and an infinitely wide bin puts every pixel in the first one, by a conversion of
// infinity to an index that is undefined. The spacing below is the same there, and is the one taken.
void TestTheLargestFloatsKeepTheirResolution() {
    constexpr float kLargest = std::numeric_limits<float>::max();
    const float below = std::nextafter(kLargest, 0.0F);
    GrowingHistogram histogram(64);
    for (int i = 0; i < 39; ++i) {
        histogram.Add(kLargest);
    }
    histogram.Add(below);
    RequireCounts(histogram.Aggregate(4, below, kLargest), {1, 0, 0, 39}, "pixels at the largest float");

    // The most negative float seeds from the same spacing. Its pixels a float apart share one
    // provisional bin as wide as the whole target range, so they are spread as any bin is; what is
    // asked of them is that the range is finite and none is lost.
    GrowingHistogram negative(64);
    for (int i = 0; i < 39; ++i) {
        negative.Add(-kLargest);
    }
    negative.Add(-below);
    std::uint64_t total = 0;
    for (const auto count : negative.Aggregate(4, -static_cast<double>(kLargest), -static_cast<double>(below))) {
        total += count;
    }
    Require(total == 40, "pixels at the most negative float were lost");
}

// A first pixel of exactly zero has no magnitude to seed a range from, and a range of zero width
// never grows: the loop that widens it to fit the next pixel would not terminate. The spacing of the
// floats at zero is the least denormal, which costs about a hundred and fifty merges once.
//
// So this test is mostly about arriving at all. What it asserts once it does is the invariant the
// whole design rests on: no pixel is lost, only resolution.
void TestAFirstPixelOfZeroStillGrows() {
    GrowingHistogram histogram(4);
    histogram.Add(0.0F);
    histogram.Add(1.0F);

    const auto counts = histogram.Aggregate(2, 0.0, 1.0);
    Require(counts.at(0) + counts.at(1) == 2, "both pixels should still be counted somewhere");
}

// A provisional bin that straddles a target edge is split in proportion to the overlap, and the
// pixels that do not divide evenly go to the largest remainders -- so the split still adds up.
//
// Three pixels of 2^24 + 4 seed four bins of 2 around themselves, [2^24, 2^24 + 8), and all land in
// [2^24 + 4, 2^24 + 6). Aggregated onto [2^24 + 3, 2^24 + 5) and [2^24 + 5, 2^24 + 7), that bin is
// halved by the edge, so each side is owed 1.5. Each takes its whole 1, and the one leftover goes to
// the first of two equal remainders.
void TestAStraddlingBinIsSplitAndTheRemainderPlaced() {
    GrowingHistogram histogram(4);
    constexpr double kAt = 16777216.0;
    for (int i = 0; i < 3; ++i) {
        histogram.Add(static_cast<float>(kAt + 4));
    }
    RequireCounts(histogram.Aggregate(4, kAt, kAt + 8), {0, 0, 3, 0}, "all three should bin together first");

    RequireCounts(histogram.Aggregate(2, kAt + 3, kAt + 7), {2, 1}, "a bin straddling a target edge");
}

// Every pixel had the same value, so the extremes the caller aggregates over are equal and there is
// exactly one bin anything can be in. This is the shape ComputeCubeHistogram hands over for a flat
// image, and dividing by the width of that range would be a division by zero.
void TestAnEmptyTargetRangePutsEverythingInTheFirstBin() {
    GrowingHistogram histogram(8);
    for (int i = 0; i < 5; ++i) {
        histogram.Add(3.0F);
    }
    RequireCounts(histogram.Aggregate(4, 3.0, 3.0), {5, 0, 0, 0}, "a range of one value");
}

}  // namespace

int main() {
    try {
        TestTheRangeSeedsAndGrowsBothWays();
        TestANarrowDistributionFarFromZeroKeepsItsResolution();
        TestTheLargestFloatsKeepTheirResolution();
        TestAFirstPixelOfZeroStillGrows();
        TestAStraddlingBinIsSplitAndTheRemainderPlaced();
        TestAnEmptyTargetRangePutsEverythingInTheFirstBin();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "growing histogram test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
