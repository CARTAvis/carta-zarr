/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The one-pass histogram's arithmetic, against chosen values and exact counts.
//
// Every assertion here used to be a tolerance over a fixture: the class was file-local, so reaching
// it meant ComputeCubeHistogram over a directory tree, compared against a two-pass oracle to within
// a fraction of a bin. Its two fixes were each a defect in one of the four things below and neither
// could be stated in those terms: 9c21397 gave a straddling bin whole to one side, and 6ba0ce2
// seeded the range at one rather than at the pixel it saw.
//
// This target links nothing. A growing histogram takes floats and gives back counts.

#include "reduce/growing_histogram.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "support/check.h"

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
// provisional ones -- which is what makes these exact rather than approximate.
void TestTheRangeSeedsAndGrowsBothWays() {
    GrowingHistogram histogram(4);

    // Seeded around 8: [0, 16) in four bins of 4, so this lands in [8, 12).
    histogram.Add(8.0F);
    RequireCounts(histogram.Aggregate(4, 0.0, 16.0), {0, 0, 1, 0}, "the pixel the range was seeded around");

    // Below the range: it doubles downward to [-16, 16) in bins of 8, and the 8 that was in bin 2 of
    // four is now in bin 3 of four. -1 lands in [-8, 0).
    histogram.Add(-1.0F);
    RequireCounts(histogram.Aggregate(4, -16.0, 16.0), {0, 1, 0, 1}, "after growing downward");

    // Above it: the lower bound stays and the range doubles to [-16, 48) in bins of 16, merging the
    // two halves of the old range into the first two bins. 20 lands in [16, 32).
    histogram.Add(20.0F);
    RequireCounts(histogram.Aggregate(4, -16.0, 48.0), {1, 1, 1, 0}, "after growing upward");
}

// A first pixel of exactly zero has no magnitude to seed a range from, and a range of zero width
// never grows: the loop that widens it to fit the next pixel would not terminate. The floor is one
// std::numeric_limits<float>::min(), which costs about a hundred and thirty merges once.
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
// Three pixels in the provisional bin [8, 12), aggregated onto [6, 10) and [10, 14): the bin is
// halved by the edge, so each side is owed 1.5. Each takes its whole 1, and the one leftover goes to
// the first of two equal remainders.
void TestAStraddlingBinIsSplitAndTheRemainderPlaced() {
    GrowingHistogram histogram(4);
    histogram.Add(8.0F);
    histogram.Add(9.0F);
    histogram.Add(10.0F);
    RequireCounts(histogram.Aggregate(4, 0.0, 16.0), {0, 0, 3, 0}, "all three should bin together first");

    RequireCounts(histogram.Aggregate(2, 6.0, 14.0), {2, 1}, "a bin straddling a target edge");
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
        TestAFirstPixelOfZeroStillGrows();
        TestAStraddlingBinIsSplitAndTheRemainderPlaced();
        TestAnEmptyTargetRangePutsEverythingInTheFirstBin();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "growing histogram test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
