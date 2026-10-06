/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The layout a spectral reduction accumulates in, and the block it hands that over as.
//
// None of this needs a pixel: which slot a statistic lands in, what a slot starts as, how two
// partials combine, and what an untouched extremum is reported as are all questions about doubles
// in a buffer. So this links nothing but the header, and reads every answer through the one block a
// sink is ever given.

#include "reduce/statistic_slots.h"

#include "reduce/deviations.h"
#include "reduce/tuning.h"
#include "support/check.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using carta::zarr::SpectralBlock;
using carta::zarr::SpectralTotals;
using carta::zarr::Statistic;
using carta::zarr::StatisticSet;
using carta::zarr::internal::kStatisticOrder;
using carta::zarr::internal::RowTotals;
using carta::zarr::internal::Spread;
using carta::zarr::internal::StatisticLayout;
using carta::zarr::internal::StatisticSlots;

using carta::zarr::testing::Require;

// A set holds only statistics there are: no integer becomes one, however it is asked.
static_assert(!std::is_convertible_v<unsigned, StatisticSet> && !std::is_constructible_v<StatisticSet, unsigned>,
              "a StatisticSet must not be made from bits");
static_assert((Statistic::sum | Statistic::min).Contains(Statistic::min) &&
                  !StatisticSet(Statistic::sum).Contains(Statistic::min) && StatisticSet{}.empty(),
              "a set contains what was joined into it and nothing else");

constexpr StatisticSet kEverything = Statistic::num_pixels | Statistic::nan_count | Statistic::sum | Statistic::sum_sq |
                                     Statistic::min | Statistic::max | Statistic::sum_sq_dev;

RowTotals Row(std::uint64_t good, std::uint64_t bad, double sum, double sum_sq, double smallest, double largest) {
    RowTotals row;
    row.good = good;
    row.bad = bad;
    row.sum = sum;
    row.sum_sq = sum_sq;
    row.smallest = smallest;
    row.largest = largest;
    return row;
}

// What a sink saw, copied out, since the block's values are only valid for the duration of the call.
struct Seen {
    std::uint64_t first_channel = 0;
    std::uint64_t channel_count = 0;
    std::size_t region_count = 0;
    bool complete = false;
    double completeness = 0.0;
    std::vector<std::vector<SpectralTotals>> totals;  // [region][channel]
};

Seen HandOver(StatisticSlots& slots, bool complete = true, std::uint64_t first_channel = 0) {
    Seen seen;
    const bool kept_going =
        slots.HandOver(first_channel, complete, complete ? 1.0 : 0.5, [&](const SpectralBlock& block) {
            seen.first_channel = block.first_channel;
            seen.channel_count = block.channel_count;
            seen.region_count = block.region_count;
            seen.complete = block.complete;
            seen.completeness = block.completeness;
            seen.totals.resize(block.region_count);
            for (std::size_t r = 0; r < block.region_count; ++r) {
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    seen.totals.at(r).push_back(block.Totals(r, c));
                }
            }
            return true;
        });
    Require(kept_going, "a sink that says yes lets the reduction go on");
    return seen;
}

bool Same(double expected, double actual) {
    return std::isnan(expected) ? std::isnan(actual) : expected == actual;
}

void RequireTotals(const SpectralTotals& expected, const SpectralTotals& actual, const std::string& where) {
    Require(Same(expected.num_pixels, actual.num_pixels), where + ": num_pixels");
    Require(Same(expected.nan_count, actual.nan_count), where + ": nan_count");
    Require(Same(expected.sum, actual.sum), where + ": sum");
    Require(Same(expected.sum_sq, actual.sum_sq), where + ": sum_sq");
    Require(Same(expected.min, actual.min), where + ": min " + std::to_string(actual.min));
    Require(Same(expected.max, actual.max), where + ": max " + std::to_string(actual.max));
}

SpectralTotals Totals(double num_pixels, double nan_count, double sum, double sum_sq, double min, double max) {
    return SpectralTotals{num_pixels, nan_count, sum, sum_sq, min, max};
}

const SpectralTotals kNothing{};

// Every one of the 127 sets a request can make, in the one order a block reports them. Each is made
// from the statistics in it, bit i of `members` standing for kStatisticOrder[i]: a StatisticSet is
// never made from bits. A set with sum_sq_dev in it is laid out with num_pixels and sum as well, and
// keeps two more values per region than it reports: the mean its merges need.
void TestTheLayoutFollowsTheOneOrder() {
    for (unsigned members = 1; members < (1u << kStatisticOrder.size()); ++members) {
        StatisticSet requested;
        std::vector<Statistic> expected;
        const bool deviations = (members & (1u << (kStatisticOrder.size() - 1))) != 0;
        for (std::size_t i = 0; i < kStatisticOrder.size(); ++i) {
            if ((members & (1u << i)) != 0) {
                requested |= kStatisticOrder.at(i);
            }
            const bool brought = deviations && (kStatisticOrder.at(i) == Statistic::num_pixels ||
                                                kStatisticOrder.at(i) == Statistic::sum);
            if ((members & (1u << i)) != 0 || brought) {
                expected.push_back(kStatisticOrder.at(i));
            }
        }
        Require(kEverything.Contains(requested), "every set is part of the set of everything");
        const auto layout = StatisticLayout::Of(requested);
        const std::string set = "set " + std::to_string(members);
        Require(layout.count() == expected.size(), set + ": one slot per statistic asked for");
        for (std::size_t i = 0; i < expected.size(); ++i) {
            Require(layout.statistics()[i] == expected.at(i), set + ": slot " + std::to_string(i));
        }
        Require(layout.Deviations() == deviations, set + ": whether the loop takes deviations");
        Require(layout.BytesPerChannel(3) == 3 * (expected.size() + (deviations ? 2 : 0)) * sizeof(double),
                set + ": a channel costs one double per statistic per region, and two for a mean");
    }
}

// Nothing folded in yet reads as nothing counted: zero for the counts and sums, NaN for the extrema,
// because there is no smallest value of nothing.
void TestAResetHoldsNothing() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(kEverything), 2, 3);
    Require(slots.regions() == 2 && slots.channels() == 3, "a reset takes the shape it was given");

    const auto seen = HandOver(slots);
    for (std::size_t r = 0; r < 2; ++r) {
        for (std::uint64_t c = 0; c < 3; ++c) {
            RequireTotals(kNothing, seen.totals.at(r).at(c),
                          "untouched " + std::to_string(r) + "," + std::to_string(c));
        }
    }

    // And a reset to another shape starts again rather than keeping what was there.
    slots.Fold<true>(1, 2, Row(4, 1, 10.0, 30.0, 1.0, 4.0));
    slots.Reset(StatisticLayout::Of(kEverything), 3, 2);
    const auto again = HandOver(slots);
    Require(again.region_count == 3 && again.channel_count == 2, "the block has the new shape");
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::uint64_t c = 0; c < 2; ++c) {
            RequireTotals(kNothing, again.totals.at(r).at(c), "after a second reset");
        }
    }
}

void TestAFoldAddsTheCountsAndKeepsTheExtrema() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(kEverything), 2, 3);
    slots.Fold<true>(1, 2, Row(3, 1, 6.0, 14.0, 1.0, 3.0));
    slots.Fold<true>(1, 2, Row(2, 0, 9.0, 41.0, 4.0, 5.0));
    slots.Fold<true>(1, 2, Row(1, 2, -2.0, 4.0, -2.0, -2.0));

    const auto seen = HandOver(slots);
    RequireTotals(Totals(6, 3, 13.0, 59.0, -2.0, 5.0), seen.totals.at(1).at(2), "three rows into one channel");
    RequireTotals(kNothing, seen.totals.at(0).at(2), "the other region");
    RequireTotals(kNothing, seen.totals.at(1).at(1), "the next channel down");
}

// A statistic nobody asked for is neither accumulated nor reported, and asking for fewer does not
// move the ones that are.
void TestOnlyWhatWasAskedForIsCarried() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(Statistic::num_pixels | Statistic::max), 1, 2);
    slots.Fold<false>(0, 1, Row(3, 1, 6.0, 14.0, 1.0, 3.0));

    bool called = false;
    slots.HandOver(0, true, 1.0, [&](const SpectralBlock& block) {
        called = true;
        Require(block.Carries(Statistic::num_pixels | Statistic::max), "carries what was asked for");
        Require(block.Carries(Statistic::max), "and each part of it");
        Require(!block.Carries(Statistic::num_pixels | Statistic::sum), "but not a set with one more in it");
        Require(block.Carries(StatisticSet{}), "everything in the empty set is carried");
        Require(block.Series(0, Statistic::sum) == nullptr, "no series for a statistic not asked for");
        Require(block.Series(0, Statistic::min) == nullptr, "nor for an extremum not asked for");
        Require(block.Series(0, Statistic::num_pixels) != nullptr, "a series for one that was");

        // Read as nothing counted, which is what the absent value of each is worth.
        RequireTotals(Totals(3, 0, 0.0, 0.0, kNothing.min, 3.0), block.Totals(0, 1), "a partial set");
        return true;
    });
    Require(called, "the sink is called");
}

// Partials are merged in the order they are given, at the offset they are given, and nowhere else.
// The sums are compared exactly: the order is the determinism rule, so the same order must give the
// same bits.
void TestPartialsMergeInOrderAtTheirOffset() {
    const auto layout = StatisticLayout::Of(kEverything);
    std::vector<StatisticSlots> partials(3);
    for (auto& partial : partials) {
        partial.Reset(layout, 2, 4);
    }
    const std::array<double, 3> parts{0.1, 0.2, 0.3};
    for (std::size_t task = 0; task < 3; ++task) {
        for (std::uint64_t c = 0; c < 4; ++c) {
            const double v = parts.at(task) + static_cast<double>(c);
            partials.at(task).Fold<true>(1, c, Row(task + 1, task, v, v * v, v - 1.0, v + 1.0));
        }
    }

    StatisticSlots block;
    block.Reset(layout, 2, 12);
    block.MergeInOrder(partials.data(), partials.size(), 5);

    const auto seen = HandOver(block);
    for (std::uint64_t c = 0; c < 12; ++c) {
        RequireTotals(kNothing, seen.totals.at(0).at(c), "region 0 channel " + std::to_string(c));
        if (c < 5 || c >= 9) {
            RequireTotals(kNothing, seen.totals.at(1).at(c), "outside the merged channels " + std::to_string(c));
            continue;
        }
        const auto k = static_cast<double>(c - 5);
        double sum = 0.0;
        double sum_sq = 0.0;
        for (const double part : parts) {
            // Squared into a name first, as the fold above was handed it: written as one expression
            // the compiler may fuse the multiply into the add, which rounds once instead of twice.
            const double v = part + k;
            const double squared = v * v;
            sum += v;
            sum_sq += squared;
        }
        RequireTotals(Totals(6, 3, sum, sum_sq, parts.at(0) + k - 1.0, parts.at(2) + k + 1.0), seen.totals.at(1).at(c),
                      "merged channel " + std::to_string(c));
    }
}

// A block handed over before it is finished is handed over again, refined. The NaN an untouched
// extremum is reported as must not survive into the refinement, or the smallest value of anything
// would be NaN forever after -- std::min keeps whichever argument it was given first.
void TestAHandOverPutsTheIdentitiesBack() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(kEverything), 1, 2);
    slots.Fold<true>(0, 0, Row(1, 0, 5.0, 25.0, 5.0, 5.0));

    const auto early = HandOver(slots, false);
    Require(!early.complete && early.completeness == 0.5, "the block says it is unfinished");
    Require(std::isnan(early.totals.at(0).at(1).min), "an untouched minimum is reported as NaN");

    slots.Fold<true>(0, 1, Row(1, 0, 3.0, 9.0, 3.0, 3.0));
    slots.Fold<true>(0, 0, Row(1, 0, 7.0, 49.0, 7.0, 7.0));
    const auto late = HandOver(slots, true);
    RequireTotals(Totals(2, 0, 12.0, 74.0, 5.0, 7.0), late.totals.at(0).at(0), "refined");
    RequireTotals(Totals(1, 0, 3.0, 9.0, 3.0, 3.0), late.totals.at(0).at(1), "first touched after a hand-over");

    // Every hand-over puts them back, not only an unfinished one: what a finished block leaves
    // behind is the same accumulator an unfinished one does.
    slots.Fold<true>(0, 1, Row(1, 0, 1.0, 1.0, 1.0, 1.0));
    const auto after = HandOver(slots, true);
    RequireTotals(Totals(2, 0, 4.0, 10.0, 1.0, 3.0), after.totals.at(0).at(1), "folded into after a finished block");
}

// Rows folded into one channel, and partials merged into a block, come to the sum of squared
// deviations the pixels have all together -- the merges are the only place two sets of pixels meet.
namespace {

// What the per-pixel loop hands a fold for these pixels, taking its distances from the first.
RowTotals RowOf(const std::vector<double>& pixels) {
    double sum = 0.0;
    double sum_sq = 0.0;
    double distance_sum = 0.0;
    double distance_sum_sq = 0.0;
    double smallest = std::numeric_limits<double>::infinity();
    double largest = -std::numeric_limits<double>::infinity();
    for (const double pixel : pixels) {
        sum += pixel;
        sum_sq += pixel * pixel;
        distance_sum += pixel - pixels.front();
        distance_sum_sq += (pixel - pixels.front()) * (pixel - pixels.front());
        smallest = std::min(smallest, pixel);
        largest = std::max(largest, pixel);
    }
    auto row = Row(pixels.size(), 0, sum, sum_sq, smallest, largest);
    const auto spread =
        Spread::OfSpan(static_cast<double>(pixels.size()), pixels.front(), distance_sum, distance_sum_sq);
    row.base = spread.base;
    row.offset = spread.offset;
    row.sum_sq_dev = spread.sum_sq_dev;
    return row;
}

double DeviationsOf(const std::vector<std::vector<double>>& rows) {
    double count = 0.0;
    double sum = 0.0;
    for (const auto& row : rows) {
        for (const double pixel : row) {
            count += 1.0;
            sum += pixel;
        }
    }
    const double mean = sum / count;
    double deviations = 0.0;
    for (const auto& row : rows) {
        for (const double pixel : row) {
            deviations += (pixel - mean) * (pixel - mean);
        }
    }
    return deviations;
}

void RequireDeviations(double actual, double expected, const std::string& where) {
    Require(std::abs(actual - expected) <= 1e-12 * std::max(1.0, expected),
            where + ": sum_sq_dev " + std::to_string(actual) + ", the pixels' own " + std::to_string(expected));
}

const std::vector<std::vector<double>> kRows{
    {1.0e7, 1.0e7 + 1.0, 1.0e7 - 2.0}, {1.0e7 + 3.0, 1.0e7 + 3.5}, {9.0e6}, {}};

}  // namespace

void TestDeviationsFoldAsOneSetOfPixels() {
    StatisticSlots slots;
    // Asked for alone: the count and the sum it is merged with come with it.
    slots.Reset(StatisticLayout::Of(Statistic::sum_sq_dev), 1, 2);
    for (const auto& pixels : kRows) {
        slots.Fold<true>(0, 1, pixels.empty() ? RowTotals{} : RowOf(pixels));
    }
    const auto seen = HandOver(slots);
    RequireDeviations(seen.totals.at(0).at(1).sum_sq_dev, DeviationsOf(kRows), "four rows, one of them empty");
    Require(seen.totals.at(0).at(1).num_pixels == 6.0, "the count came with it");
    Require(seen.totals.at(0).at(0).sum_sq_dev == 0.0, "nothing folded is no spread");
}

void TestDeviationsMergeInOrderAsOneSetOfPixels() {
    const auto layout = StatisticLayout::Of(Statistic::sum_sq_dev | Statistic::max);
    std::vector<StatisticSlots> partials(3);
    for (std::size_t task = 0; task < partials.size(); ++task) {
        partials.at(task).Reset(layout, 1, 2);
        partials.at(task).Fold<true>(0, 0, RowOf(kRows.at(task)));
    }
    // A partial that saw nothing at a channel, and one that saw something only there.
    partials.at(2).Fold<true>(0, 1, RowOf({4.0, 5.0}));

    StatisticSlots block;
    block.Reset(layout, 1, 6);
    block.Fold<true>(0, 3, RowOf({1.0e7 + 10.0}));
    block.MergeInOrder(partials.data(), partials.size(), 3);
    const auto seen = HandOver(block);
    RequireDeviations(seen.totals.at(0).at(3).sum_sq_dev,
                      DeviationsOf({{1.0e7 + 10.0}, kRows.at(0), kRows.at(1), kRows.at(2)}),
                      "a channel the block already held, and three partials");
    RequireDeviations(seen.totals.at(0).at(4).sum_sq_dev, 0.5, "a channel only the last partial saw");
    Require(seen.totals.at(0).at(2).sum_sq_dev == 0.0 && seen.totals.at(0).at(5).sum_sq_dev == 0.0,
            "and nothing outside the merged channels");
}

void TestASinkThatSaysNoIsReported() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(kEverything), 1, 1);
    const bool kept_going = slots.HandOver(0, false, 0.5, [](const SpectralBlock&) { return false; });
    Require(!kept_going, "a sink that says no stops the reduction");

    slots.Fold<true>(0, 0, Row(1, 0, 2.0, 4.0, 2.0, 2.0));
    RequireTotals(Totals(1, 0, 2.0, 4.0, 2.0, 2.0), HandOver(slots).totals.at(0).at(0), "and still leaves identities");
}

void TestTheBlockSaysWhereItIs() {
    StatisticSlots slots;
    slots.Reset(StatisticLayout::Of(kEverything), 2, 3);
    const auto seen = HandOver(slots, false, 40);
    Require(seen.first_channel == 40, "the first channel is the caller's");
    Require(seen.channel_count == 3 && seen.region_count == 2, "the block's shape is the reset's");
    Require(!seen.complete && seen.completeness == 0.5, "and so is how finished it is");
}

}  // namespace

int main() {
    try {
        TestTheLayoutFollowsTheOneOrder();
        TestAResetHoldsNothing();
        TestAFoldAddsTheCountsAndKeepsTheExtrema();
        TestOnlyWhatWasAskedForIsCarried();
        TestPartialsMergeInOrderAtTheirOffset();
        TestDeviationsFoldAsOneSetOfPixels();
        TestDeviationsMergeInOrderAsOneSetOfPixels();
        TestAHandOverPutsTheIdentitiesBack();
        TestASinkThatSaysNoIsReported();
        TestTheBlockSaysWhereItIs();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "statistic slots test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
