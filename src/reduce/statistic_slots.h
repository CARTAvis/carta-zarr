/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_STATISTIC_SLOTS_H_
#define CARTA_ZARR_SRC_REDUCE_STATISTIC_SLOTS_H_

// Where a spectral reduction keeps what it has counted, and the block it hands that over as.
//
// The layout is [region][statistic][channel], with only the statistics that were asked for, in
// kStatisticOrder, and nothing outside this file indexes into it. The block accumulator and each
// task's private partial are the same thing at two extents, so they are one type here.
//
// Not the cube histogram's accumulator, although it counts the same statistics. It keeps them
// in registers for one region and one channel, decides an extremum is untouched by num_pixels being
// zero rather than by it still being infinite, and is padded to a cache line beside its histogram,
// which is a number ADR 0005 measured. Folding it into this would move that number for nothing.

#include "carta-zarr/reduce.h"

#include "reduce/deviations.h"
#include "reduce/tuning.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace carta::zarr::internal {

// One row of one region inside one chunk, accumulated in registers so that the per-pixel loop never
// asks which statistics were requested. The fold into the slots happens once per row.
struct RowTotals {
    std::uint64_t good = 0;
    std::uint64_t bad = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double smallest = std::numeric_limits<double>::infinity();
    double largest = -std::numeric_limits<double>::infinity();
    // The good pixels' mean, as a base and an offset, and their squared deviations from it; left at
    // zero by a loop not asked for them. See reduce/deviations.h.
    double base = 0.0;
    double offset = 0.0;
    double sum_sq_dev = 0.0;
};

// Which statistics a reduction accumulates, and the slot each one has.
class StatisticLayout {
public:
    // The statistics in `requested`, in kStatisticOrder. A request for none has been refused before
    // anything asks for its layout.
    //
    // sum_sq_dev brings num_pixels and sum with it: two sets' deviations are put together from their
    // counts and means, so a fold or a merge without them would have nothing to do it with.
    static StatisticLayout Of(StatisticSet requested) {
        if (requested.Contains(Statistic::sum_sq_dev)) {
            requested |= Statistic::num_pixels | Statistic::sum;
        }
        StatisticLayout layout;
        layout._slot_of.fill(-1);
        for (std::size_t i = 0; i < kStatisticOrder.size(); ++i) {
            if (requested.Contains(kStatisticOrder.at(i))) {
                layout._slot_of.at(i) = static_cast<int>(layout._count);
                layout._statistics.at(layout._count++) = kStatisticOrder.at(i);
            }
        }
        return layout;
    }

    std::size_t count() const noexcept { return _count; }
    // count() of them, in the order a block reports them.
    const Statistic* statistics() const noexcept { return _statistics.data(); }
    // Whether the per-pixel loop has deviations to take.
    bool Deviations() const noexcept { return _slot_of.back() >= 0; }
    // How many values one region keeps at one channel: the statistics, and with sum_sq_dev the base
    // and offset of the mean its merges need, which no block reports. See reduce/deviations.h.
    std::size_t Stored() const noexcept { return _count + (Deviations() ? 2 : 0); }
    // What one channel of a block over `regions` regions occupies.
    std::size_t BytesPerChannel(std::size_t regions) const noexcept { return regions * Stored() * sizeof(double); }

    bool operator==(const StatisticLayout& other) const noexcept { return _slot_of == other._slot_of; }

private:
    friend class StatisticSlots;

    // Only the first count() of these mean anything; the rest are whatever the order left there.
    std::array<Statistic, kStatisticOrder.size()> _statistics{kStatisticOrder};
    std::size_t _count = 0;
    // By position in kStatisticOrder; -1 for a statistic not accumulated.
    std::array<int, kStatisticOrder.size()> _slot_of{};
};

// A block's accumulator, or one task's private partial of it: the same layout at two extents.
//
// The order of use is Reset, then any number of Fold and MergeInOrder, then HandOver -- and after a
// HandOver, more of the same. Nothing here fails at run time; the preconditions are the caller's,
// and a Debug build checks them.
class StatisticSlots {
public:
    // Reshape to `regions` x `channels` and start from nothing: zero for the counts and sums, and
    // the identities of the extrema, which is what lets a fold never ask whether it is the first.
    // Keeps its storage, so a partial reused slab after slab does not allocate.
    void Reset(const StatisticLayout& layout, std::size_t regions, std::uint64_t channels) {
        _layout = layout;
        _regions = regions;
        _channels = static_cast<std::size_t>(channels);
        _region_stride = layout.Stored() * _channels;
        _base_at = layout.count() * _channels;
        _offset_at = _base_at + _channels;
        // Worked out once here rather than from the slot table on every fold: the fold runs once per
        // row, and these are loop invariants for all of them.
        for (std::size_t i = 0; i < kStatisticOrder.size(); ++i) {
            const int slot = layout._slot_of.at(i);
            _offset.at(i) = slot < 0 ? kAbsent : static_cast<std::size_t>(slot) * _channels;
        }

        _values.assign(_regions * _region_stride, 0.0);
        for (std::size_t r = 0; r < _regions; ++r) {
            for (const auto& [i, identity] : kExtrema) {
                if (_offset.at(i) != kAbsent) {
                    double* base = _values.data() + (r * _region_stride) + _offset.at(i);
                    std::fill(base, base + _channels, identity);
                }
            }
        }
    }

    // One row's totals into one region at one channel. The hot path: it runs once for every row of
    // every region in every chunk, so it stays inline and asks only the questions Reset answered.
    //
    // kDeviations says whether the layout has sum_sq_dev, which the caller knows once for every row
    // it folds: asked here, of every row, it would cost a reduction of one-pixel strips 2% whether or
    // not it was asked for the spread.
    template <bool kDeviations>
    void Fold(std::size_t region, std::uint64_t channel, const RowTotals& row) noexcept {
        assert(region < _regions && channel < _channels);
        double* at = _values.data() + (region * _region_stride) + static_cast<std::size_t>(channel);
        // First, while the count is still that of what came before the row.
        assert(kDeviations == (std::get<kSumSqDev>(_offset) != kAbsent));
        if constexpr (kDeviations) {
            double& deviations = at[std::get<kSumSqDev>(_offset)];
            Spread spread{at[std::get<kNumPixels>(_offset)], at[_base_at], at[_offset_at], deviations};
            spread.Merge({static_cast<double>(row.good), row.base, row.offset, row.sum_sq_dev});
            at[_base_at] = spread.base;
            at[_offset_at] = spread.offset;
            deviations = spread.sum_sq_dev;
        }
        if (std::get<kNumPixels>(_offset) != kAbsent) {
            at[std::get<kNumPixels>(_offset)] += static_cast<double>(row.good);
        }
        if (std::get<kNanCount>(_offset) != kAbsent) {
            at[std::get<kNanCount>(_offset)] += static_cast<double>(row.bad);
        }
        if (std::get<kSum>(_offset) != kAbsent) {
            at[std::get<kSum>(_offset)] += row.sum;
        }
        if (std::get<kSumSq>(_offset) != kAbsent) {
            at[std::get<kSumSq>(_offset)] += row.sum_sq;
        }
        if (std::get<kMin>(_offset) != kAbsent) {
            double& current = at[std::get<kMin>(_offset)];
            current = std::min(current, row.smallest);
        }
        if (std::get<kMax>(_offset) != kAbsent) {
            double& current = at[std::get<kMax>(_offset)];
            current = std::max(current, row.largest);
        }
    }

    // Adds partials[0], partials[1], ... into channels [offset, offset + their channels), in that
    // order. The order is the determinism rule -- sums added in another association are different
    // bits -- which is why this takes them all rather than one at a time. Every partial has this
    // block's layout and regions, and fits at `offset`.
    void MergeInOrder(const StatisticSlots* partials, std::size_t count, std::uint64_t offset) {
        for (std::size_t task = 0; task < count; ++task) {
            const StatisticSlots& partial = partials[task];
            assert(partial._layout == _layout && partial._regions == _regions &&
                   offset + partial._channels <= _channels);
            for (std::size_t r = 0; r < _regions; ++r) {
                // The deviations first, as in Fold: they need the counts from before.
                if (std::get<kSumSqDev>(_offset) != kAbsent) {
                    const double* from = partial._values.data() + (r * partial._region_stride);
                    double* to = _values.data() + (r * _region_stride) + static_cast<std::size_t>(offset);
                    const std::size_t from_count = std::get<kNumPixels>(partial._offset);
                    const std::size_t from_deviations = std::get<kSumSqDev>(partial._offset);
                    const std::size_t to_count = std::get<kNumPixels>(_offset);
                    const std::size_t to_deviations = std::get<kSumSqDev>(_offset);
                    for (std::size_t c = 0; c < partial._channels; ++c) {
                        Spread spread{to[to_count + c], to[_base_at + c], to[_offset_at + c], to[to_deviations + c]};
                        spread.Merge({from[from_count + c], from[partial._base_at + c], from[partial._offset_at + c],
                                      from[from_deviations + c]});
                        to[_base_at + c] = spread.base;
                        to[_offset_at + c] = spread.offset;
                        to[to_deviations + c] = spread.sum_sq_dev;
                    }
                }
                for (std::size_t slot = 0; slot < _layout.count(); ++slot) {
                    const double* from =
                        partial._values.data() + (r * partial._region_stride) + (slot * partial._channels);
                    double* to =
                        _values.data() + (r * _region_stride) + (slot * _channels) + static_cast<std::size_t>(offset);
                    const Statistic statistic = _layout._statistics.at(slot);
                    if (statistic == Statistic::sum_sq_dev) {
                        continue;
                    }
                    if (statistic == Statistic::min) {
                        for (std::size_t c = 0; c < partial._channels; ++c) {
                            to[c] = std::min(to[c], from[c]);
                        }
                    } else if (statistic == Statistic::max) {
                        for (std::size_t c = 0; c < partial._channels; ++c) {
                            to[c] = std::max(to[c], from[c]);
                        }
                    } else {
                        for (std::size_t c = 0; c < partial._channels; ++c) {
                            to[c] += from[c];
                        }
                    }
                }
            }
        }
    }

    // Hands the block to `sink` as a SpectralBlock and returns what the sink returned.
    //
    // An extremum nothing contributed to is still its identity, which is the one value it must not
    // be reported as, so it is reported as NaN. Finding those needs no bookkeeping: a finite pixel
    // can never leave an infinity behind, so only the untouched entries are still infinite -- and by
    // the same argument the NaN that replaced one is the only NaN there, so the identity can be put
    // back afterwards. It is put back every time, finished block or not: the values are only the
    // sink's for the length of the call, so nobody can tell, and the accumulator a hand-over leaves
    // is then always the one it found.
    template <typename Sink>
    bool HandOver(std::uint64_t first_channel, bool complete, double completeness, const Sink& sink) {
        Settle(true);
        SpectralBlock block;
        block.first_channel = first_channel;
        block.channel_count = _channels;
        block.region_count = _regions;
        block._values = _values.data();
        block._region_stride = _region_stride;
        block._statistics = _layout.statistics();
        block._statistic_count = _layout.count();
        block.complete = complete;
        block.completeness = completeness;
        const bool keep_going = sink(block);
        Settle(false);
        return keep_going;
    }

    std::size_t regions() const noexcept { return _regions; }
    std::uint64_t channels() const noexcept { return _channels; }

private:
    static constexpr std::size_t kAbsent = std::numeric_limits<std::size_t>::max();
    // Positions in kStatisticOrder, which the offsets are indexed by.
    static constexpr std::size_t kNumPixels = 0;
    static constexpr std::size_t kNanCount = 1;
    static constexpr std::size_t kSum = 2;
    static constexpr std::size_t kSumSq = 3;
    static constexpr std::size_t kMin = 4;
    static constexpr std::size_t kMax = 5;
    static constexpr std::size_t kSumSqDev = 6;
    static_assert(std::get<kNumPixels>(kStatisticOrder) == Statistic::num_pixels &&
                      std::get<kNanCount>(kStatisticOrder) == Statistic::nan_count &&
                      std::get<kSum>(kStatisticOrder) == Statistic::sum &&
                      std::get<kSumSq>(kStatisticOrder) == Statistic::sum_sq &&
                      std::get<kMin>(kStatisticOrder) == Statistic::min &&
                      std::get<kMax>(kStatisticOrder) == Statistic::max &&
                      std::get<kSumSqDev>(kStatisticOrder) == Statistic::sum_sq_dev,
                  "the positions above are kStatisticOrder's");

    struct Extremum {
        std::size_t position;
        double identity;
    };
    static constexpr std::array<Extremum, 2> kExtrema{Extremum{kMin, std::numeric_limits<double>::infinity()},
                                                      Extremum{kMax, -std::numeric_limits<double>::infinity()}};

    // Untouched extrema to NaN when `report`, and back to their identities when not.
    void Settle(bool report) {
        for (std::size_t r = 0; r < _regions; ++r) {
            for (const auto& [i, identity] : kExtrema) {
                if (_offset.at(i) == kAbsent) {
                    continue;
                }
                double* base = _values.data() + (r * _region_stride) + _offset.at(i);
                for (std::size_t c = 0; c < _channels; ++c) {
                    if (report) {
                        if (std::isinf(base[c])) {
                            base[c] = std::numeric_limits<double>::quiet_NaN();
                        }
                    } else if (std::isnan(base[c])) {
                        base[c] = identity;
                    }
                }
            }
        }
    }

    StatisticLayout _layout;
    std::size_t _regions = 0;
    std::size_t _channels = 0;
    std::size_t _region_stride = 0;
    std::array<std::size_t, kStatisticOrder.size()> _offset{};
    // Where in a region the mean behind sum_sq_dev is kept, after every statistic a block reports.
    std::size_t _base_at = 0;
    std::size_t _offset_at = 0;
    std::vector<double> _values;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_STATISTIC_SLOTS_H_
