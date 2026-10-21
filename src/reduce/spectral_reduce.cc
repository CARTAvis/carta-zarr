/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "spectral_reduce.h"

#include "axis_map.h"
#include "chunk_blocks.h"
#include "reduce/deviations.h"
#include "reduce/occupancy.h"
#include "reduce/pass.h"
#include "reduce/plane_selection.h"
#include "reduce/statistic_slots.h"
#include "reduce/tuning.h"
#include "zarr/pixel_selection.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace carta::zarr::internal {
namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

// The per-pixel loop, written so that a compiler can vectorize it.
//
// The first two template parameters are loop invariants rather than runtime tests, because each one
// alone is enough to stop vectorization. The u stride is 1 for every image, since the walk runs along
// the axis the store varies fastest, but the compiler cannot know that and pays a multiply per pixel
// for the possibility; and a mask costs a branch per pixel. Only a raster too fragmented for runs
// still takes that branch: every other masked region arrives here as runs, which Occupancy::Of
// makes.
//
// The finiteness test is branchless for the same reason. A non-finite value contributes zero to the
// sums and its own identity to the extrema, which is exactly what excluding it means, so there is
// nothing an if would do that arithmetic does not.
//
// kDeviations is the third: the distances from a shift that Statistic::sum_sq_dev is made from cost
// a subtraction and a multiply-add per pixel, which a reduction not asked for it does not pay. A
// non-finite pixel's distance is zero, as is its value. See reduce/deviations.h.
template <bool kUnitStride, bool kMasked, bool kDeviations>
void AccumulateRow(const float* row, std::uint64_t stride, std::uint64_t count, const std::uint8_t* selected,
                   std::uint64_t mask_stride, RowTotals& totals) {
    std::uint64_t good = 0;
    std::uint64_t selected_count = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double smallest = kInfinity;
    double largest = -kInfinity;
    const double shift = kDeviations ? DeviationShift<kMasked>(row, stride, count, selected, mask_stride) : 0.0;
    double distance_sum = 0.0;
    double distance_sum_sq = 0.0;

    for (std::uint64_t i = 0; i < count; ++i) {
        if (kMasked && selected[i * mask_stride] == 0) {
            continue;
        }
        ++selected_count;
        const float value = kUnitStride ? row[i] : row[i * stride];
        const bool finite = std::isfinite(value);
        const double clean = finite ? static_cast<double>(value) : 0.0;
        good += static_cast<std::uint64_t>(finite);
        sum += clean;
        sum_sq += clean * clean;
        smallest = std::min(smallest, finite ? clean : kInfinity);
        largest = std::max(largest, finite ? clean : -kInfinity);
        if (kDeviations) {
            const double distance = finite ? clean - shift : 0.0;
            distance_sum += distance;
            distance_sum_sq += distance * distance;
        }
    }

    if (kDeviations) {
        // Before the count below takes this span in: the row so far is what came before it.
        Spread row_so_far{static_cast<double>(totals.good), totals.base, totals.offset, totals.sum_sq_dev};
        row_so_far.Merge(Spread::OfSpan(static_cast<double>(good), shift, distance_sum, distance_sum_sq));
        totals.base = row_so_far.base;
        totals.offset = row_so_far.offset;
        totals.sum_sq_dev = row_so_far.sum_sq_dev;
    }

    // Added rather than assigned: a row described by runs is several spans, and they share one set
    // of totals. A row that is one span starts from the same zeros and reaches the same answer.
    totals.good += good;
    totals.bad += selected_count - good;
    totals.sum += sum;
    totals.sum_sq += sum_sq;
    totals.smallest = std::min(totals.smallest, smallest);
    totals.largest = std::max(totals.largest, largest);
}

// One span into a row's totals, with the loop its stride, its mask and the request call for. Every
// combination is its own instantiation, so that no loop carries the tests the others need.
template <bool kDeviations>
void AccumulateSpan(const float* pixels, std::uint64_t stride, std::uint64_t count, const std::uint8_t* selected,
                    std::uint64_t mask_step, RowTotals& totals) {
    if (stride == 1) {
        if (selected == nullptr) {
            AccumulateRow<true, false, kDeviations>(pixels, 1, count, nullptr, 1, totals);
        } else {
            AccumulateRow<true, true, kDeviations>(pixels, 1, count, selected, mask_step, totals);
        }
    } else if (selected == nullptr) {
        AccumulateRow<false, false, kDeviations>(pixels, stride, count, nullptr, 1, totals);
    } else {
        AccumulateRow<false, true, kDeviations>(pixels, stride, count, selected, mask_step, totals);
    }
}

}  // namespace

Result<void> ReduceSpectral(const ReducibleImage& image, const SpectralReduceRequest& request, const SpectralSink& sink,
                            const ReadOptions& options) {
    const auto& descriptor = image.descriptor();
    const auto& source = image.source();
    const auto& map = image.map();
    const auto& node = descriptor.id;
    if (!sink) {
        return Error{ErrorCode::invalid_argument, "A spectral reduction needs a sink", node};
    }

    // Occupancy::Of checks the regions as it places them; the plan checks everything else about the
    // request.
    if (request.statistics.empty()) {
        return Error{ErrorCode::invalid_argument, "A spectral reduction needs at least one statistic", node};
    }

    // Checked here as well as inside each slab request, so that a bad range is one error naming the
    // axis rather than a partial reduction that fails on some later slab.
    const auto planned = image.Plan(request.planes, 1, options);
    if (!planned) {
        return planned.error();
    }
    const auto& plan = planned.value();

    // The caller's regions, placed on the walk's axes and indexed by the chunks they occupy. Which
    // spatial axis the store varies fastest decides where a region lands and which way its runs lie,
    // so it is taken from the plan, which applied that rule when it picked axis_u, rather than worked
    // out again from the geometry.
    const PlaneExtent plane{descriptor.axes.at(map.x).length, descriptor.axes.at(map.y).length};
    auto occupancy_result = Occupancy::Of(request.regions, plane, plan.chunk_u, plan.chunk_v,
                                          plan.SwapsSpatial() ? AxisRole::spatial_y : AxisRole::spatial_x, node);
    if (!occupancy_result) {
        return occupancy_result.error();
    }
    const auto& occupancy = occupancy_result.value();
    const auto& regions = occupancy.regions();
    // What the reduction reads, cut once against the plan's budget. They do not depend on which
    // channels a block covers, so every block walks the same ones.
    const auto footprints = occupancy.Footprints(plan.ChunksPerRead());

    // The requested statistics, in the one order a block reports them.
    const auto layout = StatisticLayout::Of(request.statistics);
    const bool deviations = layout.Deviations();

    // The budget counts the chunk data a slab decodes, not the pixels it keeps. A one-pixel region
    // asks for almost nothing and still decodes every chunk it touches, so sizing by the region's own
    // area would let a cursor-sized request pull through an unbounded amount. The plan holds this
    // rule.

    // The pass walks the footprints, and counts a block's progress in the chunks they occupy
    // together. Not the plan's layer, which is the whole plane: a reduction spends its emit budget
    // against the region set it was given.
    //
    // The pass works out how many channels one emitted block may hold. The hint is the caller's, the
    // budgets are the library's and the chunk alignment is the walk's; the smallest wins, and the
    // block reports what it used.
    //
    // Without a hint, a block costs one budget of decoded bytes, as a piece of Read does, and spends
    // whatever the spatial walk left over. A small region leaves almost all of it, so the block spans
    // many chunks along the spectrum. A region covering the image spends the budget spatially, what is
    // left allows one chunk along the spectrum, and the block is the single chunk layer the walk is
    // already reading: results are handed over a layer at a time, as soon as they are finished.
    //
    // Emitting only at the end would be silent for the whole reduction. One layer of a 7763x4742
    // image is 160 MiB and 70 ms; a thousand channels is a minute of work with no partial answer and
    // no chance to cancel.
    auto pass =
        PassOverFootprints(source, plan, options, footprints, "The spectral reduction was cancelled by its sink");

    StatisticSlots accumulator;
    // One private accumulator per task, reused across slabs; see the dispatch below. Called
    // "partials" as in the plane histogram, not "sinks": a sink in this library is where a finished
    // block goes, such as the caller's SpectralSink in this same function.
    std::vector<StatisticSlots> partials;

    const auto reset_block = [&](std::uint64_t length) { accumulator.Reset(layout, request.regions.size, length); };

    // One read's worth of pixels, accumulated into the block's own totals.
    //
    // Named rather than written into the call below, for the reason plane_histogram.cc gives for
    // bin_slab, and because it holds another lambda.
    //
    // The pass hands it the footprint each slab was read over, unpacked into locals below. It stays a
    // lambda passed as a template parameter, never a std::function: the per-pixel loop inlines
    // through it. ADR 0005.
    const auto accumulate_slab = [&](const OccupiedFootprint& footprint, const Slab& slab) {
        const std::uint64_t chunk_cv_begin = footprint.chunk_cv_begin;
        const std::uint64_t chunk_cv_end = footprint.chunk_cv_end;
        const std::uint64_t chunk_cu_begin = footprint.chunk_cu_begin;
        const std::uint64_t chunk_cu_end = footprint.chunk_cu_end;
        const std::uint64_t u_begin = footprint.slab.u_start;
        const std::uint64_t u_end = footprint.slab.u_start + footprint.slab.u_count;
        const std::uint64_t v_begin = footprint.slab.v_start;
        const std::uint64_t v_end = footprint.slab.v_start + footprint.slab.v_count;

        // Into locals for the accumulation below.
        const float* const slab_pixels = slab.pixels;
        const std::uint64_t stride_u = slab.stride_u;
        const std::uint64_t stride_v = slab.stride_v;
        const std::uint64_t stride_z = slab.stride_z;
        const std::uint64_t slab_length = slab.channel_count;

        const std::uint64_t cv_span = chunk_cv_end - chunk_cv_begin;
        const std::uint64_t cu_span = chunk_cu_end - chunk_cu_begin;
        const std::uint64_t cells = std::max<std::uint64_t>(1, cv_span * cu_span);
        const std::uint64_t units = slab_length * cells;
        const std::size_t partial_bytes = std::max<std::size_t>(
            1, layout.BytesPerChannel(request.regions.size) * static_cast<std::size_t>(slab_length));

        // One (channel, chunk cell) unit of the accumulation, into a private partial the length of
        // this slab. Private because two units of the same channel can touch the same region (a
        // region wider than a chunk spans several cells), and would otherwise add to one double.
        //
        // `deviations` arrives as a type, so that neither of the two instantiations tests it for
        // every row and span below: see AccumulateSpan and Fold.
        const auto accumulate_unit = [&](auto deviations, std::uint64_t channel, std::uint64_t cell_cv,
                                         std::uint64_t cell_cu, StatisticSlots& partial) {
            constexpr bool kDeviations = decltype(deviations)::value;
            const float* plane = slab_pixels + (channel * stride_z);

            const std::uint64_t chunk_cv = cell_cv;
            const std::uint64_t cell_v0 = std::max(v_begin, chunk_cv * plan.chunk_v);
            const std::uint64_t cell_v1 = std::min(v_end, (chunk_cv + 1) * plan.chunk_v);
            if (cell_v0 >= cell_v1) {
                return;
            }
            const std::uint64_t chunk_cu = cell_cu;
            const std::uint64_t cell_u0 = std::max(u_begin, chunk_cu * plan.chunk_u);
            const std::uint64_t cell_u1 = std::min(u_end, (chunk_cu + 1) * plan.chunk_u);
            if (cell_u0 >= cell_u1) {
                return;
            }

            for (const auto index : occupancy.RegionsTouching(chunk_cu, chunk_cv)) {
                const std::size_t r = index;
                const auto& region = regions.at(r);
                const auto rows = region.RowsWithin(cell_u0, cell_u1, cell_v0, cell_v1);

                for (std::uint64_t y = rows.first; y < rows.last; ++y) {
                    // One row's totals, however many spans it has, folded once, so the sums are taken
                    // in row order.
                    RowTotals totals;
                    const float* row = plane + ((y - v_begin) * stride_v);
                    region.ForEachSpan(y, cell_u0, cell_u1,
                                       [&](std::uint64_t first, std::uint64_t last, const std::uint8_t* selected,
                                           std::uint64_t mask_step) {
                                           // A span whose every pixel is selected -- a run, or a region
                                           // that is its whole box -- takes the loop with no test per pixel.
                                           const float* pixels = row + ((first - u_begin) * stride_u);
                                           AccumulateSpan<kDeviations>(pixels, stride_u, last - first, selected,
                                                                       mask_step, totals);
                                       });
                    partial.Fold<kDeviations>(r, channel, totals);
                }
            }
        };

        // The unit order is channel, then chunk row, then chunk column, so a single task sums in the
        // order a plain nested loop would. Several tasks do not: each sums its own contiguous run of
        // units, and the partials are added in task order, which is deterministic but a different
        // association. The statistics are doubles over as many as 5x10^11 values, and partial sums
        // are if anything more accurate; the tests compare to 1e-9 relative for this reason. The
        // counts and extrema are unaffected, since integers and min/max do not depend on order. There
        // is one partial per task, so the split also costs an allocation and a memset of this size per
        // slab, which is why the split is made here rather than once: a partial is as long as the
        // slab. See kSpectralPartialBudgetBytes.
        const auto split = image.Split(kSpectralPartialBudgetBytes, partial_bytes);
        const std::size_t tasks = split.Tasks(plan.chunk_u * plan.chunk_v, units);

        if (partials.size() < tasks) {
            partials.resize(tasks);
        }
        for (std::size_t task = 0; task < tasks; ++task) {
            partials.at(task).Reset(layout, request.regions.size, slab_length);
        }

        // Contiguous runs of units, cut the way the histograms cut their rows. Dealing them out
        // round-robin measured 5-13% slower at ten threads on a 7763x4742 plane, and never faster.
        split.Run(tasks, units, [&](std::size_t task, std::uint64_t first, std::uint64_t last) {
            StatisticSlots& partial = partials.at(task);
            for (std::uint64_t unit = first; unit < last; ++unit) {
                const std::uint64_t channel = unit / cells;
                const std::uint64_t cell = unit % cells;
                const std::uint64_t cell_cv = chunk_cv_begin + (cell / cu_span);
                const std::uint64_t cell_cu = chunk_cu_begin + (cell % cu_span);
                if (deviations) {
                    accumulate_unit(std::true_type{}, channel, cell_cv, cell_cu, partial);
                } else {
                    accumulate_unit(std::false_type{}, channel, cell_cv, cell_cu, partial);
                }
            }
        });

        accumulator.MergeInOrder(partials.data(), tasks, slab.first_channel.index);
    };

    // From the typed index into the public block's plain one, here and nowhere else.
    const auto hand_over = [&](SelectionChannel first_channel, [[maybe_unused]] std::uint64_t length, bool complete,
                               double completeness) {
        assert(length == accumulator.channels());
        return accumulator.HandOver(first_channel.index, complete, completeness, sink);
    };

    return pass.InBlocks(layout.BytesPerChannel(request.regions.size), request.emit_every_channels, reset_block,
                         accumulate_slab, hand_over);
}

}  // namespace carta::zarr::internal
