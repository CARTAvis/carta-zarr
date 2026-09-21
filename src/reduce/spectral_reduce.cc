/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "spectral_reduce.h"

#include "chunk_blocks.h"
#include "axis_map.h"
#include "reduce/tuning.h"
#include "reduce/block_emit.h"
#include "reduce/occupancy.h"
#include "reduce/pass.h"
#include "reduce/plane_selection.h"
#include "zarr/pixel_selection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace carta::zarr::internal {
namespace {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

Result<void> ValidateRequest(const ImageDescriptor& descriptor, const AxisMap& axes,
                             const SpectralReduceRequest& request) {
    const auto& node = descriptor.id;
    if (request.region_count == 0 || request.regions == nullptr) {
        return Error{ErrorCode::invalid_argument, "A spectral reduction needs at least one region", node};
    }
    if (request.region_count > kMaxSpectralRegions) {
        return Error{ErrorCode::invalid_argument,
                     "A spectral reduction accepts at most " + std::to_string(kMaxSpectralRegions) +
                         " regions, not " + std::to_string(request.region_count),
                     node};
    }
    if (request.statistics == 0) {
        return Error{ErrorCode::invalid_argument, "A spectral reduction needs at least one statistic", node};
    }

    const auto width = descriptor.axes.at(axes.x).length;
    const auto height = descriptor.axes.at(axes.y).length;
    for (std::size_t i = 0; i < request.region_count; ++i) {
        const auto& region = request.regions[i];
        if (region.width == 0 || region.height == 0) {
            return Error{ErrorCode::invalid_argument, "Region " + std::to_string(i) + " is empty", node};
        }
        if (region.x_start >= width || region.width > width - region.x_start || region.y_start >= height ||
            region.height > height - region.y_start) {
            return Error{ErrorCode::invalid_argument,
                         "Region " + std::to_string(i) + " falls outside the image", node};
        }
        if ((region.row_runs == nullptr) != (region.row_run_offsets == nullptr)) {
            return Error{ErrorCode::invalid_argument,
                         "Region " + std::to_string(i) + " supplies one run array without the other", node};
        }
    }

    return {};
}

// One row of one region inside one chunk, accumulated in registers so that the per-pixel loop never
// asks which statistics were requested. The merge into the output happens once per row.
struct RowTotals {
    std::uint64_t good = 0;
    std::uint64_t bad = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double smallest = kInfinity;
    double largest = -kInfinity;
};

// The per-pixel loop, written so that a compiler can vectorise it.
//
// Both template parameters are loop invariants that used to be runtime tests, and each one on its
// own was enough to stop vectorisation: the x stride is 1 for every image whose fastest logical
// axis is x, which is every XRADIO image, but the compiler cannot know that and pays a multiply per
// pixel for the possibility; and a mask that most regions do not have cost a branch per pixel.
//
// The finiteness test is branchless for the same reason. A non-finite value contributes zero to the
// sums and its own identity to the extrema, which is exactly what excluding it means, so there is
// nothing an if would do that arithmetic does not.
template <bool kUnitStride, bool kMasked>
void AccumulateRow(const float* row, std::uint64_t stride, std::uint64_t count, const std::uint8_t* selected,
                   std::uint64_t mask_stride, RowTotals& totals) {
    std::uint64_t good = 0;
    std::uint64_t selected_count = 0;
    double sum = 0.0;
    double sum_sq = 0.0;
    double smallest = kInfinity;
    double largest = -kInfinity;

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

// The chunk cells and pixel bounds of one spatial footprint.
//
// Everything else the accumulation needs belongs to the reduction, so this is the whole of what
// changes from one footprint to the next -- and the whole reason the accumulation had to be written
// inside the loop that cuts them.
struct FootprintBounds {
    std::uint64_t chunk_cv_begin = 0;
    std::uint64_t chunk_cv_end = 0;
    std::uint64_t chunk_cu_begin = 0;
    std::uint64_t chunk_cu_end = 0;
    std::uint64_t u_begin = 0;
    std::uint64_t u_end = 0;
    std::uint64_t v_begin = 0;
    std::uint64_t v_end = 0;
};

}  // namespace

Result<void> ReduceSpectral(const ReadableImage& image, const SpectralReduceRequest& request,
                            const SpectralSink& sink, const ReadOptions& options) {
    const auto& descriptor = image.descriptor();
    const auto& geometry = image.geometry();
    const auto& source = image.source();
    const auto& map = image.map();
    auto& workers = image.workers();
    const auto& node = descriptor.id;
    if (!sink) {
        return Error{ErrorCode::invalid_argument, "A spectral reduction needs a sink", node};
    }

    if (auto valid = ValidateRequest(descriptor, map, request); !valid) {
        return valid.error();
    }

    // Checked here as well as inside each slab request, so that a bad range is one error naming the
    // axis rather than a partial reduction that fails on some later slab.
    const auto checked = CheckedPlanes::Of(descriptor, map, request.planes);
    if (!checked) {
        return checked.error();
    }
    const auto& planes = checked.value();

    const auto plan = PlanPass(descriptor, geometry, map, planes, 1, options);

    // The caller's regions, placed on the walk's axes and indexed by the chunks they occupy. Which
    // spatial axis the store varies fastest decides both -- where a region lands, and whether its
    // runs are the kind worth taking -- so the plan is asked for it once and the occupancy is told.
    // Asked of the plan rather than of the geometry: the plan applied that rule when it picked
    // axis_u, and working it out again here would be a second place for it to be got wrong.
    auto occupancy_result =
        Occupancy::Of(request.regions, request.region_count, plan.chunk_u, plan.chunk_v,
                      plan.SwapsSpatial() ? AxisRole::spatial_y : AxisRole::spatial_x, node);
    if (!occupancy_result) {
        return occupancy_result.error();
    }
    const auto& occupancy = occupancy_result.value();
    const auto& regions = occupancy.regions();
    const auto& runs_per_row = occupancy.runs_per_row();

    // The requested statistics, in the one order a block reports them.
    std::vector<Statistic> statistics;
    std::array<int, kStatisticOrder.size()> slot_of{};
    slot_of.fill(-1);
    for (std::size_t i = 0; i < kStatisticOrder.size(); ++i) {
        if (Contains(request.statistics, kStatisticOrder.at(i))) {
            slot_of.at(i) = static_cast<int>(statistics.size());
            statistics.push_back(kStatisticOrder.at(i));
        }
    }
    const std::size_t statistic_count = statistics.size();
    const int slot_num_pixels = slot_of.at(0);
    const int slot_nan_count = slot_of.at(1);
    const int slot_sum = slot_of.at(2);
    const int slot_sum_sq = slot_of.at(3);
    const int slot_min = slot_of.at(4);
    const int slot_max = slot_of.at(5);

    // The budget is over the chunk data a slab decodes, not the pixels it keeps. A one-pixel region
    // asks for almost nothing and still decodes an entire chunk per chunk it touches, so sizing by
    // the region's own area would let a cursor-sized request pull an unbounded amount through. That
    // is the plan's rule as much as it is this one's, which is why the plan holds it.

    // The chunks one spectral layer of the whole region set occupies. Not the plan's layer, which
    // is the whole plane: a reduction spends its emit budget against the region set it was given.
    const std::uint64_t layer_chunks = occupancy.LayerChunks();

    // How many channels one emitted block may hold. The hint is the caller's, the budgets are the
    // library's and the chunk alignment is the walk's; the smallest wins, and the block reports
    // what it used.
    //
    // Without a hint a block costs one budget of decoded bytes, the same invariant a piece of Read
    // carries, which is why it is free: the block spends whatever the spatial walk left over. A
    // small region leaves almost all of it, so the block spans many chunks along the spectrum. A
    // region covering the image spends the budget spatially, `spectral_chunks` below collapses to
    // one, and the block becomes the single chunk layer the walk is already reading -- the results
    // are handed over a layer at a time because that is when they are finished, not sooner.
    //
    // Emitting only at the end instead, which is what a zero used to mean, is silent for as long as
    // the whole reduction takes. One layer of a 7763x4742 image is 160 MiB and 70 ms; a thousand
    // channels of it is a minute of work with no partial answer and nowhere to cancel.
    const std::size_t bytes_per_channel = request.region_count * statistic_count * sizeof(double);

    const BlockEmitter emitter(plan, layer_chunks, bytes_per_channel, request.emit_every_channels,
                               "The spectral reduction was cancelled by its sink");

    std::vector<double> accumulator;
    // One private accumulator per task, reused across slabs. See the dispatch below. "Partials" as
    // in the plane histogram, and deliberately not "sinks": a sink in this library is where a
    // finished block goes, and this one is named in the same function as the caller's SpectralSink.
    std::vector<double> partials;
    SlabWalk walk(source, plan, options);
    std::vector<ColumnRun> segments;
    // The accumulator's current shape, which is the length of the block being filled. Set by the
    // reset below and read by everything that indexes into it, including the strides a block reports.
    std::size_t statistic_stride = 0;
    std::size_t region_stride = 0;

    // An extremum nothing contributed to is still its identity, which is the one value it must
    // not be reported as. Finding those needs no extra bookkeeping: a finite pixel can never
    // leave an infinity behind, so only the untouched entries are still infinite -- and by the
    // same argument the NaN that replaced one is the only NaN there, so the identity can be put
    // back when the walk is not finished with it.
    const auto settle_extrema = [&](bool report, std::size_t length) {
        for (std::size_t r = 0; r < request.region_count; ++r) {
            for (const int slot : {slot_min, slot_max}) {
                if (slot < 0) {
                    continue;
                }
                const double identity = slot == slot_min ? kInfinity : -kInfinity;
                auto* base = accumulator.data() + (r * region_stride) +
                             (static_cast<std::size_t>(slot) * statistic_stride);
                for (std::size_t c = 0; c < length; ++c) {
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
    };

    const auto reset_block = [&](std::uint64_t length) {
        statistic_stride = static_cast<std::size_t>(length);
        region_stride = statistic_count * statistic_stride;

        accumulator.assign(request.region_count * region_stride, 0.0);
        for (std::size_t r = 0; r < request.region_count; ++r) {
            if (slot_min >= 0) {
                auto* base = accumulator.data() + (r * region_stride) +
                             (static_cast<std::size_t>(slot_min) * statistic_stride);
                std::fill(base, base + statistic_stride, kInfinity);
            }
            if (slot_max >= 0) {
                auto* base = accumulator.data() + (r * region_stride) +
                             (static_cast<std::size_t>(slot_max) * statistic_stride);
                std::fill(base, base + statistic_stride, -kInfinity);
            }
        }
    };

    // One read's worth of pixels, accumulated into the block's own totals.
    //
    // Named rather than written into the call below, for the reason plane_histogram.cc gives for
    // bin_slab: a visitor nested inside the walk inside the emitter is three lambdas deep before the
    // first loop, and this one carries a fourth inside it.
    //
    // The footprint's bounds arrive as an argument and are unpacked into the names the body already
    // used, so that the accumulation is the same text it has been since it stopped doing its own
    // reading. It stays a lambda passed as a template parameter, never a std::function: the
    // per-pixel loop inlines through it. ADR 0005.
    const auto accumulate_slab = [&](const FootprintBounds& footprint_bounds, const Slab& slab) {
        const std::uint64_t chunk_cv_begin = footprint_bounds.chunk_cv_begin;
        const std::uint64_t chunk_cv_end = footprint_bounds.chunk_cv_end;
        const std::uint64_t chunk_cu_begin = footprint_bounds.chunk_cu_begin;
        const std::uint64_t chunk_cu_end = footprint_bounds.chunk_cu_end;
        const std::uint64_t u_begin = footprint_bounds.u_begin;
        const std::uint64_t u_end = footprint_bounds.u_end;
        const std::uint64_t v_begin = footprint_bounds.v_begin;
        const std::uint64_t v_end = footprint_bounds.v_end;

        // Into locals so that the accumulation below is the text it was when this
        // function did its own reading.
        const float* const slab_pixels = slab.pixels;
        const std::uint64_t stride_u = slab.stride_u;
        const std::uint64_t stride_v = slab.stride_v;
        const std::uint64_t stride_z = slab.stride_z;
        const std::uint64_t slab_length = slab.channel_count;

        const std::uint64_t cv_span = chunk_cv_end - chunk_cv_begin;
        const std::uint64_t cu_span = chunk_cu_end - chunk_cu_begin;
        const std::uint64_t cells = std::max<std::uint64_t>(1, cv_span * cu_span);
        const std::uint64_t units = slab_length * cells;
        const std::size_t partial_region_stride =
            statistic_count * static_cast<std::size_t>(slab_length);
        const std::size_t partial_stride =
            std::max<std::size_t>(1, request.region_count * partial_region_stride);

        // One (channel, chunk cell) unit of the accumulation, into a private partial
        // laid out [region][statistic][channel within this slab]. Private because two
        // units of the same channel can touch the same region -- a region wider than a
        // chunk spans several cells -- so they would otherwise be adding to one double.
        const auto accumulate_unit = [&](std::uint64_t channel, std::uint64_t cell_cv,
                                         std::uint64_t cell_cu, double* partial) {
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

                const std::uint64_t ru0 = std::max(cell_u0, region.u_start);
                const std::uint64_t ru1 = std::min(cell_u1, region.u_start + region.u_size);
                const std::uint64_t rv0 = std::max(cell_v0, region.v_start);
                const std::uint64_t rv1 = std::min(cell_v1, region.v_start + region.v_size);
                if (ru0 >= ru1 || rv0 >= rv1) {
                    continue;
                }

                double* out = partial + (r * partial_region_stride);
                for (std::uint64_t y = rv0; y < rv1; ++y) {
                    RowTotals totals;
                    if (region.runs != nullptr) {
                        // Every pixel of a run is selected, so each one goes
                        // through the same loop an unmasked region uses and the
                        // per-pixel test never happens.
                        const auto r = static_cast<std::size_t>(y - region.v_start);
                        for (auto k = region.run_offsets[r];
                             k < region.run_offsets[r + 1]; ++k) {
                            const std::uint64_t run_begin =
                                region.u_start + region.runs[2 * k];
                            if (run_begin >= ru1) {
                                break;  // runs ascend, so the rest are past the cell
                            }
                            const std::uint64_t run_end =
                                region.u_start + region.runs[(2 * k) + 1];
                            const std::uint64_t first = std::max(ru0, run_begin);
                            const std::uint64_t last = std::min(ru1, run_end);
                            if (first >= last) {
                                continue;
                            }
                            const float* span = plane + ((y - v_begin) * stride_v) +
                                                ((first - u_begin) * stride_u);
                            if (stride_u == 1) {
                                AccumulateRow<true, false>(span, 1, last - first, nullptr, 1,
                                                           totals);
                            } else {
                                AccumulateRow<false, false>(span, stride_u, last - first,
                                                            nullptr, 1, totals);
                            }
                        }
                    } else {
                        const float* pixel_row =
                            plane + ((y - v_begin) * stride_v) + ((ru0 - u_begin) * stride_u);
                        const std::uint8_t* selected =
                            region.mask == nullptr
                                ? nullptr
                                : region.mask +
                                      ((y - region.v_start) * region.mask_v_stride) +
                                      ((ru0 - region.u_start) * region.mask_u_stride);
                        const std::uint64_t run = ru1 - ru0;
                        const std::uint64_t mask_step = region.mask_u_stride;
                        if (stride_u == 1) {
                            if (selected == nullptr) {
                                AccumulateRow<true, false>(pixel_row, 1, run, nullptr, 1, totals);
                            } else {
                                AccumulateRow<true, true>(pixel_row, 1, run, selected, mask_step,
                                                          totals);
                            }
                        } else if (selected == nullptr) {
                            AccumulateRow<false, false>(pixel_row, stride_u, run, nullptr, 1,
                                                        totals);
                        } else {
                            AccumulateRow<false, true>(pixel_row, stride_u, run, selected,
                                                       mask_step, totals);
                        }
                    }

                    if (slot_num_pixels >= 0) {
                        out[(static_cast<std::size_t>(slot_num_pixels) * slab_length) + channel] += static_cast<double>(totals.good);
                    }
                    if (slot_nan_count >= 0) {
                        out[(static_cast<std::size_t>(slot_nan_count) * slab_length) + channel] += static_cast<double>(totals.bad);
                    }
                    if (slot_sum >= 0) {
                        out[(static_cast<std::size_t>(slot_sum) * slab_length) + channel] += totals.sum;
                    }
                    if (slot_sum_sq >= 0) {
                        out[(static_cast<std::size_t>(slot_sum_sq) * slab_length) + channel] += totals.sum_sq;
                    }
                    if (slot_min >= 0) {
                        double& current =
                            out[(static_cast<std::size_t>(slot_min) * slab_length) + channel];
                        current = std::min(current, totals.smallest);
                    }
                    if (slot_max >= 0) {
                        double& current =
                            out[(static_cast<std::size_t>(slot_max) * slab_length) + channel];
                        current = std::max(current, totals.largest);
                    }
                }
            }
        };

        // The unit order is the order the nested loops used to run in -- channel, then
        // chunk row, then chunk column -- so a single task reproduces the old sums bit
        // for bit. Several tasks do not: each one sums its own units and the partials
        // are added in task order, which is deterministic but a different association.
        // The statistics are doubles over as many as 5x10^11 values, and partial sums
        // are if anything the more accurate arrangement; the tests compare to 1e-9
        // relative for exactly this reason, and the counts and extrema are unaffected
        // because integers and min/max do not care what order they arrive in.
        // One per task, so the split is also an allocation and a memset of this size
        // once per slab. Capped so that a reduction over thousands of regions does not
        // spend more on the split than on the pixels.
        constexpr std::size_t kSpectralPartialBudgetBytes = 16U << 20U;
        const std::size_t tasks_by_memory =
            std::max<std::size_t>(1, kSpectralPartialBudgetBytes / (partial_stride * sizeof(double)));
        const std::size_t max_tasks = std::min(workers.size(), tasks_by_memory);
        const std::size_t tasks =
            PlanRowTasks(plan.chunk_u * plan.chunk_v, units, max_tasks, kLeastPixelsPerTask);

        partials.assign(tasks * partial_stride, 0.0);
        for (std::size_t task = 0; task < tasks; ++task) {
            for (std::size_t r = 0; r < request.region_count; ++r) {
                double* base =
                    partials.data() + (task * partial_stride) + (r * partial_region_stride);
                if (slot_min >= 0) {
                    auto* from = base + (static_cast<std::size_t>(slot_min) * slab_length);
                    std::fill(from, from + slab_length, kInfinity);
                }
                if (slot_max >= 0) {
                    auto* from = base + (static_cast<std::size_t>(slot_max) * slab_length);
                    std::fill(from, from + slab_length, -kInfinity);
                }
            }
        }

        workers.Run(tasks, [&](std::size_t task, std::size_t) {
            double* partial = partials.data() + (task * partial_stride);
            for (std::uint64_t unit = task; unit < units; unit += tasks) {
                const std::uint64_t channel = unit / cells;
                const std::uint64_t cell = unit % cells;
                accumulate_unit(channel, chunk_cv_begin + (cell / cu_span),
                                chunk_cu_begin + (cell % cu_span), partial);
            }
        });

        for (std::size_t task = 0; task < tasks; ++task) {
            const double* partial = partials.data() + (task * partial_stride);
            for (std::size_t r = 0; r < request.region_count; ++r) {
                for (std::size_t slot = 0; slot < statistic_count; ++slot) {
                    const double* from =
                        partial + (r * partial_region_stride) + (slot * slab_length);
                    double* to = accumulator.data() + (r * region_stride) +
                                 (slot * statistic_stride) +
                                 static_cast<std::size_t>(slab.first_channel);
                    const bool is_min = slot_min >= 0 && slot == static_cast<std::size_t>(slot_min);
                    const bool is_max = slot_max >= 0 && slot == static_cast<std::size_t>(slot_max);
                    for (std::uint64_t channel = 0; channel < slab_length; ++channel) {
                        if (is_min) {
                            to[channel] = std::min(to[channel], from[channel]);
                        } else if (is_max) {
                            to[channel] = std::max(to[channel], from[channel]);
                        } else {
                            to[channel] += from[channel];
                        }
                    }
                }
            }
        }
    };

    // Which chunk runs this block's regions occupy, cut into bands of identical rows and then into
    // segments a budget wide. Two footprints of the same block share `reads_done`, so a block taken
    // in a single read still reports once -- at the end, through the emitter.
    const auto walk_block = [&](EmitBlock& block, const auto& report) -> Result<void> {
        for (std::uint64_t row = 0; row < occupancy.rows();) {
            const auto& runs = runs_per_row.at(static_cast<std::size_t>(row));
            if (runs.empty()) {
                ++row;
                continue;
            }

            // Rows below that repeat this row's runs exactly are read with it, so that a solid
            // rectangle becomes a few large requests while a diagonal stays one chunk per row.
            std::uint64_t widest = 0;
            for (const auto& run : runs) {
                widest = std::max(widest, run.last - run.first + 1);
            }
            const std::uint64_t band_limit = plan.UnitsAffordable(widest);
            std::uint64_t band_end = row + 1;
            while (band_end < occupancy.rows() && (band_end - row) < band_limit &&
                   runs_per_row.at(static_cast<std::size_t>(band_end)) == runs) {
                ++band_end;
            }

            const std::uint64_t chunk_cv_begin = occupancy.chunk_cv0() + row;
            const std::uint64_t chunk_cv_end = occupancy.chunk_cv0() + band_end;
            const std::uint64_t v_begin = std::max(occupancy.v0(), chunk_cv_begin * plan.chunk_v);
            const std::uint64_t v_end = std::min(occupancy.v1(), chunk_cv_end * plan.chunk_v);

            // A run wider than the budget is read in pieces, not in one request. Without this the
            // smallest request is a whole chunk row of the region: 157 chunks of a 80000-pixel-wide
            // region is 628 MiB, ten times the budget it was supposed to respect, and nothing to
            // report or cancel from until all of it lands.
            const std::uint64_t band_rows = std::max<std::uint64_t>(1, band_end - row);
            const std::uint64_t segment_limit = plan.UnitsAffordable(band_rows);
            segments.clear();
            for (const auto& whole : runs) {
                for (std::uint64_t first = whole.first; first <= whole.last;) {
                    const std::uint64_t width = std::min(segment_limit, whole.last - first + 1);
                    segments.push_back(ColumnRun{first, first + width - 1});
                    first += width;
                }
            }

            for (const auto& run : segments) {
                const std::uint64_t chunk_cu_begin = occupancy.chunk_cu0() + run.first;
                const std::uint64_t chunk_cu_end = occupancy.chunk_cu0() + run.last + 1;
                const std::uint64_t u_begin = std::max(occupancy.u0(), chunk_cu_begin * plan.chunk_u);
                const std::uint64_t u_end = std::min(occupancy.u1(), chunk_cu_end * plan.chunk_u);
                const std::uint64_t run_chunks =
                    std::max<std::uint64_t>(1, (run.last - run.first + 1) * (band_end - row));
                SlabFootprint footprint;
                footprint.u_start = u_begin;
                footprint.u_count = u_end - u_begin;
                footprint.v_start = v_begin;
                footprint.v_count = v_end - v_begin;
                footprint.chunks = run_chunks;

                FootprintBounds bounds;
                bounds.chunk_cv_begin = chunk_cv_begin;
                bounds.chunk_cv_end = chunk_cv_end;
                bounds.chunk_cu_begin = chunk_cu_begin;
                bounds.chunk_cu_end = chunk_cu_end;
                bounds.u_begin = u_begin;
                bounds.u_end = u_end;
                bounds.v_begin = v_begin;
                bounds.v_end = v_end;

                const auto walked =
                    walk.Over(footprint, block.begin, block.end, block.reads_done, block.chunks_done, report,
                              [&](const Slab& slab) { accumulate_slab(bounds, slab); });
                if (!walked) {
                    return walked.error();
                }
            }
            row = band_end;
        }
        return {};
    };

    const auto hand_over = [&](std::uint64_t first_channel, std::uint64_t length, bool complete,
                               double completeness) {
        settle_extrema(true, static_cast<std::size_t>(length));
        SpectralBlock block;
        block.first_channel = first_channel;
        block.channel_count = length;
        block.values = accumulator.data();
        block.value_count = accumulator.size();
        block.region_stride = region_stride;
        block.statistic_stride = statistic_stride;
        block.statistics = statistics.data();
        block.statistic_count = statistic_count;
        block.complete = complete;
        block.completeness = completeness;
        const bool keep_going = sink(block);
        if (!complete) {
            settle_extrema(false, static_cast<std::size_t>(length));
        }
        return keep_going;
    };

    return emitter.Over(reset_block, walk_block, hand_over);
}

}  // namespace carta::zarr::internal
