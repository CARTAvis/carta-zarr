/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "plane_histogram.h"

#include "axis_map.h"
#include "chunk_blocks.h"
#include "reduce/pass.h"
#include "reduce/plane_selection.h"
#include "reduce/provisional_histograms.h"
#include "reduce/tuning.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

namespace carta::zarr::internal {

// The one place a HistogramBlock is handed its counts. HistogramBlock names this its friend, so the
// layout Counts reads is set here and nowhere a caller can reach.
class HistogramBlocks {
public:
    // `counts` is [channel][bin], block.bin_count wide.
    static void Hold(HistogramBlock& block, const std::uint64_t* counts) noexcept { block._counts = counts; }
};

namespace {

// The one check both histograms make. It takes the count rather than a request because their two
// requests are different types.
Result<void> ValidateBins(const std::string& node, std::uint32_t bins) {
    if (bins == 0 || bins > kMaxHistogramBins) {
        return Error{ErrorCode::invalid_argument,
                     "A histogram needs between 1 and " + std::to_string(kMaxHistogramBins) + " bins", node};
    }
    return {};
}

// The bounds a fixed-range histogram bins against. A cube histogram has none: its range comes from
// the data's own extremes, so there is nothing here for it to be checked against.
Result<void> ValidateRange(const std::string& node, const HistogramRequest& request) {
    // A zero-width range would divide by zero on every pixel. The caller decides what an image with
    // no finite pixel should look like -- there is more than one defensible answer -- so this says
    // no rather than inventing one.
    if (!(request.lower < request.upper) || !std::isfinite(request.lower) || !std::isfinite(request.upper)) {
        return Error{ErrorCode::invalid_argument,
                     "A histogram needs a finite range with a lower bound below its upper bound", node};
    }
    // Pixels are float, so the range is narrowed once and every pixel is binned against the
    // narrowed copy. Bounds that do not fit a float leave nothing to bin against; the caller's own
    // bounds are its pixels' extremes, which always fit. Everything else that narrowing does to a
    // range -- two bounds narrowed to one float, a width underflowing to zero or overflowing to
    // infinity -- the caller bins all the same, and ComputeHistogram bins as it does.
    if (!std::isfinite(static_cast<float>(request.lower)) || !std::isfinite(static_cast<float>(request.upper))) {
        return Error{ErrorCode::invalid_argument, "A histogram needs bounds that fit in a float, as its pixels do",
                     node};
    }
    return {};
}

}  // namespace

Result<void> ComputeHistogram(const ReducibleImage& image, const HistogramRequest& request, const HistogramSink& sink,
                              const ReadOptions& options) {
    const auto& descriptor = image.descriptor();
    const auto& source = image.source();
    const auto& node = descriptor.id;
    if (!sink) {
        return Error{ErrorCode::invalid_argument, "A histogram needs a sink", node};
    }
    if (auto valid = ValidateBins(node, request.bins); !valid) {
        return valid.error();
    }
    if (auto valid = ValidateRange(node, request); !valid) {
        return valid.error();
    }
    const auto planned = image.Plan(request.planes, 1, options);
    if (!planned) {
        return planned.error();
    }
    const auto& plan = planned.value();

    auto pass = PassOverPlane(source, plan, options, "The histogram was cancelled by its sink");

    // The caller's own sequence: divide in double, narrow the width, compare against the narrowed
    // bounds. Doing any one of those in the other type moves pixels across bin edges.
    const float width = static_cast<float>((request.upper - request.lower) / request.bins);
    const float lower = static_cast<float>(request.lower);
    const float upper = static_cast<float>(request.upper);
    // Over a range whose span and width are positive and finite in float, every offset is finite and
    // at least zero, and the float sequence above is the whole of it. Over any other, an offset found
    // in float can be NaN or infinite, and converting one to an index is undefined; there the caller
    // bins every pixel in double, against a width found from the narrowed bounds, and a range of no
    // width puts every pixel it admits in the first bin. Decided once for the range, as the caller
    // decides it: deciding per pixel binned the lower half of a range wider than a float one way and
    // the upper half the other, which is neither of the caller's answers.
    const float span = upper - lower;
    const bool finite_offsets = span > 0.0F && std::isfinite(span) && width > 0.0F && std::isfinite(width);
    const double wide_lower = lower;
    const double wide_width = (static_cast<double>(upper) - wide_lower) / request.bins;
    const double wide_last = static_cast<double>(request.bins - 1);
    const auto bins = static_cast<std::size_t>(request.bins);
    std::vector<std::uint64_t> counts;

    // How many private histograms the binning may split a plane into. Capped three ways: by the
    // pool, by memory -- see kHistogramPartialBudgetBytes -- and, inside the visit where the plane's
    // size is known, by whether there is enough work to be worth waking anyone for.
    const auto split = image.Split(kHistogramPartialBudgetBytes, bins * sizeof(std::uint64_t));
    // One allocation for the whole plan. Each task owns one row of it, so no two of them ever touch
    // the same bin and the sum at the end is the only place they meet.
    std::vector<std::uint64_t> partials;
    if (split.most() > 1) {
        partials.resize(split.most() * bins);
    }

    // One read's worth of pixels binned into the block's counts. Named rather than written
    // into the call below, because it is the longest of the three lambdas the pass is handed
    // and written in place it would bury the other two.
    const auto bin_slab = [&](const Slab& slab) {
        // Hoisted into locals so that the loops below are the same text they were when the
        // pass handed these over as eight separate arguments.
        const std::uint64_t stride_u = slab.stride_u;
        const std::uint64_t stride_v = slab.stride_v;
        const std::uint64_t stride_z = slab.stride_z;
        const std::uint64_t u_count = slab.u_count;
        const std::uint64_t v_count = slab.v_count;
        for (std::uint64_t offset = 0; offset < slab.channel_count; ++offset) {
            const float* plane = slab.pixels + (offset * stride_z);
            std::uint64_t* into =
                counts.data() + (static_cast<std::size_t>((slab.first_channel + offset).index) * bins);

            // One loop per way of binning, each its own instantiation, so that the common one is the
            // loop it always was. The caller's own rule in both: a pixel outside the range is not
            // counted, and NaN fails both comparisons.
            const auto bin_rows_as = [&](auto finite, std::uint64_t v_first, std::uint64_t v_last,
                                         std::uint64_t* destination) {
                // Copied into locals whose addresses never escape. Reached through the captures, each
                // was loaded again for every pixel -- the store into a bin may, as far as the compiler
                // knows, have changed it -- and through two lambdas' captures, twice over.
                const float* const pixels = plane;
                const std::uint64_t row_stride = stride_v;
                const std::uint64_t column_stride = stride_u;
                const std::uint64_t columns = u_count;
                const float low = lower;
                const float high = upper;
                const float step = width;
                const std::size_t count = bins;
                const double wide_low = wide_lower;
                const double wide_step = wide_width;
                const double last = wide_last;
                for (std::uint64_t v = v_first; v < v_last; ++v) {
                    const float* row = pixels + (v * row_stride);
                    for (std::uint64_t u = 0; u < columns; ++u) {
                        const float value = row[u * column_stride];
                        if (low <= value && value <= high) {
                            if constexpr (decltype(finite)::value) {
                                auto bin = static_cast<std::size_t>((value - low) / step);
                                if (bin >= count) {
                                    bin = count - 1;
                                }
                                ++destination[bin];
                            } else {
                                const double offset =
                                    wide_step > 0.0 ? (static_cast<double>(value) - wide_low) / wide_step : 0.0;
                                ++destination[static_cast<std::size_t>(std::min(offset, last))];
                            }
                        }
                    }
                }
            };
            const auto bin_rows = [&](std::uint64_t v_first, std::uint64_t v_last, std::uint64_t* destination) {
                if (finite_offsets) {
                    bin_rows_as(std::true_type{}, v_first, v_last, destination);
                } else {
                    bin_rows_as(std::false_type{}, v_first, v_last, destination);
                }
            };

            // Rows, not planes: a read holding one plane is the common case for a large image,
            // so splitting by plane would leave the split with nothing to divide.
            const std::size_t tasks = split.Tasks(u_count, v_count);
            if (tasks <= 1) {
                bin_rows(0, v_count, into);
                continue;
            }

            std::fill(partials.begin(), partials.begin() + static_cast<std::ptrdiff_t>(tasks * bins), 0);
            split.Run(tasks, v_count, [&](std::size_t task, std::uint64_t first, std::uint64_t last) {
                bin_rows(first, last, partials.data() + (task * bins));
            });
            // Integer counts, so this sum is the serial loop's answer exactly -- which is what
            // lets histogram_test keep comparing against an oracle rather than a tolerance.
            for (std::size_t task = 0; task < tasks; ++task) {
                const std::uint64_t* from = partials.data() + (task * bins);
                for (std::size_t bin = 0; bin < bins; ++bin) {
                    into[bin] += from[bin];
                }
            }
        }
    };

    return pass.InBlocks(
        bins * sizeof(std::uint64_t), request.emit_every_channels,
        [&](std::uint64_t length) { counts.assign(static_cast<std::size_t>(length) * bins, 0); }, bin_slab,
        [&](SelectionChannel first_channel, std::uint64_t length, bool complete, double completeness) {
            HistogramBlock block;
            // Out of the type and into the public block, which is the one place it happens.
            block.first_channel = first_channel.index;
            block.channel_count = length;
            block.bin_count = bins;
            HistogramBlocks::Hold(block, counts.data());
            block.complete = complete;
            block.completeness = completeness;
            return sink(block);
        });
}

Result<CubeHistogramResult> ComputeCubeHistogram(const ReducibleImage& image, const CubeHistogramRequest& request,
                                                 const ReadOptions& options,
                                                 const CubeHistogramProgressCallback& progress) {
    const auto& descriptor = image.descriptor();
    const auto& source = image.source();
    const auto& node = descriptor.id;
    if (auto valid = ValidateBins(node, request.bins); !valid) {
        return valid.error();
    }
    if (request.spatial_sample == 0) {
        return Error{ErrorCode::invalid_argument, "A spatial sample of zero selects nothing", node};
    }
    const auto planned = image.Plan(request.planes, request.spatial_sample, options);
    if (!planned) {
        return planned.error();
    }
    const auto& plan = planned.value();

    ProvisionalHistograms histograms(image, request);

    auto pass = PassOverPlane(source, plan, options, "The histogram was cancelled by its caller");
    const auto walked = pass.Whole(
        [&](double fraction) {
            if (!progress) {
                return true;
            }
            CubeHistogramProgress update;
            update.progress = fraction;
            // By reference and lazily: re-aggregating on every read would cost more than the
            // binning does on a cube with thousands of them, and a caller that only draws a bar
            // never asks. Safe because the pass reports between reads, with no Add under way.
            update.snapshot = [&histograms]() { return histograms.Collect(); };
            return progress(update);
        },
        [&](const Slab& slab) { histograms.Add(slab); });
    if (!walked) {
        return walked.error();
    }

    return histograms.Collect();
}

}  // namespace carta::zarr::internal
