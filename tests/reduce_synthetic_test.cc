/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// A whole reduction, over pixels that were never written down.
//
// This target links no TensorStore and no Store. That is the point of it as much as the assertions
// are: until the pass took a PixelSource, a reduction could only be reached through a directory tree,
// so its arithmetic was pinned against the one committed fixture large enough to reach the parallel
// paths -- 512x520x4 -- and the oracle had to be written out beside it each time.
//
// Here the oracle is the formula the pixels come from, so the image can be whatever size the
// question needs.

#include "reduce/pass.h"
#include "reduce/plane_histogram.h"
#include "reduce/spectral_reduce.h"
#include "support/check.h"
#include "support/synthetic_pixel_source.h"
#include "work_pool.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ChunkGeometry;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadOptions;
using carta::zarr::internal::MapAxes;
using carta::zarr::internal::PlanPass;
using carta::zarr::internal::PlanRowTasks;
using carta::zarr::internal::ReducibleImage;
using carta::zarr::internal::WorkPool;
using carta::zarr::testing::SyntheticPixelSource;

using carta::zarr::testing::Require;

ReducibleImage Reducible(const SyntheticPixelSource& source, const ImageDescriptor& image,
                         const ChunkGeometry& geometry, WorkPool& workers) {
    auto reducible = ReducibleImage::Of(source, image, geometry, geometry, workers);
    Require(static_cast<bool>(reducible), "the synthetic image's axes could not be mapped");
    return reducible.value();
}

ImageDescriptor MakeImage(std::uint64_t x, std::uint64_t y, std::uint64_t channels) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, x},
             {"m", AxisRole::spatial_y, y},
             {"frequency", AxisRole::spectral, channels},
             {"polarization", AxisRole::polarization, 1},
             {"time", AxisRole::time, 1}};
    const std::size_t storage[]{3, 4, 1, 2, 0};
    for (std::size_t i = 0; i < 5; ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = axes[i].name;
        axis.role = axes[i].role;
        axis.length = axes[i].length;
        axis.storage_index = storage[i];
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

ChunkGeometry MakeGeometry(std::uint64_t cx, std::uint64_t cy, std::uint64_t cz) {
    ChunkGeometry geometry;
    geometry.fastest_spatial_axis = AxisRole::spatial_y;
    geometry.chunk_shape = {cx, cy, cz, 1, 1};
    return geometry;
}

// Values in [0, 1000), spread so that every bin of a coarse histogram gets something.
float Value(const std::vector<std::uint64_t>& logical) {
    return static_cast<float>(((logical.at(0) * 7) + (logical.at(1) * 13) + (logical.at(2) * 29)) % 1000);
}

constexpr std::uint64_t kX = 256;
constexpr std::uint64_t kY = 260;
constexpr std::uint64_t kZ = 8;

// A plane large enough for the binning to be worth splitting across workers.
//
// Both histograms hand a plane to PlanRowTasks and bin it in place when it comes back with one
// task, and the split is where a private accumulator per task is merged -- 9c21397 and 6ba0ce2 were
// both in that neighbourhood. kX by kY is 66,560 pixels, one under the 65,536 a task has to be
// worth, so every reduction in this file until now took the serial branch. The committed 512x520
// fixture was the only thing in the repository that reached the other one.
constexpr std::uint64_t kSplitX = 384;
constexpr std::uint64_t kSplitY = 390;
constexpr std::uint64_t kSplitZ = 4;
// Enough that one band covers the plane and one slab covers the channels, so what varies between
// the two runs below is the split and nothing else.
constexpr std::size_t kRoomyBudget = 16U << 20U;

// Every pixel of every plane binned, at a size no fixture reaches, against the formula rather than
// against a recorded answer.
void TestAHistogramCountsEveryPixel() {
    const auto image = MakeImage(kX, kY, kZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticPixelSource source(image, geometry, Value);

    carta::zarr::HistogramRequest request;
    request.planes.spectral = {0, kZ, 1};
    request.bins = 50;
    request.lower = 0.0;
    request.upper = 1000.0;

    ReadOptions options;
    options.read_budget_bytes = 4 * 64 * 65 * 2 * 4;
    WorkPool workers(4);

    std::vector<std::uint64_t> counts(request.bins * kZ, 0);
    std::uint64_t blocks = 0;
    const auto reducible = Reducible(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ComputeHistogram(
        reducible, request,
        [&](const carta::zarr::HistogramBlock& block) {
            if (!block.complete) {
                return true;
            }
            ++blocks;
            for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                for (std::size_t bin = 0; bin < block.bin_count; ++bin) {
                    counts.at(((block.first_channel + c) * request.bins) + bin) = block.Counts(c)[bin];
                }
            }
            return true;
        },
        options);
    Require(static_cast<bool>(outcome),
            std::string("the histogram failed: ") + (outcome ? "" : outcome.error().message));
    Require(blocks > 1,
            "this selection should emit more than one block; a pass that starts at a channel other than "
            "the first is the case that read the wrong channels before");

    // The oracle: the same rule, over the same formula.
    std::vector<std::uint64_t> expected(request.bins * kZ, 0);
    const float width = static_cast<float>((request.upper - request.lower) / request.bins);
    const float lower = static_cast<float>(request.lower);
    const float upper = static_cast<float>(request.upper);
    for (std::uint64_t z = 0; z < kZ; ++z) {
        for (std::uint64_t l = 0; l < kX; ++l) {
            for (std::uint64_t m = 0; m < kY; ++m) {
                const float value = Value({l, m, z, 0, 0});
                if (lower <= value && value <= upper) {
                    auto bin = static_cast<std::size_t>((value - lower) / width);
                    if (bin >= request.bins) {
                        bin = request.bins - 1;
                    }
                    ++expected.at((z * request.bins) + bin);
                }
            }
        }
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        Require(counts.at(i) == expected.at(i), "bin " + std::to_string(i) + ": expected " +
                                                    std::to_string(expected.at(i)) + ", got " +
                                                    std::to_string(counts.at(i)));
    }
    Require(source.most_hits_on_one_chunk() == 1, "and it decoded each chunk once getting there");
}

// Two overlapping regions over every channel, which is what the multi-region pass exists for.
// A range wider than a float can hold. Both bounds and the bin width survive narrowing, so the
// request is valid, but value - lower for a pixel in the upper half of it is not: 1e38 - (-3e38)
// rounds to infinity, and infinity converted to a bin index was undefined -- in practice every such
// pixel was counted in the last bin. Three columns: one at the lower bound, whose offset is zero;
// one in the middle, whose offset overflows but whose bin is 6; one at the upper bound, whose
// offset overflows too and whose bin is the last.
void TestARangeWiderThanAFloatStillBinsItsPixels() {
    constexpr std::uint64_t kColumns = 3;
    constexpr std::uint64_t kRows = 4;
    const float values[kColumns]{-3e38F, 1e38F, 3e38F};
    const auto image = MakeImage(kColumns, kRows, 1);
    const auto geometry = MakeGeometry(kColumns, kRows, 1);
    SyntheticPixelSource source(image, geometry,
                                [&](const std::vector<std::uint64_t>& logical) { return values[logical.at(0)]; });

    carta::zarr::HistogramRequest request;
    request.planes.spectral = {0, 1, 1};
    request.bins = 10;
    request.lower = -3e38;
    request.upper = 3e38;

    WorkPool workers(1);
    std::vector<std::uint64_t> counts;
    const auto reducible = Reducible(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ComputeHistogram(
        reducible, request,
        [&](const carta::zarr::HistogramBlock& block) {
            if (block.complete) {
                counts.assign(block.Counts(0), block.Counts(0) + block.bin_count);
            }
            return true;
        },
        ReadOptions{});
    Require(static_cast<bool>(outcome),
            std::string("the histogram failed: ") + (outcome ? "" : outcome.error().message));

    std::vector<std::uint64_t> expected(request.bins, 0);
    expected.at(0) = kRows;
    expected.at(6) = kRows;
    expected.at(9) = kRows;
    for (std::size_t bin = 0; bin < expected.size(); ++bin) {
        Require(counts.at(bin) == expected.at(bin), "bin " + std::to_string(bin) + ": expected " +
                                                        std::to_string(expected.at(bin)) + ", got " +
                                                        std::to_string(counts.at(bin)));
    }
}

// carta-backend's Histogram::Fill, which is the answer the library's histogram stands in for. Over a
// range whose float offsets are all finite -- positive span and bin width, neither overflowing -- a
// pixel is binned in float, against the narrowed bounds and width. Over any other range, every
// pixel is binned in double, against a width found from the narrowed bounds, and a range of no width
// puts every pixel it admits in the first bin. The choice is made once for the range, not per pixel.
std::vector<std::uint64_t> CallerCounts(const std::vector<float>& values, double range_lower, double range_upper,
                                        std::size_t bins) {
    const float lower = static_cast<float>(range_lower);
    const float upper = static_cast<float>(range_upper);
    const float width = static_cast<float>((range_upper - range_lower) / static_cast<double>(bins));
    const float span = upper - lower;
    const bool finite_offsets = span > 0 && std::isfinite(span) && width > 0 && std::isfinite(width);
    const double wide_width = (static_cast<double>(upper) - static_cast<double>(lower)) / static_cast<double>(bins);
    std::vector<std::uint64_t> counts(bins, 0);
    for (const float value : values) {
        if (!(lower <= value && value <= upper)) {
            continue;
        }
        std::size_t bin = 0;
        if (finite_offsets) {
            bin = std::min(static_cast<std::size_t>((value - lower) / width), bins - 1);
        } else {
            const double offset =
                wide_width > 0 ? (static_cast<double>(value) - static_cast<double>(lower)) / wide_width : 0.0;
            bin = static_cast<std::size_t>(std::min(offset, static_cast<double>(bins - 1)));
        }
        ++counts.at(bin);
    }
    return counts;
}

// Every range the caller bins, binned as the caller bins it. The library refused the two whose float
// width does not survive narrowing -- one underflows to zero, one overflows to infinity -- although
// the caller counts both; and over a range wider than a float holds, it switched to double only for
// the pixels whose float offset overflowed, so pixels below the middle were binned one way and those
// above it another, and the counts differed from the caller's.
void TestEveryRangeIsBinnedAsTheCallerBinsIt() {
    struct Case {
        const char* what;
        std::vector<float> values;
        double lower;
        double upper;
        std::uint32_t bins;
    };
    const float tiny = std::numeric_limits<float>::denorm_min();
    const std::vector<Case> cases{
        {"a bin width that underflows a float", {0.0F, tiny, 2 * tiny, 7 * tiny}, 0.0, 1e-44, 100},
        {"a bin width that overflows a float", {-3e38F, 0.0F, 3e38F}, -3e38, 3e38, 1},
        {"a range wider than a float, in three bins",
         {-3e38F, -2.5e38F, -1e38F, -1e30F, 0.0F, 1e38F, 1.0000001e38F, 2e38F, 2.9e38F, 3e38F},
         -3e38,
         3e38,
         3},
        {"a range wider than a float, in ten bins", {-3e38F, -1.2e38F, 6e37F, 1e38F, 3e38F}, -3e38, 3e38, 10},
        {"a range narrower than a float", {1.0F, 2.0F}, 1.0, 1.0 + 1e-9, 8},
        {"an ordinary range", {0.0F, 0.1F, 0.3F, 0.7F, 0.99999994F, 1.0F, 1.5F}, 0.0, 1.0, 10},
    };
    for (const auto& c : cases) {
        const auto columns = static_cast<std::uint64_t>(c.values.size());
        const auto image = MakeImage(columns, 1, 1);
        const auto geometry = MakeGeometry(columns, 1, 1);
        SyntheticPixelSource source(image, geometry, [&](const std::vector<std::uint64_t>& logical) {
            return c.values.at(static_cast<std::size_t>(logical.at(0)));
        });
        carta::zarr::HistogramRequest request;
        request.planes.spectral = {0, 1, 1};
        request.bins = c.bins;
        request.lower = c.lower;
        request.upper = c.upper;

        WorkPool workers(1);
        std::vector<std::uint64_t> counts;
        const auto reducible = Reducible(source, image, geometry, workers);
        const auto outcome = carta::zarr::internal::ComputeHistogram(
            reducible, request,
            [&](const carta::zarr::HistogramBlock& block) {
                if (block.complete) {
                    counts.assign(block.Counts(0), block.Counts(0) + block.bin_count);
                }
                return true;
            },
            ReadOptions{});
        Require(static_cast<bool>(outcome),
                std::string(c.what) + ": the histogram failed: " + (outcome ? "" : outcome.error().message));
        const auto expected = CallerCounts(c.values, c.lower, c.upper, c.bins);
        for (std::size_t bin = 0; bin < expected.size(); ++bin) {
            Require(counts.at(bin) == expected.at(bin), std::string(c.what) + ", bin " + std::to_string(bin) +
                                                            ": the caller counts " + std::to_string(expected.at(bin)) +
                                                            ", the library " + std::to_string(counts.at(bin)));
        }
    }
}

void TestASpectralReductionAgreesWithTheFormula() {
    const auto image = MakeImage(kX, kY, kZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticPixelSource source(image, geometry, Value);

    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, kX, kY},
        {40, 50, 100, 120},
    };
    carta::zarr::SpectralReduceRequest request;
    request.planes.spectral = {0, kZ, 1};
    request.regions = {regions.data(), regions.size()};
    request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::sum |
                         carta::zarr::Statistic::min | carta::zarr::Statistic::max;

    ReadOptions options;
    options.read_budget_bytes = 4 * 64 * 65 * 2 * 4;
    WorkPool workers(4);

    std::vector<double> pixels(regions.size() * kZ, 0.0);
    std::vector<double> sums(regions.size() * kZ, 0.0);
    std::vector<double> minima(regions.size() * kZ, 0.0);
    std::vector<double> maxima(regions.size() * kZ, 0.0);
    const auto reducible = Reducible(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ReduceSpectral(
        reducible, request,
        [&](const carta::zarr::SpectralBlock& block) {
            if (!block.complete) {
                return true;
            }
            for (std::size_t r = 0; r < regions.size(); ++r) {
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    const auto totals = block.Totals(r, c);
                    const auto at = (r * kZ) + block.first_channel + c;
                    pixels.at(at) = totals.num_pixels;
                    sums.at(at) = totals.sum;
                    minima.at(at) = totals.min;
                    maxima.at(at) = totals.max;
                }
            }
            return true;
        },
        options);
    Require(static_cast<bool>(outcome),
            std::string("the reduction failed: ") + (outcome ? "" : outcome.error().message));

    for (std::size_t r = 0; r < regions.size(); ++r) {
        const auto& region = regions.at(r);
        for (std::uint64_t z = 0; z < kZ; ++z) {
            double count = 0.0;
            double sum = 0.0;
            double smallest = std::numeric_limits<double>::infinity();
            double largest = -std::numeric_limits<double>::infinity();
            for (std::uint64_t l = region.x_start; l < region.x_start + region.width; ++l) {
                for (std::uint64_t m = region.y_start; m < region.y_start + region.height; ++m) {
                    const double value = Value({l, m, z, 0, 0});
                    count += 1.0;
                    sum += value;
                    smallest = std::min(smallest, value);
                    largest = std::max(largest, value);
                }
            }
            const auto at = (r * kZ) + z;
            Require(pixels.at(at) == count, "region " + std::to_string(r) + " pixel count");
            Require(std::abs(sums.at(at) - sum) <= 1e-9 * (1.0 + std::abs(sum)),
                    "region " + std::to_string(r) + " sum");
            Require(minima.at(at) == smallest, "region " + std::to_string(r) + " minimum");
            Require(maxima.at(at) == largest, "region " + std::to_string(r) + " maximum");
        }
    }
    Require(source.most_hits_on_one_chunk() == 1, "and it decoded each chunk once getting there");
}

// A region occupies the chunks its mask touches, not the chunks its bounding box covers.
//
// This is what the reduction's chunk index is for, and until now the only thing asserting it was a
// pair of cases on a 4 x 5 x 2 x 3 fixture that ran a reduction, saw that no hand-over was left
// unfinished, and reasoned backwards to "so it read one chunk" -- both saying in their own comments
// that they stop testing anything if the chunk shape ever changes. Here the source is asked what it
// was given instead, at a size where the box and the occupancy are four times apart.
//
// Occupancy's own tests say the index is right. This says the reduction walks it: an index that is
// correct and then ignored reads all sixteen chunks and still returns the right statistics.
void TestAMaskedRegionReadsOnlyTheChunksItOccupies() {
    constexpr std::uint64_t kSide = 256;
    constexpr std::uint64_t kChunk = 64;
    constexpr std::uint64_t kChannels = 4;
    constexpr std::uint64_t kGrid = kSide / kChunk;  // 4 x 4 chunks in the plane

    const auto image = MakeImage(kSide, kSide, kChannels);
    const auto geometry = MakeGeometry(kChunk, kChunk, kChannels);

    // Set only where the chunk is on the diagonal of the chunk grid. The bounding box is the whole
    // plane; the selection is four of its sixteen chunks.
    std::vector<std::uint8_t> raster(static_cast<std::size_t>(kSide * kSide), 0);
    for (std::uint64_t y = 0; y < kSide; ++y) {
        for (std::uint64_t x = 0; x < kSide; ++x) {
            if (x / kChunk == y / kChunk) {
                raster.at(static_cast<std::size_t>((y * kSide) + x)) = 1;
            }
        }
    }

    // What the diagonal comes to, from the formula the pixels come from.
    std::vector<double> expected_sum(kChannels, 0.0);
    std::vector<double> expected_count(kChannels, 0.0);
    for (std::uint64_t z = 0; z < kChannels; ++z) {
        for (std::uint64_t y = 0; y < kSide; ++y) {
            for (std::uint64_t x = 0; x < kSide; ++x) {
                if (x / kChunk == y / kChunk) {
                    expected_sum.at(z) += Value({x, y, z, 0, 0});
                    expected_count.at(z) += 1.0;
                }
            }
        }
    }

    const auto reduce = [&](const carta::zarr::RegionMask& region, const char* described) {
        SyntheticPixelSource source(image, geometry, Value);
        const std::vector<carta::zarr::RegionMask> regions{region};

        carta::zarr::SpectralReduceRequest request;
        request.planes.spectral = {0, kChannels, 1};
        request.regions = {regions.data(), regions.size()};
        request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::sum;

        ReadOptions options;
        options.read_budget_bytes = kRoomyBudget;
        WorkPool workers(4);

        std::vector<double> sums(kChannels, 0.0);
        std::vector<double> counts(kChannels, 0.0);
        const auto reducible = Reducible(source, image, geometry, workers);
        const auto outcome = carta::zarr::internal::ReduceSpectral(
            reducible, request,
            [&](const carta::zarr::SpectralBlock& block) {
                if (!block.complete) {
                    return true;
                }
                const double* block_sums = block.Series(0, carta::zarr::Statistic::sum);
                const double* block_counts = block.Series(0, carta::zarr::Statistic::num_pixels);
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    sums.at(block.first_channel + c) = block_sums[c];
                    counts.at(block.first_channel + c) = block_counts[c];
                }
                return true;
            },
            options);
        Require(static_cast<bool>(outcome),
                std::string(described) + ": the reduction failed: " + (outcome ? "" : outcome.error().message));

        for (std::uint64_t z = 0; z < kChannels; ++z) {
            Require(counts.at(z) == expected_count.at(z), std::string(described) + ": the diagonal's pixel count");
            Require(std::abs(sums.at(z) - expected_sum.at(z)) <= 1e-9 * (1.0 + std::abs(expected_sum.at(z))),
                    std::string(described) + ": the diagonal's sum");
        }

        // The assertion the fixture-driven pair could only reach by inference. Four chunks of the
        // grid, one channel chunk, each decoded once.
        Require(source.chunks_touched() == kGrid,
                std::string(described) + ": read " + std::to_string(source.chunks_touched()) +
                    " chunks rather than the " + std::to_string(kGrid) + " the mask occupies");
        Require(source.most_hits_on_one_chunk() == 1, std::string(described) + ": a chunk was decoded more than once");
    };

    // This image varies y fastest, so the runs the reduction makes from the raster are ranges of
    // rows, one line per column -- the case that used to cost a column walk down the raster.
    carta::zarr::RegionMask rastered{0, 0, kSide, kSide, {raster.data(), raster.size()}};
    reduce(rastered, "as a raster");
}

// The same counts whether the plane is binned in place or in four pieces that are added up.
//
// Integer counts, so this is exact rather than close: a split that divided the rows wrongly, or a
// merge that dropped or double-counted a worker's partial, cannot agree with the serial answer to
// any tolerance at all.
void TestAPlaneHistogramSplitAcrossWorkers() {
    Require(PlanRowTasks(kSplitY, kSplitX, 4, 1U << 16U) > 1,
            "this plane is meant to be large enough to split; if it is not, this test proves nothing");

    const auto image = MakeImage(kSplitX, kSplitY, kSplitZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticPixelSource source(image, geometry, Value);

    carta::zarr::HistogramRequest request;
    request.planes.spectral = {0, kSplitZ, 1};
    request.bins = 37;
    request.lower = 0.0;
    request.upper = 1000.0;

    ReadOptions options;
    options.read_budget_bytes = kRoomyBudget;

    const auto counts_from = [&](std::size_t worker_count) {
        WorkPool workers(worker_count);
        const auto reducible = Reducible(source, image, geometry, workers);
        std::vector<std::uint64_t> counts(request.bins * kSplitZ, 0);
        const auto outcome = carta::zarr::internal::ComputeHistogram(
            reducible, request,
            [&](const carta::zarr::HistogramBlock& block) {
                if (!block.complete) {
                    return true;
                }
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    for (std::size_t bin = 0; bin < block.bin_count; ++bin) {
                        counts.at(((block.first_channel + c) * request.bins) + bin) = block.Counts(c)[bin];
                    }
                }
                return true;
            },
            options);
        Require(static_cast<bool>(outcome),
                std::string("the histogram failed: ") + (outcome ? "" : outcome.error().message));
        return counts;
    };

    // One worker cannot split -- PlanRowTasks refuses before the plane is even measured -- so this
    // is the serial branch, and the oracle the other one has to match.
    const auto in_place = counts_from(1);
    const auto split = counts_from(4);
    for (std::size_t i = 0; i < in_place.size(); ++i) {
        Require(in_place.at(i) == split.at(i), "bin " + std::to_string(i) + ": binning in place gave " +
                                                   std::to_string(in_place.at(i)) + ", four workers gave " +
                                                   std::to_string(split.at(i)));
    }

    std::uint64_t total = 0;
    for (const auto count : split) {
        total += count;
    }
    Require(total == kSplitX * kSplitY * kSplitZ, "every pixel is inside the range, so every one is counted");
}

// The cube histogram's own split, which is over rows of a whole read rather than of one plane, and
// which merges four provisional histograms of its own rather than four rows of counts.
//
// The counts depend on the thread count by design -- each accumulator re-aggregates its own
// provisional range onto the target grid -- so what is exact here is everything around them: the
// pixel count, the extremes, and that the bins still hold every pixel that was binned.
void TestACubeHistogramSplitAcrossWorkers() {
    // The budget below leaves one read holding every plane, so the rows this one splits are the
    // rows of the whole selection rather than of a plane.
    Require(PlanRowTasks(kSplitY, kSplitZ * kSplitX, 4, 1U << 16U) > 1,
            "this selection is meant to be large enough to split; if it is not, this test proves nothing");

    const auto image = MakeImage(kSplitX, kSplitY, kSplitZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticPixelSource source(image, geometry, Value);

    carta::zarr::CubeHistogramRequest request;
    request.planes.spectral = {0, kSplitZ, 1};
    request.bins = 64;

    ReadOptions options;
    options.read_budget_bytes = kRoomyBudget;
    WorkPool workers(4);

    const auto reducible = Reducible(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ComputeCubeHistogram(reducible, request, options, {});
    Require(static_cast<bool>(outcome),
            std::string("the cube histogram failed: ") + (outcome ? "" : outcome.error().message));
    const auto& result = outcome.value();

    double sum = 0.0;
    double smallest = std::numeric_limits<double>::infinity();
    double largest = -std::numeric_limits<double>::infinity();
    for (std::uint64_t z = 0; z < kSplitZ; ++z) {
        for (std::uint64_t l = 0; l < kSplitX; ++l) {
            for (std::uint64_t m = 0; m < kSplitY; ++m) {
                const double value = Value({l, m, z, 0, 0});
                sum += value;
                smallest = std::min(smallest, value);
                largest = std::max(largest, value);
            }
        }
    }

    const auto pixels = static_cast<double>(kSplitX * kSplitY * kSplitZ);
    Require(result.totals.num_pixels == pixels, "every pixel is finite, so every one counts");
    Require(result.totals.nan_count == 0.0, "and none of them is a NaN");
    Require(result.totals.min == smallest && result.totals.max == largest,
            "the extremes are tracked exactly, whatever the split");
    Require(std::abs(result.totals.sum - sum) <= 1e-9 * (1.0 + std::abs(sum)), "the sum over four accumulators");

    std::uint64_t binned = 0;
    for (const auto count : result.counts) {
        binned += count;
    }
    Require(static_cast<double>(binned) == pixels,
            "four provisional histograms re-aggregated onto one grid still hold every pixel");
}

// The spread of pixels far from zero against it, which is what sum and sum_sq cannot give: the cases
// ADR 0018 was measured on, and one where the pixel a span takes its distances from is the odd one out.

// About normal, of unit variance, and the same every time for the same pixel.
double Noise(const std::vector<std::uint64_t>& logical) {
    std::uint64_t state = (logical.at(0) * 0x9E3779B97F4A7C15ULL) ^ (logical.at(1) * 0xC2B2AE3D27D4EB4FULL) ^
                          (logical.at(2) * 0x165667B19E3779F9ULL);
    double total = 0.0;
    for (int draw = 0; draw < 4; ++draw) {
        state += 0x9E3779B97F4A7C15ULL;
        std::uint64_t mixed = state;
        mixed = (mixed ^ (mixed >> 30U)) * 0xBF58476D1CE4E5B9ULL;
        mixed = (mixed ^ (mixed >> 27U)) * 0x94D049BB133111EBULL;
        mixed ^= mixed >> 31U;
        total += static_cast<double>(mixed >> 11U) / 9007199254740992.0;
    }
    return (total - 2.0) * std::sqrt(3.0);
}

struct Spread {
    const char* name;
    SyntheticPixelSource::Formula value;
};

std::vector<Spread> Spreads() {
    return {
        {"all 1e8", [](const std::vector<std::uint64_t>&) { return 1.0e8F; }},
        {"1e6 + 0.5 noise",
         [](const std::vector<std::uint64_t>& at) { return static_cast<float>(1.0e6 + 0.5 * Noise(at)); }},
        {"1e7 + 0.5 noise",
         [](const std::vector<std::uint64_t>& at) { return static_cast<float>(1.0e7 + 0.5 * Noise(at)); }},
        {"1e8 + 20 noise",
         [](const std::vector<std::uint64_t>& at) { return static_cast<float>(1.0e8 + 20.0 * Noise(at)); }},
        // The walk's rows run along m, so every pixel at the start of a chunk's row is the outlier.
        {"1e8 + 20 noise, the first of each row far out",
         [](const std::vector<std::uint64_t>& at) {
             return static_cast<float>(1.0e8 + 20.0 * Noise(at) + (at.at(1) % 65 == 0 ? 1.0e5 : 0.0));
         }},
    };
}

// The sum of squared deviations over the pixels `in` selects, worked out the slow way: in long
// double, from a mean found first.
template <typename Selected>
double CentredReference(const SyntheticPixelSource::Formula& value, std::uint64_t channels_from,
                        std::uint64_t channels_to, std::uint64_t x, std::uint64_t y, const Selected& in) {
    long double count = 0.0L;
    long double sum = 0.0L;
    for (std::uint64_t z = channels_from; z < channels_to; ++z) {
        for (std::uint64_t l = 0; l < x; ++l) {
            for (std::uint64_t m = 0; m < y; ++m) {
                if (in(l, m)) {
                    count += 1.0L;
                    sum += value({l, m, z, 0, 0});
                }
            }
        }
    }
    const long double mean = sum / count;
    long double deviations = 0.0L;
    for (std::uint64_t z = channels_from; z < channels_to; ++z) {
        for (std::uint64_t l = 0; l < x; ++l) {
            for (std::uint64_t m = 0; m < y; ++m) {
                if (in(l, m)) {
                    const long double distance = value({l, m, z, 0, 0}) - mean;
                    deviations += distance * distance;
                }
            }
        }
    }
    return static_cast<double>(deviations);
}

void RequireCentred(double actual, double expected, const std::string& where) {
    Require(std::abs(actual - expected) <= 1e-10 * std::max(1.0, std::abs(expected)),
            where + ": sum_sq_dev " + std::to_string(actual) + ", the pixels' own " + std::to_string(expected));
}

void TestASpectralSpreadIsCentredWhateverTheMagnitude() {
    const auto image = MakeImage(kX, kY, kZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    // A region with holes in it, so that its rows are several spans and take the masked loop.
    constexpr std::uint64_t kHolesX = 150;
    constexpr std::uint64_t kHolesY = 170;
    std::vector<std::uint8_t> holes(kHolesX * kHolesY);
    for (std::uint64_t i = 0; i < holes.size(); ++i) {
        holes[i] = static_cast<std::uint8_t>(((i / kHolesX) % 7 != 3) && ((i % kHolesX) % 11 != 5));
    }
    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, kX, kY},
        {40, 50, 100, 120},
        {30, 20, kHolesX, kHolesY, {holes.data(), holes.size()}},
    };
    const auto inside = [&](std::size_t r, std::uint64_t l, std::uint64_t m) {
        const auto& region = regions.at(r);
        if (l < region.x_start || l >= region.x_start + region.width || m < region.y_start ||
            m >= region.y_start + region.height) {
            return false;
        }
        return region.mask.size == 0 ||
               region.mask.data[((m - region.y_start) * region.width) + (l - region.x_start)] != 0;
    };

    for (const auto& spread : Spreads()) {
        SyntheticPixelSource source(image, geometry, spread.value);
        carta::zarr::SpectralReduceRequest request;
        request.planes.spectral = {0, kZ, 1};
        request.regions = {regions.data(), regions.size()};
        // Alone: the count and the sum come with it.
        request.statistics = carta::zarr::Statistic::sum_sq_dev;

        ReadOptions options;
        options.read_budget_bytes = 4 * 64 * 65 * 2 * 4;
        WorkPool workers(4);
        std::vector<double> deviations(regions.size() * kZ, -1.0);
        const auto reducible = Reducible(source, image, geometry, workers);
        const auto outcome = carta::zarr::internal::ReduceSpectral(
            reducible, request,
            [&](const carta::zarr::SpectralBlock& block) {
                Require(block.Carries(carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::sum |
                                      carta::zarr::Statistic::sum_sq_dev),
                        "a block asked for sum_sq_dev carries the count and the sum it was made with");
                if (block.complete) {
                    for (std::size_t r = 0; r < regions.size(); ++r) {
                        for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                            deviations.at((r * kZ) + block.first_channel + c) = block.Totals(r, c).sum_sq_dev;
                        }
                    }
                }
                return true;
            },
            options);
        Require(static_cast<bool>(outcome),
                std::string("the reduction failed: ") + (outcome ? "" : outcome.error().message));

        for (std::size_t r = 0; r < regions.size(); ++r) {
            for (std::uint64_t z = 0; z < kZ; ++z) {
                const auto expected = CentredReference(
                    spread.value, z, z + 1, kX, kY, [&](std::uint64_t l, std::uint64_t m) { return inside(r, l, m); });
                RequireCentred(
                    deviations.at((r * kZ) + z), expected,
                    std::string(spread.name) + ", region " + std::to_string(r) + ", channel " + std::to_string(z));
            }
        }
    }
}

void TestACubeSpreadIsCentredWhateverTheMagnitude() {
    const auto image = MakeImage(kSplitX, kSplitY, kSplitZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    for (const auto& spread : Spreads()) {
        SyntheticPixelSource source(image, geometry, spread.value);
        carta::zarr::CubeHistogramRequest request;
        request.planes.spectral = {0, kSplitZ, 1};
        request.bins = 64;
        ReadOptions options;
        options.read_budget_bytes = 4 * 64 * 65 * 2 * 4;
        WorkPool workers(4);
        const auto reducible = Reducible(source, image, geometry, workers);
        const auto outcome = carta::zarr::internal::ComputeCubeHistogram(reducible, request, options, {});
        Require(static_cast<bool>(outcome),
                std::string("the cube histogram failed: ") + (outcome ? "" : outcome.error().message));
        RequireCentred(outcome.value().totals.sum_sq_dev,
                       CentredReference(spread.value, 0, kSplitZ, kSplitX, kSplitY,
                                        [](std::uint64_t, std::uint64_t) { return true; }),
                       std::string(spread.name) + ", the cube");
    }
}

// Progress is a count of chunks read against a count of chunks to read, and the two are worked out
// in different places: the walk adds up what each read covered, the total is taken once over the
// whole run. Every report is made before a read that is still to come, so none of them may say the
// whole run is done.
//
// Each spectral selection below starts off a chunk boundary or steps across one, which is where a
// count taken over the whole run and a sum taken read by read can disagree. A one-byte budget makes
// every chunk layer its own read, so there are reports to look at.
void TestProgressNeverClaimsTheWholeRunBeforeItsLastRead() {
    constexpr std::uint64_t kSide = 8;
    constexpr std::uint64_t kDepth = 24;
    const auto image = MakeImage(kSide, kSide, kDepth);
    const auto geometry = MakeGeometry(4, 4, 4);
    SyntheticPixelSource source(image, geometry, Value);
    WorkPool workers(1);
    const auto reducible = Reducible(source, image, geometry, workers);

    ReadOptions options;
    options.read_budget_bytes = 1;

    // Run to the end, and stopped short of it: a run whose ends both fall inside a chunk is the one
    // that covers a chunk more than its length suggests.
    const Range runs[]{{0, 24, 1}, {1, 23, 1}, {1, 8, 1}, {3, 6, 1}, {0, 8, 3},
                       {1, 8, 3},  {3, 7, 3},  {1, 4, 3}, {0, 3, 9}};
    std::string wrong;
    for (const auto& spectral : runs) {
        const auto start = spectral.start;
        const auto stride = spectral.stride;
        const auto count = spectral.count;
        const std::string which = " (start " + std::to_string(start) + ", count " + std::to_string(count) +
                                  ", stride " + std::to_string(stride) + ")";

        carta::zarr::CubeHistogramRequest cube;
        cube.planes.spectral = spectral;
        cube.bins = 16;
        std::vector<double> reported;
        const auto histogram = carta::zarr::internal::ComputeCubeHistogram(
            reducible, cube, options, [&](const carta::zarr::CubeHistogramProgress& update) {
                reported.push_back(update.progress);
                return true;
            });
        Require(static_cast<bool>(histogram), "the cube histogram failed" + which);
        Require(!reported.empty(), "a one-byte budget should make the cube histogram report" + which);
        for (const double progress : reported) {
            if (progress >= 1.0) {
                wrong +=
                    "\n  a cube histogram reported " + std::to_string(progress) + " with a read still to come" + which;
                break;
            }
        }

        // One block for the whole run, so that a block is several reads.
        const carta::zarr::RegionMask whole{0, 0, kSide, kSide};
        carta::zarr::SpectralReduceRequest reduce;
        reduce.planes.spectral = spectral;
        reduce.regions = {&whole, 1};
        reduce.statistics = carta::zarr::Statistic::sum;
        reduce.emit_every_channels = static_cast<std::uint32_t>(count);
        const auto reduced = carta::zarr::internal::ReduceSpectral(
            reducible, reduce,
            [&](const carta::zarr::SpectralBlock& block) {
                if (!block.complete && block.completeness >= 1.0) {
                    wrong +=
                        "\n  an unfinished block claimed " + std::to_string(block.completeness) + " of itself" + which;
                }
                return true;
            },
            options);
        Require(static_cast<bool>(reduced), "the spectral reduction failed" + which);
    }
    Require(wrong.empty(), "progress ran ahead of the reads:" + wrong);
}

// ---------------------------------------------------------------------------------------------
// When a reduction hands over, and what it says about how far along it is, exactly.
//
// The tests above bound these: progress never claims the whole run early, a block is complete only
// at its end. These pin them, read by read, at sizes small enough to count by hand, because the
// numbers come from three places that have to agree -- the walk counts the chunks it has read, the
// block or the run says what that is a fraction of, and which layer it is a fraction of depends on
// whether the walk is over a plane or over a region set. Every value below is worked out in the
// comment above it.

// One hand-over, as a sink saw it.
struct HandOver {
    std::uint64_t first_channel = 0;
    std::uint64_t channel_count = 0;
    bool complete = false;
    double completeness = 0.0;

    bool operator==(const HandOver& other) const {
        return first_channel == other.first_channel && channel_count == other.channel_count &&
               complete == other.complete && completeness == other.completeness;
    }
};

std::string Describe(const std::vector<HandOver>& handed) {
    std::string text;
    for (const auto& one : handed) {
        text += "\n  [" + std::to_string(one.first_channel) + ", +" + std::to_string(one.channel_count) + ") " +
                (one.complete ? "complete" : "partial") + " " + std::to_string(one.completeness);
    }
    return text;
}

void RequireHandOvers(const std::vector<HandOver>& handed, const std::vector<HandOver>& expected,
                      const std::string& what) {
    Require(handed == expected, what + ": handed over" + Describe(handed) + "\nrather than" + Describe(expected));
}

std::vector<HandOver> ReduceAndRecord(const SyntheticPixelSource& source, const ImageDescriptor& image,
                                      const ChunkGeometry& geometry,
                                      const std::vector<carta::zarr::RegionMask>& regions, const Range& spectral,
                                      const ReadOptions& options) {
    WorkPool workers(1);
    const auto reducible = Reducible(source, image, geometry, workers);
    carta::zarr::SpectralReduceRequest request;
    request.planes.spectral = spectral;
    request.regions = {regions.data(), regions.size()};
    request.statistics = carta::zarr::Statistic::sum;
    std::vector<HandOver> handed;
    const auto outcome = carta::zarr::internal::ReduceSpectral(
        reducible, request,
        [&](const carta::zarr::SpectralBlock& block) {
            handed.push_back({block.first_channel, block.channel_count, block.complete, block.completeness});
            return true;
        },
        options);
    Require(static_cast<bool>(outcome),
            std::string("the reduction failed: ") + (outcome ? "" : outcome.error().message));
    return handed;
}

// An 8 x 8 plane of 4 x 4 chunks, eight channels deep in chunks of four, read one chunk at a time.
//
// The plane histogram's blocks are the plan's: its layer is the plane's four chunks, and a budget
// below one chunk affords one layer of them, which is one spectral chunk -- four channels. So two
// blocks, each read as four one-chunk reads, each handed over before every read after its first,
// at the chunks it has of the four it covers, and once more finished.
void TestAPlaneHistogramHandsOverAtEveryReadOfItsBlock() {
    const auto image = MakeImage(8, 8, 8);
    const auto geometry = MakeGeometry(4, 4, 4);
    SyntheticPixelSource source(image, geometry, Value);
    WorkPool workers(1);
    const auto reducible = Reducible(source, image, geometry, workers);

    carta::zarr::HistogramRequest request;
    request.planes.spectral = {0, 8, 1};
    request.bins = 16;
    request.lower = 0.0;
    request.upper = 1000.0;
    ReadOptions options;
    options.read_budget_bytes = 1;

    std::vector<HandOver> handed;
    const auto outcome = carta::zarr::internal::ComputeHistogram(
        reducible, request,
        [&](const carta::zarr::HistogramBlock& block) {
            handed.push_back({block.first_channel, block.channel_count, block.complete, block.completeness});
            return true;
        },
        options);
    Require(static_cast<bool>(outcome), "the histogram failed");
    RequireHandOvers(handed,
                     {{0, 4, false, 0.25},
                      {0, 4, false, 0.5},
                      {0, 4, false, 0.75},
                      {0, 4, true, 1.0},
                      {4, 4, false, 0.25},
                      {4, 4, false, 0.5},
                      {4, 4, false, 0.75},
                      {4, 4, true, 1.0}},
                     "a plane histogram");
    Require(source.pixel_reads() == 8, "eight one-chunk reads, not " + std::to_string(source.pixel_reads()));
}

// The same image and budget, reduced over a region covering the whole plane. A region set's layer
// is the chunks it occupies, which here is the plane's four, and the occupancy cuts them into four
// one-chunk footprints because that is all a read affords. A block is walked footprint by footprint,
// and the four footprints of one block share one count of reads: they report as the plane's four
// bands do, not once per footprint.
void TestARegionCoveringThePlaneHandsOverAsThePlaneDoes() {
    const auto image = MakeImage(8, 8, 8);
    const auto geometry = MakeGeometry(4, 4, 4);
    SyntheticPixelSource source(image, geometry, Value);
    ReadOptions options;
    options.read_budget_bytes = 1;

    const auto handed = ReduceAndRecord(source, image, geometry, {{0, 0, 8, 8}}, Range{0, 8, 1}, options);
    RequireHandOvers(handed,
                     {{0, 4, false, 0.25},
                      {0, 4, false, 0.5},
                      {0, 4, false, 0.75},
                      {0, 4, true, 1.0},
                      {4, 4, false, 0.25},
                      {4, 4, false, 0.5},
                      {4, 4, false, 0.75},
                      {4, 4, true, 1.0}},
                     "a region covering the plane");
    Require(source.pixel_reads() == 8, "eight one-chunk reads, not " + std::to_string(source.pixel_reads()));
}

// Two regions in opposite corners of a 4 x 4 grid of chunks one channel deep, under a budget that
// affords everything. They occupy two chunks of the sixteen, in two chunk rows, so they are two
// footprints of one chunk each, and the four channels are one block that each footprint takes in a
// single read.
//
// The block is handed over once between the two reads, and what it says is a fraction of the two
// chunks the regions occupy -- four chunks read of eight -- not of the sixteen the plane holds, which
// would be four of sixty-four.
void TestTwoFootprintsHandOverAsAFractionOfWhatTheyOccupy() {
    const auto image = MakeImage(16, 16, 4);
    const auto geometry = MakeGeometry(4, 4, 1);
    SyntheticPixelSource source(image, geometry, Value);
    ReadOptions options;
    options.read_budget_bytes = kRoomyBudget;

    const auto handed =
        ReduceAndRecord(source, image, geometry, {{0, 0, 4, 4}, {12, 12, 4, 4}}, Range{0, 4, 1}, options);
    RequireHandOvers(handed, {{0, 4, false, 0.5}, {0, 4, true, 1.0}}, "two regions two footprints apart");
    Require(source.pixel_reads() == 2, "one read a footprint, not " + std::to_string(source.pixel_reads()));
}

// One of those regions alone is one footprint taken in one read, and a block taken in one read is
// handed over once, finished -- not once before its only read and again after it.
void TestOneFootprintInOneReadHandsOverOnce() {
    const auto image = MakeImage(16, 16, 4);
    const auto geometry = MakeGeometry(4, 4, 1);
    SyntheticPixelSource source(image, geometry, Value);
    ReadOptions options;
    options.read_budget_bytes = kRoomyBudget;

    const auto handed = ReduceAndRecord(source, image, geometry, {{0, 0, 4, 4}}, Range{0, 4, 1}, options);
    RequireHandOvers(handed, {{0, 4, true, 1.0}}, "one region in one chunk");
    Require(source.pixel_reads() == 1, "one read, not " + std::to_string(source.pixel_reads()));
}

std::vector<double> CubeProgress(const SyntheticPixelSource& source, const ImageDescriptor& image,
                                 const ChunkGeometry& geometry, const Range& spectral, std::uint64_t sample) {
    WorkPool workers(1);
    const auto reducible = Reducible(source, image, geometry, workers);
    carta::zarr::CubeHistogramRequest request;
    request.planes.spectral = spectral;
    request.bins = 16;
    request.spatial_sample = sample;
    ReadOptions options;
    options.read_budget_bytes = 1;
    std::vector<double> reported;
    const auto outcome = carta::zarr::internal::ComputeCubeHistogram(
        reducible, request, options, [&](const carta::zarr::CubeHistogramProgress& update) {
            reported.push_back(update.progress);
            return true;
        });
    Require(static_cast<bool>(outcome), "the cube histogram failed");
    return reported;
}

std::string Describe(const std::vector<double>& reported) {
    std::string text;
    for (const double one : reported) {
        text += " " + std::to_string(one);
    }
    return text;
}

// A cube histogram is one run rather than blocks, and reports before every read after its first
// and never after its last.
//
// The 8 x 8 x 8 image above in one-chunk reads is eight reads over a run of eight chunks -- four to
// a layer, two layers deep -- so it reports seven times, at one through seven eighths.
//
// Sampled every eighth pixel, a 16 x 16 plane of 4 x 4 chunks keeps pixels only in the chunks whose
// first row and column are multiples of eight: four of sixteen. Those four are the run, and the
// chunks the sample steps over are neither read nor counted, so it reports at one, two and three
// quarters. It used to count the stepped-over chunks as it passed them, which made the fraction jump
// over what it did not read -- 2/16, 8/16, 10/16 -- and end three quarters short of the whole.
void TestACubeHistogramReportsEveryReadButItsFirst() {
    {
        const auto image = MakeImage(8, 8, 8);
        const auto geometry = MakeGeometry(4, 4, 4);
        SyntheticPixelSource source(image, geometry, Value);
        const auto reported = CubeProgress(source, image, geometry, Range{0, 8, 1}, 1);
        const std::vector<double> expected{1.0 / 8, 2.0 / 8, 3.0 / 8, 4.0 / 8, 5.0 / 8, 6.0 / 8, 7.0 / 8};
        Require(reported == expected,
                "a whole cube reported" + Describe(reported) + " rather than" + Describe(expected));
        Require(source.pixel_reads() == 8, "eight one-chunk reads, not " + std::to_string(source.pixel_reads()));
    }
    {
        const auto image = MakeImage(16, 16, 4);
        const auto geometry = MakeGeometry(4, 4, 4);
        SyntheticPixelSource source(image, geometry, Value);
        const auto reported = CubeProgress(source, image, geometry, Range{0, 4, 1}, 8);
        const std::vector<double> expected{1.0 / 4, 2.0 / 4, 3.0 / 4};
        Require(reported == expected,
                "a sampled cube reported" + Describe(reported) + " rather than" + Describe(expected));
        Require(source.pixel_reads() == 4, "four reads, not " + std::to_string(source.pixel_reads()));
    }
}

}  // namespace

int main() {
    try {
        TestAHistogramCountsEveryPixel();
        TestARangeWiderThanAFloatStillBinsItsPixels();
        TestEveryRangeIsBinnedAsTheCallerBinsIt();
        TestASpectralReductionAgreesWithTheFormula();
        TestAMaskedRegionReadsOnlyTheChunksItOccupies();
        TestAPlaneHistogramSplitAcrossWorkers();
        TestACubeHistogramSplitAcrossWorkers();
        TestASpectralSpreadIsCentredWhateverTheMagnitude();
        TestACubeSpreadIsCentredWhateverTheMagnitude();
        TestProgressNeverClaimsTheWholeRunBeforeItsLastRead();
        TestAPlaneHistogramHandsOverAtEveryReadOfItsBlock();
        TestARegionCoveringThePlaneHandsOverAsThePlaneDoes();
        TestTwoFootprintsHandOverAsAFractionOfWhatTheyOccupy();
        TestOneFootprintInOneReadHandsOverOnce();
        TestACubeHistogramReportsEveryReadButItsFirst();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "reduce synthetic test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
