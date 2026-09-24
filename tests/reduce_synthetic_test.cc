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
#include "work_pool.h"

#include "support/check.h"

#include "support/synthetic_pixel_source.h"

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
using carta::zarr::internal::ReadableImage;
using carta::zarr::internal::WorkPool;
using carta::zarr::testing::SyntheticPixelSource;

using carta::zarr::testing::Require;

ReadableImage Readable(const SyntheticPixelSource& source, const ImageDescriptor& image,
                       const ChunkGeometry& geometry, WorkPool& workers) {
    auto readable = ReadableImage::Of(source, image, geometry, workers);
    Require(static_cast<bool>(readable), "the synthetic image's axes could not be mapped");
    return readable.value();
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
    options.temporary_memory_limit_bytes = 4 * 64 * 65 * 2 * 4;
    WorkPool workers(4);

    std::vector<std::uint64_t> counts(request.bins * kZ, 0);
    std::uint64_t blocks = 0;
    const auto readable = Readable(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ComputeHistogram(
        readable, request, [&](const carta::zarr::HistogramBlock& block) {
            if (!block.complete) {
                return true;
            }
            ++blocks;
            for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                for (std::size_t bin = 0; bin < block.bin_count; ++bin) {
                    counts.at(((block.first_channel + c) * request.bins) + bin) =
                        block.counts[(c * block.bin_count) + bin];
                }
            }
            return true;
        }, options);
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
        Require(counts.at(i) == expected.at(i),
                "bin " + std::to_string(i) + ": expected " + std::to_string(expected.at(i)) + ", got " +
                    std::to_string(counts.at(i)));
    }
    Require(source.most_hits_on_one_chunk() == 1, "and it decoded each chunk once getting there");
}

// Two overlapping regions over every channel, which is what the multi-region pass exists for.
void TestASpectralReductionAgreesWithTheFormula() {
    const auto image = MakeImage(kX, kY, kZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticPixelSource source(image, geometry, Value);

    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, kX, kY, nullptr},
        {40, 50, 100, 120, nullptr},
    };
    carta::zarr::SpectralReduceRequest request;
    request.planes.spectral = {0, kZ, 1};
    request.regions = regions.data();
    request.region_count = regions.size();
    request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::sum |
                         carta::zarr::Statistic::min | carta::zarr::Statistic::max;

    ReadOptions options;
    options.temporary_memory_limit_bytes = 4 * 64 * 65 * 2 * 4;
    WorkPool workers(4);

    std::vector<double> pixels(regions.size() * kZ, 0.0);
    std::vector<double> sums(regions.size() * kZ, 0.0);
    std::vector<double> minima(regions.size() * kZ, 0.0);
    std::vector<double> maxima(regions.size() * kZ, 0.0);
    const auto readable = Readable(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ReduceSpectral(
        readable, request, [&](const carta::zarr::SpectralBlock& block) {
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
        }, options);
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
        request.regions = regions.data();
        request.region_count = regions.size();
        request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::sum;

        ReadOptions options;
        options.temporary_memory_limit_bytes = kRoomyBudget;
        WorkPool workers(4);

        std::vector<double> sums(kChannels, 0.0);
        std::vector<double> counts(kChannels, 0.0);
        const auto readable = Readable(source, image, geometry, workers);
        const auto outcome = carta::zarr::internal::ReduceSpectral(
            readable, request, [&](const carta::zarr::SpectralBlock& block) {
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
            }, options);
        Require(static_cast<bool>(outcome), std::string(described) + ": the reduction failed: " +
                                                (outcome ? "" : outcome.error().message));

        for (std::uint64_t z = 0; z < kChannels; ++z) {
            Require(counts.at(z) == expected_count.at(z),
                    std::string(described) + ": the diagonal's pixel count");
            Require(std::abs(sums.at(z) - expected_sum.at(z)) <= 1e-9 * (1.0 + std::abs(expected_sum.at(z))),
                    std::string(described) + ": the diagonal's sum");
        }

        // The assertion the fixture-driven pair could only reach by inference. Four chunks of the
        // grid, one channel chunk, each decoded once.
        Require(source.chunks_touched() == kGrid,
                std::string(described) + ": read " + std::to_string(source.chunks_touched()) +
                    " chunks rather than the " + std::to_string(kGrid) + " the mask occupies");
        Require(source.most_hits_on_one_chunk() == 1,
                std::string(described) + ": a chunk was decoded more than once");
    };

    // This image varies y fastest, so the runs the reduction makes from the raster are ranges of
    // rows, one line per column -- the case that used to cost a column walk down the raster.
    carta::zarr::RegionMask rastered{0, 0, kSide, kSide, raster.data()};
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
    options.temporary_memory_limit_bytes = kRoomyBudget;

    const auto counts_from = [&](std::size_t worker_count) {
        WorkPool workers(worker_count);
        const auto readable = Readable(source, image, geometry, workers);
        std::vector<std::uint64_t> counts(request.bins * kSplitZ, 0);
        const auto outcome = carta::zarr::internal::ComputeHistogram(
            readable, request, [&](const carta::zarr::HistogramBlock& block) {
                if (!block.complete) {
                    return true;
                }
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    for (std::size_t bin = 0; bin < block.bin_count; ++bin) {
                        counts.at(((block.first_channel + c) * request.bins) + bin) =
                            block.counts[(c * block.bin_count) + bin];
                    }
                }
                return true;
            }, options);
        Require(static_cast<bool>(outcome),
                std::string("the histogram failed: ") + (outcome ? "" : outcome.error().message));
        return counts;
    };

    // One worker cannot split -- PlanRowTasks refuses before the plane is even measured -- so this
    // is the serial branch, and the oracle the other one has to match.
    const auto in_place = counts_from(1);
    const auto split = counts_from(4);
    for (std::size_t i = 0; i < in_place.size(); ++i) {
        Require(in_place.at(i) == split.at(i),
                "bin " + std::to_string(i) + ": binning in place gave " + std::to_string(in_place.at(i)) +
                    ", four workers gave " + std::to_string(split.at(i)));
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
    options.temporary_memory_limit_bytes = kRoomyBudget;
    WorkPool workers(4);

    const auto readable = Readable(source, image, geometry, workers);
    const auto outcome = carta::zarr::internal::ComputeCubeHistogram(readable, request, options);
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
    Require(result.num_pixels == pixels, "every pixel is finite, so every one counts");
    Require(result.nan_count == 0.0, "and none of them is a NaN");
    Require(result.minimum == smallest && result.maximum == largest,
            "the extremes are tracked exactly, whatever the split");
    Require(std::abs(result.sum - sum) <= 1e-9 * (1.0 + std::abs(sum)), "the sum over four accumulators");

    std::uint64_t binned = 0;
    for (const auto count : result.counts) {
        binned += count;
    }
    Require(static_cast<double>(binned) == pixels,
            "four provisional histograms re-aggregated onto one grid still hold every pixel");
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
    const auto readable = Readable(source, image, geometry, workers);

    ReadOptions options;
    options.temporary_memory_limit_bytes = 1;

    // Run to the end, and stopped short of it: a run whose ends both fall inside a chunk is the one
    // that covers a chunk more than its length suggests.
    const Range runs[]{{0, 24, 1}, {1, 23, 1}, {1, 8, 1}, {3, 6, 1},
                       {0, 8, 3}, {1, 8, 3}, {3, 7, 3}, {1, 4, 3}, {0, 3, 9}};
    std::string wrong;
    for (const auto& spectral : runs) {
        const auto start = spectral.start;
        const auto stride = spectral.stride;
        const auto count = spectral.count;
        const std::string which = " (start " + std::to_string(start) + ", count " +
                                  std::to_string(count) + ", stride " + std::to_string(stride) + ")";

        carta::zarr::CubeHistogramRequest cube;
        cube.planes.spectral = spectral;
        cube.bins = 16;
        std::vector<double> reported;
        cube.progress = [&](const carta::zarr::CubeHistogramProgress& update) {
            reported.push_back(update.progress);
            return true;
        };
        const auto histogram = carta::zarr::internal::ComputeCubeHistogram(readable, cube, options);
        Require(static_cast<bool>(histogram), "the cube histogram failed" + which);
        Require(!reported.empty(), "a one-byte budget should make the cube histogram report" + which);
        for (const double progress : reported) {
            if (progress >= 1.0) {
                wrong += "\n  a cube histogram reported " + std::to_string(progress) +
                         " with a read still to come" + which;
                break;
            }
        }

        // One block for the whole run, so that a block is several reads.
        const carta::zarr::RegionMask whole{0, 0, kSide, kSide, nullptr};
        carta::zarr::SpectralReduceRequest reduce;
        reduce.planes.spectral = spectral;
        reduce.regions = &whole;
        reduce.region_count = 1;
        reduce.statistics = static_cast<carta::zarr::StatisticSet>(carta::zarr::Statistic::sum);
        reduce.emit_every_channels = static_cast<std::uint32_t>(count);
        const auto reduced = carta::zarr::internal::ReduceSpectral(
            readable, reduce,
            [&](const carta::zarr::SpectralBlock& block) {
                if (!block.complete && block.completeness >= 1.0) {
                    wrong += "\n  an unfinished block claimed " + std::to_string(block.completeness) +
                             " of itself" + which;
                }
                return true;
            },
            options);
        Require(static_cast<bool>(reduced), "the spectral reduction failed" + which);
    }
    Require(wrong.empty(), "progress ran ahead of the reads:" + wrong);
}

}  // namespace

int main() {
    try {
        TestAHistogramCountsEveryPixel();
        TestASpectralReductionAgreesWithTheFormula();
        TestAMaskedRegionReadsOnlyTheChunksItOccupies();
        TestAPlaneHistogramSplitAcrossWorkers();
        TestACubeHistogramSplitAcrossWorkers();
        TestProgressNeverClaimsTheWholeRunBeforeItsLastRead();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "reduce synthetic test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
