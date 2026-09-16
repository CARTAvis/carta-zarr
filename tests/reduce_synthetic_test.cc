/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// A whole reduction, over pixels that were never written down.
//
// This target links no TensorStore and no Store. That is the point of it as much as the assertions
// are: until the pass took a SlabSource, a reduction could only be reached through a directory tree,
// so its arithmetic was pinned against the one committed fixture large enough to reach the parallel
// paths -- 512x520x4 -- and the oracle had to be written out beside it each time.
//
// Here the oracle is the formula the pixels come from, so the image can be whatever size the
// question needs.

#include "reduce/pass.h"
#include "reduce/plane_histogram.h"
#include "reduce/spectral_reduce.h"
#include "work_pool.h"

#include "support/synthetic_slab_source.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
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
using carta::zarr::internal::WorkPool;
using carta::zarr::testing::SyntheticSlabSource;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
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

// Every pixel of every plane binned, at a size no fixture reaches, against the formula rather than
// against a recorded answer.
void TestAHistogramCountsEveryPixel() {
    const auto image = MakeImage(kX, kY, kZ);
    const auto geometry = MakeGeometry(64, 65, 2);
    SyntheticSlabSource source(image, geometry, Value);

    carta::zarr::HistogramRequest request;
    request.spectral = {0, kZ, 1};
    request.bins = 50;
    request.lower = 0.0;
    request.upper = 1000.0;

    ReadOptions options;
    options.temporary_memory_limit_bytes = 4 * 64 * 65 * 2 * 4;
    WorkPool workers(4);

    std::vector<std::uint64_t> counts(request.bins * kZ, 0);
    std::uint64_t blocks = 0;
    const auto outcome = carta::zarr::internal::ComputeHistogram(
        source, image, geometry, request, [&](const carta::zarr::HistogramBlock& block) {
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
        }, options, workers);
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
    SyntheticSlabSource source(image, geometry, Value);

    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, kX, kY, nullptr},
        {40, 50, 100, 120, nullptr},
    };
    carta::zarr::SpectralReduceRequest request;
    request.spectral = {0, kZ, 1};
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
    const auto outcome = carta::zarr::internal::ReduceSpectral(
        source, image, geometry, request, [&](const carta::zarr::SpectralBlock& block) {
            if (!block.complete) {
                return true;
            }
            for (std::size_t r = 0; r < regions.size(); ++r) {
                for (std::size_t slot = 0; slot < block.statistic_count; ++slot) {
                    const double* from =
                        block.values + (r * block.region_stride) + (slot * block.statistic_stride);
                    auto* into = block.statistics[slot] == carta::zarr::Statistic::num_pixels ? &pixels
                                 : block.statistics[slot] == carta::zarr::Statistic::sum      ? &sums
                                 : block.statistics[slot] == carta::zarr::Statistic::min      ? &minima
                                                                                             : &maxima;
                    for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                        into->at((r * kZ) + block.first_channel + c) = from[c];
                    }
                }
            }
            return true;
        }, options, workers);
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

}  // namespace

int main() {
    try {
        TestAHistogramCountsEveryPixel();
        TestASpectralReductionAgreesWithTheFormula();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "reduce synthetic test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
