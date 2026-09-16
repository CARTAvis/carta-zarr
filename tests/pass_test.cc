/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a pass decides before it reads anything.
//
// This is the half of the pass with an answer worth checking on its own: how deep a slab goes, how
// wide a band is, which spatial axis is the inner one, and how much a read may decode. Every one of
// those has moved at least once, and until now each could only be observed through a reduction's
// output against a directory tree on disk -- so the tests that cover them assert how many times a
// result was handed over and reason backwards to what the walk must have read.
//
// None of it needs a store, a transport or a fixture, which is the point.

#include "reduce/axis_map.h"
#include "reduce/pass.h"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ChunkGeometry;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadOptions;
using carta::zarr::internal::MapAxes;
using carta::zarr::internal::PassPlan;
using carta::zarr::internal::PlanPass;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// An image in the logical order XRADIO reports: x, y, spectral, polarization, time.
ImageDescriptor MakeImage(std::uint64_t x, std::uint64_t y, std::uint64_t channels, bool has_mask = false) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    descriptor.has_pixel_mask = has_mask;
    descriptor.pixel_mask_id = has_mask ? "FLAG" : "";
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, x},
             {"m", AxisRole::spatial_y, y},
             {"frequency", AxisRole::spectral, channels},
             {"polarization", AxisRole::polarization, 1},
             {"time", AxisRole::time, 1}};
    for (std::size_t i = 0; i < 5; ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = axes[i].name;
        axis.role = axes[i].role;
        axis.length = axes[i].length;
        // The stored order XRADIO writes: time, frequency, polarization, l, m -- so m is last and
        // varies fastest.
        axis.storage_index = std::vector<std::size_t>{3, 4, 1, 2, 0}.at(i);
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

ChunkGeometry MakeGeometry(std::uint64_t chunk_x, std::uint64_t chunk_y, std::uint64_t chunk_z,
                           AxisRole fastest) {
    ChunkGeometry geometry;
    geometry.fastest_spatial_axis = fastest;
    geometry.chunk_shape = {chunk_x, chunk_y, chunk_z, 1, 1};
    return geometry;
}

PassPlan Plan(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const Range& spectral,
              const ReadOptions& options, std::uint64_t sample = 1) {
    const auto map = MapAxes(descriptor);
    Require(static_cast<bool>(map), "MapAxes failed on a well-formed image");
    return PlanPass(descriptor, geometry, map.value(), spectral, 0, 0, sample, options);
}

// Which spatial axis the pass walks along is the store's decision, not the image's. Reading a plane
// with the other one fastest transposes every chunk on the way into the destination.
void TestTheInnerAxisFollowsTheStore() {
    const auto image = MakeImage(512, 520, 32);
    const Range spectral{0, 32, 1};
    const ReadOptions options;

    const auto m_fastest = Plan(image, MakeGeometry(256, 260, 2, AxisRole::spatial_y), spectral, options);
    Require(m_fastest.axis_u == 1 && m_fastest.axis_v == 0,
            "when the store varies y fastest, y is the pass's inner axis");
    Require(m_fastest.u_length == 520 && m_fastest.v_length == 512, "the lengths follow the axes");
    Require(m_fastest.chunk_u == 260 && m_fastest.chunk_v == 256, "so do the chunk extents");

    const auto l_fastest = Plan(image, MakeGeometry(256, 260, 2, AxisRole::spatial_x), spectral, options);
    Require(l_fastest.axis_u == 0 && l_fastest.axis_v == 1,
            "when the store varies x fastest, x is the pass's inner axis");
    Require(l_fastest.u_length == 512 && l_fastest.v_length == 520, "the lengths follow the axes");
}

// A budget is in decoded chunk bytes, and a masked read decodes the flag too -- at the flag's own
// cost, which is the thing this used to get wrong. RequireUsableFlag holds a flag to bool over the
// image's shape, so beside a float32 chunk it is a quarter of one and not a second one. Counting it
// as a second one shrank every masked read by the difference.
void TestAMaskCostsWhatTheFlagCosts() {
    const Range spectral{0, 32, 1};
    const auto geometry = MakeGeometry(256, 260, 2, AxisRole::spatial_y);
    const std::uint64_t chunk_elements = 256ULL * 260ULL * 2ULL;
    ReadOptions options;
    options.apply_pixel_mask = true;

    const auto plain = Plan(MakeImage(512, 520, 32, false), geometry, spectral, options);
    const auto masked = Plan(MakeImage(512, 520, 32, true), geometry, spectral, options);
    Require(masked.apply_mask && !plain.apply_mask, "only an image with a flag applies one");
    Require(plain.chunk_bytes == chunk_elements * 4, "a float32 chunk is four bytes an element");
    Require(masked.chunk_bytes == plain.chunk_bytes + chunk_elements,
            "the flag beside it is one byte an element, so a masked float32 chunk costs a quarter "
            "more and not twice as much");

    ReadOptions declined;
    declined.apply_pixel_mask = false;
    const auto not_applied = Plan(MakeImage(512, 520, 32, true), geometry, spectral, declined);
    Require(!not_applied.apply_mask && not_applied.chunk_bytes == plain.chunk_bytes,
            "an image with a flag the caller declined costs what an unmasked one costs");
}

// The caller's ceiling is taken as given; without one the pass asks chunk_blocks for the policy.
void TestTheCallersCeilingWins() {
    const auto image = MakeImage(512, 520, 32);
    const auto geometry = MakeGeometry(256, 260, 2, AxisRole::spatial_y);
    const Range spectral{0, 32, 1};

    ReadOptions limited;
    limited.temporary_memory_limit_bytes = 1u << 20;
    const auto small = Plan(image, geometry, spectral, limited);
    Require(small.slab_budget_bytes == (1u << 20), "a stated limit is the budget");

    const auto def = Plan(image, geometry, spectral, ReadOptions{});
    Require(def.slab_budget_bytes ==
                carta::zarr::internal::DefaultReadBytes(
                    carta::zarr::internal::DecodedChunkBytes(image, geometry)),
            "without a limit the budget is the one chunk_blocks measured");
}

// least_channels is how many selected channels one chunk of the spectral axis holds, which is what
// keeps a slab from ending inside a chunk and making one decode serve two slabs.
void TestASlabIsCountedInChunksOfTheSpectralAxis() {
    const auto image = MakeImage(512, 520, 64);
    const auto geometry = MakeGeometry(256, 260, 8, AxisRole::spatial_y);

    Require(Plan(image, geometry, Range{0, 64, 1}, ReadOptions{}).least_channels == 8,
            "eight channels to a chunk, read every one");
    Require(Plan(image, geometry, Range{0, 32, 2}, ReadOptions{}).least_channels == 4,
            "eight channels to a chunk, every second one selected, is four");
    Require(Plan(image, geometry, Range{0, 8, 8}, ReadOptions{}).least_channels == 1,
            "a stride of a whole chunk selects one channel from each");
    Require(Plan(image, geometry, Range{0, 4, 16}, ReadOptions{}).least_channels == 1,
            "a stride wider than a chunk still selects one, never none");
}

// layer_chunks is the chunks in one spectral layer, and band_rows is how many chunk rows of it one
// read may hold. Both feed the progress a caller sees, so a wrong one is a bar that lies.
void TestALayerIsCountedInWholeChunks() {
    const auto image = MakeImage(512, 520, 32);
    // Two chunks along each spatial axis.
    const auto geometry = MakeGeometry(256, 260, 2, AxisRole::spatial_y);
    const auto plan = Plan(image, geometry, Range{0, 32, 1}, ReadOptions{});
    Require(plan.layer_chunks == 4, "two chunks each way is four to a layer");

    // An axis that does not divide by its chunk still counts the partial chunk.
    const auto ragged = Plan(MakeImage(513, 520, 32), geometry, Range{0, 32, 1}, ReadOptions{});
    Require(ragged.layer_chunks == 6, "513 pixels over a 256-wide chunk is three, not two");
}

// A budget below one chunk row does not produce a band of zero rows, which would read nothing and
// never advance.
void TestABandIsNeverEmpty() {
    const auto image = MakeImage(512, 520, 32);
    const auto geometry = MakeGeometry(256, 260, 2, AxisRole::spatial_y);
    ReadOptions tiny;
    tiny.temporary_memory_limit_bytes = 1;
    const auto plan = Plan(image, geometry, Range{0, 32, 1}, tiny);
    Require(plan.band_rows >= 1, "a band holds at least one chunk row however small the budget");
    Require(plan.least_channels >= 1, "and a slab at least one channel");
}

// Sampling of zero would select nothing and divide by nothing; the plan floors it at one.
void TestSamplingFloorsAtOne() {
    const auto image = MakeImage(512, 520, 32);
    const auto geometry = MakeGeometry(256, 260, 2, AxisRole::spatial_y);
    Require(Plan(image, geometry, Range{0, 32, 1}, ReadOptions{}, 0).sample == 1,
            "a sample of zero reads every pixel rather than none");
    Require(Plan(image, geometry, Range{0, 32, 1}, ReadOptions{}, 4).sample == 4,
            "a sample the caller stated is kept");
}

// The sampled range is the one the pass asks the store for. Its edges are where an off-by-one costs
// a row of the image.
void TestSampledRangePicksTheMultiplesInside() {
    std::uint64_t start = 0;
    std::uint64_t count = 0;
    carta::zarr::internal::SampledRange(0, 10, 1, start, count);
    Require(start == 0 && count == 10, "a stride of one selects the whole range");

    carta::zarr::internal::SampledRange(0, 10, 3, start, count);
    Require(start == 0 && count == 4, "0, 3, 6 and 9 fall inside [0, 10)");

    carta::zarr::internal::SampledRange(4, 10, 3, start, count);
    Require(start == 6 && count == 2, "the first multiple of three at or after 4 is 6, then 9");

    carta::zarr::internal::SampledRange(7, 9, 3, start, count);
    Require(count == 0, "a range holding no multiple selects nothing rather than one");
}

}  // namespace

int main() {
    try {
        TestTheInnerAxisFollowsTheStore();
        TestAMaskCostsWhatTheFlagCosts();
        TestTheCallersCeilingWins();
        TestASlabIsCountedInChunksOfTheSpectralAxis();
        TestALayerIsCountedInWholeChunks();
        TestABandIsNeverEmpty();
        TestSamplingFloorsAtOne();
        TestSampledRangePicksTheMultiplesInside();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "pass test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
