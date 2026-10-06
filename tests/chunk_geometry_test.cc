/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The stored layout permuted into the order the library reads.
//
// This rule lived in the facade, where the only way to it was opening a dataset on disk -- so the
// two cases it exists for, a sharded array and one whose stored order is not the logical one, were
// checked by whichever committed fixture happened to have them. It takes a descriptor and a layout
// and answers, so this target links nothing.

#include "schema/chunk_geometry.h"

#include "support/check.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ImageDescriptor;
using carta::zarr::internal::BuildChunkGeometry;
using carta::zarr::internal::zarr::StorageLayout;

using carta::zarr::testing::Require;

// Five axes in logical order, with the stored index of each one given: that permutation is the
// whole question this answers.
ImageDescriptor MakeImage(const std::vector<std::size_t>& storage, const std::vector<std::uint64_t>& lengths) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    const AxisRole roles[]{AxisRole::spatial_x, AxisRole::spatial_y, AxisRole::spectral, AxisRole::polarization,
                           AxisRole::time};
    const char* names[]{"l", "m", "frequency", "polarization", "time"};
    for (std::size_t i = 0; i < storage.size(); ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = names[i];
        axis.role = roles[i];
        axis.length = lengths.at(i);
        axis.storage_index = storage.at(i);
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

// The order XRADIO writes: time, frequency, polarization, l, m -- so m is the last stored dimension
// and the one a plane is contiguous along, and every axis of the result has to be looked up through
// the permutation rather than taken positionally.
void TestTheLayoutIsReadThroughThePermutation() {
    const auto image = MakeImage({3, 4, 1, 2, 0}, {512, 520, 32, 4, 2});
    StorageLayout layout;
    // In stored order: time, frequency, polarization, l, m.
    layout.chunk_shape = {1, 2, 1, 256, 260};
    layout.compressor = "zstd";

    const auto geometry = BuildChunkGeometry(image, layout);

    Require(geometry.chunk_shape == std::vector<std::uint64_t>{256, 260, 2, 1, 1},
            "the chunk shape should come back in logical order");
    Require(geometry.grid_shape == std::vector<std::uint64_t>{2, 2, 16, 4, 2},
            "and the grid should be each axis's length over its own chunk");
    Require(geometry.fastest_spatial_axis == AxisRole::spatial_y, "m is stored last, so a plane is contiguous along m");
    Require(!geometry.sharded && geometry.shard_shape == geometry.chunk_shape,
            "an unsharded array's shard is its chunk");
    Require(geometry.compressor == "zstd", "and the compressor is carried through");
}

// The same facts when the store writes the logical order: nothing to permute.
void TestAnUntransposedLayout() {
    const auto image = MakeImage({0, 1, 2, 3, 4}, {64, 40, 6, 1, 1});
    StorageLayout layout;
    layout.chunk_shape = {16, 20, 2, 1, 1};

    const auto geometry = BuildChunkGeometry(image, layout);

    Require(geometry.chunk_shape == std::vector<std::uint64_t>{16, 20, 2, 1, 1}, "nothing to permute");
    Require(geometry.fastest_spatial_axis == AxisRole::spatial_y,
            "m is stored after l here, so a plane is still contiguous along m");
}

// The other branch of the one decision this makes that is not a permutation. Which spatial axis is
// fast is not a property of the image, it is which of the two the store wrote last -- the repository
// has a fixture of each -- and a reduction that gets it wrong transposes every chunk it reads.
void TestTheFastSpatialAxisIsWhicheverWasStoredLast() {
    StorageLayout layout;
    layout.chunk_shape = {1, 2, 1, 260, 256};

    const auto m_last = BuildChunkGeometry(MakeImage({3, 4, 1, 2, 0}, {512, 520, 32, 4, 2}), layout);
    Require(m_last.fastest_spatial_axis == AxisRole::spatial_y, "m stored last makes m the fast axis");

    const auto l_last = BuildChunkGeometry(MakeImage({4, 3, 1, 2, 0}, {512, 520, 32, 4, 2}), layout);
    Require(l_last.fastest_spatial_axis == AxisRole::spatial_x, "and l stored last makes it l");
}

// A chunk is what must be decoded to reach a byte; a shard is what one request fetches. They are
// separate granularities and a sharded array is the case where they differ, which is the reason the
// two fields exist rather than one.
void TestASharedArrayKeepsBothGranularities() {
    const auto image = MakeImage({3, 4, 1, 2, 0}, {512, 520, 32, 4, 2});
    StorageLayout layout;
    layout.sharded = true;
    layout.chunk_shape = {1, 1, 1, 64, 65};
    layout.shard_shape = {1, 4, 1, 256, 260};

    const auto geometry = BuildChunkGeometry(image, layout);

    Require(geometry.sharded, "a sharded layout is reported as one");
    Require(geometry.chunk_shape == std::vector<std::uint64_t>{64, 65, 1, 1, 1}, "the inner chunk in logical order");
    Require(geometry.shard_shape == std::vector<std::uint64_t>{256, 260, 4, 1, 1},
            "and the shard in logical order, permuted the same way");
    Require(geometry.grid_shape == std::vector<std::uint64_t>{8, 8, 32, 4, 2},
            "the grid counts inner chunks, not shards");
}

// An axis the layout does not describe is refused rather than filled in.
//
// A layout of the wrong rank is refused rather than filled in. Falling back to the axis length
// would turn "no layout for this image" into a geometry claiming one chunk covers the whole image,
// which a consumer could not tell from an image really stored that way. ParseArrayMetadata refuses
// the array that would produce it, so a layout of the wrong rank is a caller that built one some
// other way.
void TestAnAxisTheLayoutDoesNotDescribe() {
    const auto image = MakeImage({0, 1, 2, 3, 4}, {64, 40, 6, 1, 1});
    StorageLayout layout;
    // Only the two spatial dimensions are described.
    layout.chunk_shape = {16, 20};

    bool refused = false;
    try {
        BuildChunkGeometry(image, layout);
    } catch (const std::out_of_range&) {
        refused = true;
    }
    Require(refused, "a layout that does not describe every axis should throw, not be filled in");
}

}  // namespace

int main() {
    try {
        TestTheLayoutIsReadThroughThePermutation();
        TestAnUntransposedLayout();
        TestTheFastSpatialAxisIsWhicheverWasStoredLast();
        TestASharedArrayKeepsBothGranularities();
        TestAnAxisTheLayoutDoesNotDescribe();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "chunk geometry test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
