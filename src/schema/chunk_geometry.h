/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_CHUNK_GEOMETRY_H_
#define CARTA_ZARR_SRC_SCHEMA_CHUNK_GEOMETRY_H_

// The stored layout, in the order the library reads.
//
// A StorageLayout is what the array says about itself, in stored order. A ChunkGeometry is the same
// facts in logical order, which is the order a schema profile chose and the order every read and
// every pass works in. Turning one into the other is a rule about a descriptor and belongs beside
// the descriptor rather than in the facade that happened to need it first: no handle, no store, no
// context is involved, and the two cases worth checking -- a sharded array, and one whose stored
// order is not the logical one -- needed a directory tree to reach while it lived there.

#include "carta-zarr/descriptor.h"

#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

// Requires a layout of the descriptor's own rank: every storage_index must index its chunk_shape,
// and its shard_shape too when it is sharded. ParseArrayMetadata guarantees this for any array that
// parsed, so a mismatch is a caller that built a layout some other way -- `at` throws, and
// Dataset::OpenImage is wrapped in Guarded, so it reaches that caller as invalid_metadata.
//
// Said with `at` rather than filled in. Falling back to the axis length is what this used to do, and
// it turned "no layout here" into a geometry claiming one chunk covers the whole image -- a lie a
// consumer cannot tell from an image that really is stored that way.
inline ChunkGeometry BuildChunkGeometry(const ImageDescriptor& descriptor, const StorageLayout& layout) {
    ChunkGeometry geometry;
    geometry.sharded = layout.sharded;
    geometry.compressor = layout.compressor;

    const auto rank = descriptor.axes.size();
    geometry.chunk_shape.resize(rank);
    geometry.shard_shape.resize(rank);
    geometry.grid_shape.resize(rank);
    for (std::size_t logical = 0; logical < rank; ++logical) {
        const auto& axis = descriptor.axes.at(logical);
        const auto stored = axis.storage_index;
        const auto chunk = layout.chunk_shape.at(stored);
        const auto shard = layout.sharded ? layout.shard_shape.at(stored) : chunk;
        geometry.chunk_shape.at(logical) = chunk;
        geometry.shard_shape.at(logical) = shard;
        // Both extents are positive by the same guarantee, so neither zero this divides by nor the
        // shard of no extent the line above used to substitute for can arrive.
        geometry.grid_shape.at(logical) = (axis.length + chunk - 1) / chunk;
        if (stored != logical) {
            geometry.transpose_required = true;
        }
    }

    // The last stored dimension varies fastest, so of the two spatial axes the one with the larger
    // storage index is the one a plane is contiguous along.
    std::size_t x_stored = 0;
    std::size_t y_stored = 0;
    for (const auto& axis : descriptor.axes) {
        if (axis.role == AxisRole::spatial_x) {
            x_stored = axis.storage_index;
        } else if (axis.role == AxisRole::spatial_y) {
            y_stored = axis.storage_index;
        }
    }
    geometry.fastest_spatial_axis = y_stored > x_stored ? AxisRole::spatial_y : AxisRole::spatial_x;
    return geometry;
}

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_SCHEMA_CHUNK_GEOMETRY_H_
