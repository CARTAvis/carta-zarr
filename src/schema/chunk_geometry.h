/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_CHUNK_GEOMETRY_H_
#define CARTA_ZARR_SRC_SCHEMA_CHUNK_GEOMETRY_H_

// The stored layout, in the order the library reads.
//
// A StorageLayout is what the array says about itself, in stored order. A ChunkGeometry is the same
// facts in logical order, the order a schema profile chose and every read and pass works in. The
// conversion needs only a descriptor, with no handle, store or context, so it lives here where it
// can be tested on its own: a sharded array, and one whose stored order is not the logical one.

#include "carta-zarr/descriptor.h"

#include "zarr/storage_layout.h"

#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

// Requires a layout of the descriptor's own rank: every storage_index must index its chunk_shape,
// and its shard_shape too when it is sharded. ParseArrayMetadata guarantees this for any array that
// parsed, so a mismatch is a caller that built a layout some other way -- `at` throws, and
// Dataset::OpenImage is wrapped in Guarded, so it reaches that caller as invalid_metadata.
//
// It uses `at` rather than filling in a default. Falling back to the axis length would report one
// chunk covering the whole image, which a consumer cannot tell from an image really stored that way.
inline ChunkGeometry BuildChunkGeometry(const ImageDescriptor& descriptor, const zarr::StorageLayout& layout) {
    ChunkGeometry geometry;
    geometry.sharded = layout.sharded;
    geometry.compressor = layout.compressor;

    const auto rank = descriptor.axes.size();
    geometry.chunk_shape.resize(rank);
    geometry.shard_shape.resize(rank);
    for (std::size_t logical = 0; logical < rank; ++logical) {
        const auto& axis = descriptor.axes.at(logical);
        const auto stored = axis.storage_index;
        const auto chunk = layout.chunk_shape.at(stored);
        const auto shard = layout.sharded ? layout.shard_shape.at(stored) : chunk;
        geometry.chunk_shape.at(logical) = chunk;
        geometry.shard_shape.at(logical) = shard;
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
