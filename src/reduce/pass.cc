/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "reduce/pass.h"

namespace carta::zarr::internal {

PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const AxisMap& map,
                  const Range& spectral, std::uint64_t polarization, std::uint64_t time,
                  std::uint64_t sample, const ReadOptions& options) {
    PassPlan plan;
    plan.descriptor = &descriptor;
    plan.map = map;
    // The spatial axis the store varies fastest is the one to ask for first: reading a plane with
    // the other one fastest means transposing every chunk on the way into the destination.
    const bool swap_spatial = geometry.fastest_spatial_axis == AxisRole::spatial_y;
    plan.axis_u = swap_spatial ? map.y : map.x;
    plan.axis_v = swap_spatial ? map.x : map.y;
    plan.u_length = descriptor.axes.at(plan.axis_u).length;
    plan.v_length = descriptor.axes.at(plan.axis_v).length;
    plan.chunk_u = std::max<std::uint64_t>(1, geometry.chunk_shape.at(plan.axis_u));
    plan.chunk_v = std::max<std::uint64_t>(1, geometry.chunk_shape.at(plan.axis_v));
    plan.chunk_depth = std::max<std::uint64_t>(1, geometry.chunk_shape.at(map.spectral));
    plan.apply_mask = options.apply_pixel_mask && descriptor.has_pixel_mask;
    plan.chunk_bytes = DecodedChunkBytes(descriptor, geometry) * (plan.apply_mask ? 2 : 1);
    plan.slab_budget_bytes = options.temporary_memory_limit_bytes != 0 ? options.temporary_memory_limit_bytes
                                                                      : DefaultReadBytes(plan.chunk_bytes);
    plan.least_channels = ((plan.chunk_depth + spectral.stride - 1) / spectral.stride);
    const std::uint64_t row_chunks = std::max<std::uint64_t>(1, ((plan.u_length - 1) / plan.chunk_u) + 1);
    const std::uint64_t column_chunks = std::max<std::uint64_t>(1, ((plan.v_length - 1) / plan.chunk_v) + 1);
    plan.layer_chunks = std::max<std::uint64_t>(1, row_chunks * column_chunks);
    // How many chunk rows one read may hold, so that a read is a budget's worth of chunk data.
    plan.band_rows = std::max<std::uint64_t>(
        1, plan.slab_budget_bytes / std::max<std::uint64_t>(1, row_chunks * plan.chunk_bytes));
    plan.spectral = spectral;
    plan.polarization = polarization;
    plan.time = time;
    plan.sample = std::max<std::uint64_t>(1, sample);
    return plan;
}

}  // namespace carta::zarr::internal
