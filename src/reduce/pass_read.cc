/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The half of the pass that touches the reader.
//
// Separate from pass.cc so that what a pass decides stays linkable without it: PlanPass and the
// arithmetic beside it reach no further than the descriptor and the geometry, and the test that
// pins them compiles that translation unit alone, with no store, no transport and no TensorStore
// behind them.

#include "reduce/pass.h"

#include "pixel_mask.h"

namespace carta::zarr::internal {

Result<Slab> ReadSlab(const PixelSource& source, const PassPlan& plan, const ReadOptions& options,
                      const SlabRequest& request, SlabBuffers& buffers) {
    const auto& descriptor = *plan.descriptor;
    const auto rank = descriptor.axes.size();
    const Range spectral = plan.planes.spectral;

    ReadRequest read_request;
    read_request.axes.assign(rank, Range{0, 1, 1});
    read_request.axes.at(plan.axis_u) = Range{request.u_start, request.u_count, request.u_stride};
    read_request.axes.at(plan.axis_v) = Range{request.v_start, request.v_count, request.v_stride};
    read_request.axes.at(plan.map.spectral) = Range{
        spectral.start + (request.channel_index.index * spectral.stride), request.channel_count, spectral.stride};
    if (plan.map.has_polarization) {
        read_request.axes.at(plan.map.polarization) = Range{plan.planes.polarization, 1, 1};
    }
    if (plan.map.has_time) {
        read_request.axes.at(plan.map.time) = Range{plan.planes.time, 1, 1};
    }

    // Take the plane in the order the store wrote it, so that it is never transposed. Read hands
    // back logical order because its callers want a densely packed image; a pass wants whatever is
    // cheapest to read, and pays for the difference in nothing but these strides.
    auto selection = zarr::BuildSelection(descriptor, read_request, zarr::DestinationOrder::stored);
    if (!selection) {
        return selection.error();
    }
    const auto stored_stride = selection.value().DestinationStrides();
    const auto elements = static_cast<std::size_t>(selection.value().elements());

    buffers.pixels.resize(elements);
    if (auto read = source.ReadPixels(selection.value(), buffers.pixels.data(), buffers.pixels.size(), options.control);
        !read) {
        return read.error();
    }
    if (plan.apply_mask) {
        buffers.mask.resize(elements);
        if (auto read = source.ReadMask(selection.value(), buffers.mask.data(), buffers.mask.size(), options.control);
            !read) {
            return read.error();
        }
        // Folded into the pixels rather than carried into the inner loop, because a flagged pixel
        // and a NaN pixel mean the same thing to every statistic here.
        ApplyPixelMask(buffers.pixels.data(), buffers.mask.data(), elements);
    }

    Slab slab;
    slab.channel_count = request.channel_count;
    slab.pixels = buffers.pixels.data();
    slab.stride_u = stored_stride.at(descriptor.axes.at(plan.axis_u).storage_index);
    slab.stride_v = stored_stride.at(descriptor.axes.at(plan.axis_v).storage_index);
    slab.stride_z = stored_stride.at(descriptor.axes.at(plan.map.spectral).storage_index);
    slab.u_count = request.u_count;
    slab.v_count = request.v_count;
    return slab;
}

}  // namespace carta::zarr::internal
