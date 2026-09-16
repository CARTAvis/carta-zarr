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

namespace carta::zarr::internal {

Result<Slab> ReadSlab(const SlabSource& source, const PassPlan& plan, const ReadOptions& options,
                      const SlabRequest& request, SlabBuffers& buffers) {
    const auto& descriptor = *plan.descriptor;
    const auto rank = descriptor.axes.size();
    const Range spectral = plan.spectral;

    ReadRequest read_request;
    read_request.axes.assign(rank, Range{0, 1, 1});
    read_request.axes.at(plan.axis_u) = Range{request.u_start, request.u_count, request.u_stride};
    read_request.axes.at(plan.axis_v) = Range{request.v_start, request.v_count, request.v_stride};
    read_request.axes.at(plan.map.spectral) = Range{
        spectral.start + (request.channel_index * spectral.stride), request.channel_count, spectral.stride};
    if (plan.map.has_polarization) {
        read_request.axes.at(plan.map.polarization) = Range{plan.polarization, 1, 1};
    }
    if (plan.map.has_time) {
        read_request.axes.at(plan.map.time) = Range{plan.time, 1, 1};
    }

    auto selection = zarr::BuildSelection(descriptor, read_request);
    if (!selection) {
        return selection.error();
    }
    // Take the plane in the order the store wrote it. The reader's destination has its own
    // dimension 0 fastest, so asking for the stored dimensions reversed is asking for no transpose
    // at all: the last stored dimension, the one the array is contiguous along, lands fastest.
    // Read hands back logical order because its callers want a densely packed image; a pass wants
    // whatever is cheapest to read, and pays for the difference in nothing but these strides.
    for (std::size_t i = 0; i < rank; ++i) {
        selection.value().logical_to_stored.at(i) = rank - 1 - i;
    }

    // Strides of that destination, by stored dimension: the last one steps by 1 and each earlier
    // one by the product of those after it.
    std::vector<std::uint64_t> stored_stride(rank, 1);
    std::uint64_t running = 1;
    for (std::size_t stored = rank; stored-- > 0;) {
        stored_stride.at(stored) = running;
        running *= selection.value().count.at(stored);
    }
    const auto elements = static_cast<std::size_t>(running);

    buffers.pixels.resize(elements);
    if (auto read = source.ReadPixels(selection.value(), buffers.pixels.data(), buffers.pixels.size(), options);
        !read) {
        return read.error();
    }
    if (plan.apply_mask) {
        buffers.mask.resize(elements);
        if (auto read = source.ReadMask(selection.value(), buffers.mask.data(), buffers.mask.size(), options);
            !read) {
            return read.error();
        }
        // Fold the flag into the pixels rather than carry it into the inner loop: a flagged pixel
        // and a NaN pixel mean the same thing to every statistic here, and this is the rule
        // Image::Read already applies.
        for (std::size_t i = 0; i < elements; ++i) {
            if (buffers.mask[i] == 0) {
                buffers.pixels[i] = std::numeric_limits<float>::quiet_NaN();
            }
        }
    }

    Slab slab;
    slab.first_channel = request.channel_index;
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
