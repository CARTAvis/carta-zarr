/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The half of an ordinary read that touches the store. See pieces.cc for why the two are apart.

#include "read/pieces.h"

#include "chunk_blocks.h"
#include "pixel_mask.h"
#include "store.h"
#include "zarr/pixel_selection.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace carta::zarr::internal {
namespace {

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

}  // namespace

Result<std::size_t> ReadInPieces(const Store& store, const ImageDescriptor& descriptor,
                                 const ChunkGeometry& geometry, const ReadRequest& request,
                                 BufferView<float> destination, const ReadOptions& options) {
    auto selection = zarr::BuildSelection(descriptor, request);
    if (!selection) {
        return selection.error();
    }
    const auto elements = zarr::SelectionElementCount(selection.value());
    if (elements == 0 || elements > destination.size) {
        return MakeError(ErrorCode::invalid_argument, "Destination buffer is too small for the request",
                         descriptor.id);
    }

    // Do this before allocating a mask or starting any storage work. A cancelled request must not
    // consume temporary memory just to discover that it cannot proceed.
    auto control = zarr::CheckReadControl(options, descriptor.id);
    if (!control) {
        return control.error();
    }

    const bool apply_mask = options.apply_pixel_mask && descriptor.has_pixel_mask;
    const PiecePlan plan = PlanPieces(descriptor, geometry, request, options, elements, apply_mask);

    std::vector<std::uint8_t> mask;
    for (std::uint64_t begin = 0; begin < plan.units;) {
        const std::uint64_t end =
            plan.split ? AlignedBlockEnd(begin, plan.units_per_piece, plan.units,
                                         request.axes.at(plan.axis).start, request.axes.at(plan.axis).stride,
                                         plan.chunk)
                       : plan.units;

        ReadRequest piece = request;
        if (plan.split) {
            auto& range = piece.axes.at(plan.axis);
            range.start = request.axes.at(plan.axis).start + (begin * range.stride);
            range.count = end - begin;
        }
        auto piece_selection = zarr::BuildSelection(descriptor, piece);
        if (!piece_selection) {
            return piece_selection.error();
        }
        const auto piece_elements = static_cast<std::size_t>((end - begin) * plan.elements_per_unit);
        float* piece_pixels = destination.data + (begin * plan.elements_per_unit);

        if (apply_mask) {
            // The limit bounds a piece, and a read that is not split is one piece, so a request that
            // cannot be cut any further still has to say so rather than allocate.
            if (options.temporary_memory_limit_bytes != 0 &&
                piece_elements > options.temporary_memory_limit_bytes) {
                return MakeError(ErrorCode::buffer_too_small,
                                 "Pixel mask temporary buffer exceeds the configured memory limit",
                                 descriptor.id);
            }
            mask.assign(piece_elements, 0);
            // The mask is read first so that an unavailable or cancelled mask cannot leave this
            // piece of the destination updated. TensorStore still owns the pixel operation's
            // in-flight completion before it returns, so the destination remains valid for the next
            // read.
            auto mask_read = store.ReadPixelMaskBytes(descriptor.pixel_mask_id, piece_selection.value(),
                                                      mask.data(), mask.size(), options);
            if (!mask_read) {
                return mask_read.error();
            }
        }
        auto read = store.ReadPixelsFloat32(descriptor.id, piece_selection.value(), piece_pixels,
                                            piece_elements, options);
        if (!read) {
            return read.error();
        }
        if (apply_mask) {
            ApplyPixelMask(piece_pixels, mask.data(), piece_elements);
        }

        begin = end;
        if (options.progress && !options.progress(static_cast<std::size_t>(begin * plan.elements_per_unit),
                                                  static_cast<std::size_t>(elements))) {
            return MakeError(ErrorCode::cancelled, "The read was cancelled by its progress callback",
                             descriptor.id);
        }
    }
    return static_cast<std::size_t>(elements);
}

}  // namespace carta::zarr::internal
