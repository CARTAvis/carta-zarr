/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The half of an ordinary read that reads. See pieces.cc for why the two are apart.
//
// It names no Store: pixels arrive through the PixelSource seam, so everything under src/read/
// compiles without one, the same way src/reduce/ already did.

#include "read/pieces.h"

#include "chunk_blocks.h"
#include "pixel_mask.h"
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

// What a checked request comes to: where the pixels are, and how many of them there are.
struct CheckedRequest {
    zarr::PixelSelection selection;
    std::uint64_t elements = 0;
};

// What every pixel read of an image asks before it touches storage.
//
// Written out twice until this: once here and once in the facade, where the mask read built its own
// selection, sized its own buffer check and made its own control check on the way to the pixel seam.
// The three are one question -- can this image serve this request into this buffer, now -- and a
// read that answered two of them would be a read that had not checked.
Result<CheckedRequest> CheckRead(const ImageDescriptor& descriptor, const ReadRequest& request,
                                 std::size_t destination_size, const ReadOptions& options) {
    auto selection = zarr::BuildSelection(descriptor, request);
    if (!selection) {
        return selection.error();
    }
    const auto elements = zarr::SelectionElementCount(selection.value());
    if (elements == 0 || elements > destination_size) {
        return MakeError(ErrorCode::invalid_argument, "Destination buffer is too small for the request",
                         descriptor.id);
    }

    // Before allocating a mask or starting any storage work. A cancelled request must not consume
    // temporary memory, or open an array, just to discover that it cannot proceed.
    if (auto control = zarr::CheckReadControl(options, descriptor.id); !control) {
        return control.error();
    }
    return CheckedRequest{std::move(selection.value()), elements};
}

}  // namespace

Result<std::size_t> ReadInPieces(const ReadableImage& image, const ReadRequest& request,
                                 BufferView<float> destination, const ReadOptions& options) {
    const auto& descriptor = image.descriptor();
    const auto& geometry = image.geometry();
    const auto& source = image.source();
    const auto checked = CheckRead(descriptor, request, destination.size, options);
    if (!checked) {
        return checked.error();
    }
    // The whole read's element count is what a piece's share is measured against; the selection each
    // piece is actually read through is built per piece below.
    const auto elements = checked.value().elements;

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
            auto mask_read = source.ReadMask(piece_selection.value(), mask.data(), mask.size(), options);
            if (!mask_read) {
                return mask_read.error();
            }
        }
        auto read = source.ReadPixels(piece_selection.value(), piece_pixels, piece_elements, options);
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

Result<std::size_t> ReadPixelMask(const ReadableImage& image, const ReadRequest& request,
                                  BufferView<std::uint8_t> destination, const ReadOptions& options) {
    const auto& descriptor = image.descriptor();
    if (!descriptor.has_pixel_mask) {
        return MakeError(ErrorCode::not_found, "This image has no pixel mask", descriptor.id);
    }

    const auto checked = CheckRead(descriptor, request, destination.size, options);
    if (!checked) {
        return checked.error();
    }
    const auto elements = checked.value().elements;
    if (auto read = image.source().ReadMask(checked.value().selection, destination.data,
                                            static_cast<std::size_t>(elements), options);
        !read) {
        return read.error();
    }
    return static_cast<std::size_t>(elements);
}

}  // namespace carta::zarr::internal
