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

#include "pixel_mask.h"
#include "zarr/pixel_selection.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace carta::zarr::internal {
namespace {

// What every pixel read of an image asks before it touches storage.
//
// Written out twice until this: once here and once in the facade, where the mask read built its own
// selection, sized its own buffer check and made its own control check on the way to the pixel seam.
// The three are one question -- can this image serve this request into this buffer, now -- and a
// read that answered two of them would be a read that had not checked.
Result<zarr::PixelSelection> CheckRead(const ImageDescriptor& descriptor, const ReadRequest& request,
                                       std::size_t destination_size, const ReadControl& control) {
    auto selection = zarr::BuildSelection(descriptor, request, zarr::DestinationOrder::logical);
    if (!selection) {
        return selection.error();
    }
    if (selection.value().elements() > destination_size) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is too small for the request",
                     descriptor.id};
    }

    // Before allocating a mask or starting any storage work. A cancelled request must not consume
    // temporary memory, or open an array, just to discover that it cannot proceed.
    if (auto allowed = zarr::CheckReadControl(control, descriptor.id); !allowed) {
        return allowed.error();
    }
    return selection;
}

}  // namespace

Result<std::size_t> ReadInPieces(const ReadableImage& image, const ReadRequest& request,
                                 BufferView<float> destination, const ReadOptions& options,
                                 const ProgressCallback& progress) {
    const auto& descriptor = image.descriptor();
    const auto& geometry = image.geometry();
    const auto& source = image.source();
    const auto checked = CheckRead(descriptor, request, destination.size, options.control);
    if (!checked) {
        return checked.error();
    }
    // The whole read's element count is what a piece's share is measured against; the selection each
    // piece is actually read through is built per piece below.
    const auto elements = checked.value().elements();

    const bool apply_mask = AppliesPixelMask(options, descriptor);
    const auto pieces = PlanPieces(descriptor, geometry, request, options, static_cast<bool>(progress));

    std::vector<std::uint8_t> mask;
    for (const auto& piece : pieces) {
        auto piece_selection = zarr::BuildSelection(descriptor, piece.request, zarr::DestinationOrder::logical);
        if (!piece_selection) {
            return piece_selection.error();
        }
        const auto piece_elements = static_cast<std::size_t>(piece_selection.value().elements());
        float* piece_pixels = destination.data + piece.first_element;

        if (apply_mask) {
            // The limit bounds a piece, and a read that is not split is one piece, so a request that
            // cannot be cut any further still has to say so rather than allocate.
            if (options.temporary_memory_limit_bytes != 0 &&
                piece_elements > options.temporary_memory_limit_bytes) {
                return Error{ErrorCode::buffer_too_small,
                             "Pixel mask temporary buffer exceeds the configured memory limit",
                             descriptor.id};
            }
            mask.assign(piece_elements, 0);
            // The mask is read first so that an unavailable or cancelled mask cannot leave this
            // piece of the destination updated. TensorStore still owns the pixel operation's
            // in-flight completion before it returns, so the destination remains valid for the next
            // read.
            auto mask_read = source.ReadMask(piece_selection.value(), mask.data(), mask.size(), options.control);
            if (!mask_read) {
                return mask_read.error();
            }
        }
        auto read = source.ReadPixels(piece_selection.value(), piece_pixels, piece_elements, options.control);
        if (!read) {
            return read.error();
        }
        if (apply_mask) {
            ApplyPixelMask(piece_pixels, mask.data(), piece_elements);
        }

        const auto finished = static_cast<std::size_t>(piece.first_element + piece_elements);
        if (progress && !progress(finished, static_cast<std::size_t>(elements))) {
            return Error{ErrorCode::cancelled, "The read was cancelled by its progress callback",
                         descriptor.id};
        }
    }
    return static_cast<std::size_t>(elements);
}

Result<std::size_t> ReadPixelMask(const ReadableImage& image, const ReadRequest& request,
                                  BufferView<std::uint8_t> destination, const ReadControl& control) {
    const auto& descriptor = image.descriptor();
    if (!descriptor.has_pixel_mask) {
        return Error{ErrorCode::not_found, "This image has no pixel mask", descriptor.id};
    }

    const auto checked = CheckRead(descriptor, request, destination.size, control);
    if (!checked) {
        return checked.error();
    }
    const auto elements = checked.value().elements();
    if (auto read = image.source().ReadMask(checked.value(), destination.data,
                                            static_cast<std::size_t>(elements), control);
        !read) {
        return read.error();
    }
    return static_cast<std::size_t>(elements);
}

}  // namespace carta::zarr::internal
