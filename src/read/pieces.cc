/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a read decides before it reads anything.
//
// Separate from pieces_read.cc so that the strategy stays linkable without a store: PlanPieces and
// the arithmetic beside it reach no further than the descriptor and the geometry, and the test that
// pins them compiles this translation unit alone. The same division pass.cc and pass_read.cc make,
// and for the same reason.

#include "read/pieces.h"

#include "chunk_blocks.h"

#include <algorithm>
#include <optional>

namespace carta::zarr::internal {
namespace {

// The axis a split read is cut along: the slowest-varying one that selects more than a single
// element. The destination is dense in logical order with axis 0 fastest, so cutting there and
// nowhere else is what makes each finished piece extend a prefix instead of leaving holes.
std::optional<std::size_t> SlowestSelectedAxis(const ReadRequest& request) {
    for (std::size_t i = request.axes.size(); i-- > 0;) {
        if (request.axes.at(i).count > 1) {
            return i;
        }
    }
    return std::nullopt;
}

// How many elements of the split axis one piece should cover, so that the piece pulls roughly the
// budgeted amount of decompressed chunk data through. The other axes already contribute whatever
// they span, so a plane read needs far fewer rows per piece than a single-pixel column needs
// channels -- and an image with very large chunks gets pieces of one chunk rather than pieces it
// could never afford.
std::uint64_t UnitsPerPiece(const ImageDescriptor& descriptor, const ReadRequest& request,
                            const ChunkGeometry& geometry, std::size_t axis, std::size_t budget_bytes,
                            bool apply_mask) {
    std::uint64_t other_chunks = 1;
    for (std::size_t i = 0; i < request.axes.size(); ++i) {
        if (i == axis) {
            continue;
        }
        const auto chunk = i < geometry.chunk_shape.size() ? geometry.chunk_shape.at(i) : 0;
        const auto& range = request.axes.at(i);
        other_chunks *= ChunksSpanned(range.start, range.count, range.stride, chunk);
    }
    const auto row_bytes = DecodedChunkBytes(descriptor, geometry, apply_mask) * other_chunks;
    // At least one chunk: a piece smaller than that would decode the same chunk twice.
    const auto chunks = std::max<std::uint64_t>(1, budget_bytes / std::max<std::uint64_t>(1, row_bytes));
    const auto chunk = axis < geometry.chunk_shape.size() ? geometry.chunk_shape.at(axis) : 0;
    const auto stride = std::max<std::uint64_t>(1, request.axes.at(axis).stride);
    // AlignedBlockEnd rounds this out to a whole chunk, so a low estimate costs nothing.
    return std::max<std::uint64_t>(1, (chunks * std::max<std::uint64_t>(1, chunk)) / stride);
}

}  // namespace

PiecePlan PlanPieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                     const ReadRequest& request, const ReadOptions& options, bool watching,
                     std::uint64_t elements, bool apply_mask) {
    PiecePlan plan;
    plan.units = 1;
    plan.elements_per_unit = elements;

    const auto axis = SlowestSelectedAxis(request);
    if (!axis || (!watching && options.temporary_memory_limit_bytes == 0)) {
        return plan;
    }

    plan.split = true;
    plan.axis = *axis;
    plan.units = request.axes.at(*axis).count;
    plan.elements_per_unit = 1;
    for (std::size_t i = 0; i < *axis; ++i) {
        plan.elements_per_unit *= request.axes.at(i).count;
    }
    plan.chunk = *axis < geometry.chunk_shape.size() ? geometry.chunk_shape.at(*axis) : 0;

    // The flag is decoded beside the pixels when this read will apply it, so both halves of the
    // sizing count it: the budget the library chooses for itself, and the per-row cost that budget
    // is divided by. Counting it in one and not the other would size pieces against a cost the read
    // does not have.
    const auto chunk_bytes = DecodedChunkBytes(descriptor, geometry, apply_mask);
    const auto budget = options.temporary_memory_limit_bytes != 0 ? options.temporary_memory_limit_bytes
                                                                  : DefaultReadBytes(chunk_bytes);
    plan.units_per_piece = UnitsPerPiece(descriptor, request, geometry, *axis, budget, apply_mask);
    return plan;
}

}  // namespace carta::zarr::internal
