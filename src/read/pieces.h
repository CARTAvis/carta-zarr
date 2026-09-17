/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_READ_PIECES_H_
#define CARTA_ZARR_SRC_READ_PIECES_H_

#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include "pixel_source.h"

#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

/**
 * How one ordinary read is cut into pieces. See CONTEXT.md for what a piece is.
 *
 * A read that is not split is described here as one piece covering everything, so that the loop
 * reading it is the same loop either way and there are not two paths to keep in agreement.
 */
struct PiecePlan {
    // Whether the read is cut at all. False means one piece, and then `axis` and `chunk` say
    // nothing.
    bool split = false;
    std::size_t axis = 0;
    // Units of `axis` the read covers, and how many destination elements one such unit is worth.
    // Unsplit, that is one unit worth the whole destination.
    std::uint64_t units = 1;
    std::uint64_t elements_per_unit = 0;
    // Units one piece should cover, and the chunk extent along `axis` that the end of a piece is
    // rounded out to.
    std::uint64_t units_per_piece = 1;
    std::uint64_t chunk = 0;
};

/**
 * Decide how to cut a read.
 *
 * Pure: it reaches no further than the descriptor, the geometry and what the caller asked for, so
 * the strategy is checkable without a store, a transport or a directory tree.
 *
 * A read is cut when there is a reason to cut it, and either reason is enough on its own. Somebody
 * to report progress to is one. A stated memory ceiling is the other: it says how much the read may
 * hold at once, and splitting to fit is a better answer than refusing to read at all. A read with
 * neither reason, or with no axis selecting more than one element, is one piece.
 *
 * `elements` is what the whole read will produce, and `apply_mask` says whether a flag is decoded
 * beside the pixels -- both halves of the sizing count it, because a budget divided by a per-piece
 * cost that ignores the flag would size pieces against a cost the read does not have.
 */
PiecePlan PlanPieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                     const ReadRequest& request, const ReadOptions& options, std::uint64_t elements,
                     bool apply_mask);

/**
 * Read a densely packed float32 result, one piece at a time.
 *
 * Everything an ordinary read does apart from being reached through a handle: it validates the
 * request against the descriptor, plans the pieces, reads each one's flag and pixels in that order,
 * folds the flag in, and reports progress. The caller supplies the source and translates whatever
 * comes back; it does not need to know that any of this happened.
 *
 * The flag is read before the pixels, which is the opposite of what a pass does and is the reason
 * ADR 0005 gives for this staying outside one: the destination is the caller's, so a mask that
 * cannot be read must not leave a piece of it updated.
 *
 * Reports buffer_too_small when the destination cannot hold the selection, or when a piece's flag
 * buffer exceeds a ceiling that no further splitting gets under.
 */
Result<std::size_t> ReadInPieces(const PixelSource& source, const ImageDescriptor& descriptor,
                                 const ChunkGeometry& geometry, const ReadRequest& request,
                                 BufferView<float> destination, const ReadOptions& options);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_READ_PIECES_H_
