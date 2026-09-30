/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_READ_PIECES_H_
#define CARTA_ZARR_SRC_READ_PIECES_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include "pixel_source.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace carta::zarr::internal {

/**
 * One piece of an ordinary read. See CONTEXT.md for what a piece is.
 *
 * What to read -- the caller's request, narrowed along the axis the read is cut on -- and where in the
 * destination it lands. A piece fills the destination from `first_element` for as many elements as
 * its request selects, and the pieces of a read fill it end to end, so the finished part is always a
 * prefix.
 */
struct Piece {
    ReadRequest request;
    std::uint64_t first_element = 0;
};

/**
 * Cut a read into pieces.
 *
 * Pure: it reaches no further than the descriptor, the geometry and what the caller asked for, so
 * the strategy is checkable without a store, a transport or a directory tree.
 *
 * A read is cut when there is a reason to cut it, and either reason is enough on its own. Somebody
 * to report progress to is one -- `watching` says whether there is, which is the whole of what this
 * ever asked about the callback. A stated memory ceiling is the other: it says how much the read may
 * hold at once, and splitting to fit is a better answer than refusing to read at all. A read with
 * neither reason, or with no axis selecting more than one element, is one piece covering everything,
 * so that the loop reading it is the same loop either way.
 *
 * Where to cut, how much one piece may cover, and where its end is rounded out to a chunk boundary
 * all happen here. They used to be handed to the one caller as a plan to carry out -- the cut axis,
 * the unit count, the elements per unit, the chunk extent -- and the caller did the arithmetic,
 * rewrote each piece's range and worked out where it landed. Whether the flag is decoded beside the
 * pixels, which the sizing counts, is asked of AppliesPixelMask rather than of the caller.
 *
 * The request has been checked against the descriptor already: it is a selection of this image.
 */
std::vector<Piece> PlanPieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                              const ReadRequest& request, const ReadOptions& options, bool watching);

/**
 * Read a densely packed float32 result, one piece at a time.
 *
 * Everything an ordinary read does apart from being reached through a handle: it validates the
 * request against the descriptor, plans the pieces, reads each one's flag and pixels in that order,
 * folds the flag in, and reports progress. The caller supplies the source and translates whatever
 * comes back; it does not need to know that any of this happened.
 *
 * It takes the source, the descriptor and the geometry, which are all it uses, rather than the
 * ReducibleImage a reduction is handed: that also carries an axis map and a pool, and a read has no
 * use for either.
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
                                 BufferView<float> destination, const ReadOptions& options,
                                 const ProgressCallback& progress);


}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_READ_PIECES_H_
