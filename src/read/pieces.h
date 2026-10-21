/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
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

// The most chunks in all that a piece of a read with a progress callback decodes, unless one chunk
// along the cut already spans more. Such a read reports each piece as it lands, so the piece size is how long its
// caller waits between updates, as with the cursor spectrum carta-backend draws as it fills. A budget based on the
// decode threads grows with them: cut to that alone, a spectrum of the 512 x 512 x 30000 test cube came in nine pieces
// on 28 threads, 160 ms apart on a cold local NVMe disk. The wait grows with what each chunk costs to fetch, not to
// decode, since carta-backend keeps two I/O threads. Cut to sixteen chunks, it came in thirty pieces 55 ms apart, and
// took a tenth longer in total.
inline constexpr std::uint64_t kWatchedPieceChunks = 16;

/**
 * One piece of an ordinary read. See GLOSSARY.md for what a piece is.
 *
 * What to read (the caller's request, narrowed along the axis the read is cut on) and where in the
 * destination it lands. A piece fills the destination from `first_element` for as many elements as
 * its request selects, and the pieces of a read fill it end to end, so the finished part is always a
 * prefix.
 *
 * A piece is read straight into the destination unless even the smallest possible piece (one chunk
 * along the cut, and the whole of every other axis) holds more than the budget allows. Then
 * `segments` are the parts it is read in, cut along the other axes, each read straight into its place
 * in the piece's part of the destination. A segment does not fill a contiguous run of the
 * destination, so it is not a piece: progress is reported when the whole piece is in, and the prefix
 * holds.
 */
struct Piece {
    ReadRequest request;
    std::uint64_t first_element = 0;
    std::vector<ReadRequest> segments;
};

/**
 * Cut a read into pieces.
 *
 * Pure: it uses only the descriptor, the geometry and the request, so it can be tested without a
 * store, a transport or a directory tree.
 *
 * Every read is cut to fit its budget (the caller's, or the library's own when the caller gives none),
 * so what a read holds is bounded whether or not it has a progress callback. A read with one
 * (`reports_progress` true) is also cut to at most kWatchedPieceChunks, so that progress comes often.
 * The one exception is a single plane with no progress callback whose budget allows a chunk for each
 * of `decode_threads`: it is read whole, which is faster than reading segments one after another, and
 * the plane bounds what it holds beyond the budget. A read that fits, or that has no axis selecting
 * more than one element, is one piece covering everything, so the same loop reads it either way.
 *
 * This decides where to cut, how much one piece may cover, where its end is rounded out to a chunk
 * boundary, and the segments of a piece too large to read whole. Whether the flag is decoded with the
 * pixels, which the sizing counts, comes from AppliesPixelMask; what it costs comes from the flag's
 * own layout, `flag_geometry`, which may differ from the pixels'.
 *
 * The request has been checked against the descriptor already: it is a selection of this image.
 */
std::vector<Piece> PlanPieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                              const ChunkGeometry& flag_geometry, const ReadRequest& request,
                              const ReadOptions& options, std::size_t decode_threads, bool reports_progress);

/**
 * Read a densely packed float32 result, one piece at a time.
 *
 * Everything an ordinary read does apart from the handle: it validates the request against the
 * descriptor, plans the pieces, reads each piece's flag and then its pixels, applies the flag, and
 * reports progress. The caller supplies the source and translates the result.
 *
 * It takes the source, the descriptor and the two geometries, which are all it uses, rather than the
 * ReducibleImage a reduction gets, which also carries an axis map and a pool a read does not need.
 *
 * The flag is read before the pixels, the opposite of a pass, which is the reason ADR 0005 gives for
 * keeping reads outside the pass: the destination is the caller's, so a mask that cannot be read must
 * not leave part of a piece updated. A piece read in segments does the same for each segment; the
 * segments already in place stay there, past the finished prefix.
 *
 * Reports buffer_too_small when the destination cannot hold the selection.
 */
Result<std::size_t> ReadInPieces(const PixelSource& source, const ImageDescriptor& descriptor,
                                 const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                                 const ReadRequest& request, BufferView<float> destination, const ReadOptions& options,
                                 std::size_t decode_threads, const ProgressCallback& progress);

/**
 * One selected element in every chunk that `request` touches, and no more.
 *
 * Reading a single element of a chunk decodes all of it, so this request decodes exactly the chunks
 * a read of `request` would, into whatever cache the read goes through, with a destination of one
 * float per chunk rather than per pixel. A plane of a 7763 x 4742 image in 512 x 512 chunks is 160
 * floats, where reading the plane itself means allocating and filling 147 MB.
 *
 * Along an axis, the chunks a request touches are consecutive unless its stride is a chunk or more.
 * Then every element it selects is already in a chunk of its own, and the range is kept. Otherwise
 * the range becomes the first element of each chunk, from the chunk its start is in to the chunk its
 * last element is in. That element is in the chunk, though the request may not have selected it.
 *
 * Pure, and the request has been checked against the descriptor already.
 */
ReadRequest OneElementPerChunk(const ChunkGeometry& geometry, const ReadRequest& request);

/**
 * Decode the chunks a read of `request` would decode, keeping them in options.control's cache pool,
 * and return only how many there were.
 *
 * A read of OneElementPerChunk(request) through ReadInPieces, so it is checked and cancelled exactly
 * as a read is. With the mask applied, the flag's chunks are decoded too, since the read that follows
 * will want them. They are sampled by `flag_geometry` rather than the pixels' geometry: a flag
 * chunked finer than its image, sampled where the pixels are, would leave half its chunks unread.
 * The count includes both.
 */
Result<std::uint64_t> PrefetchChunks(const PixelSource& source, const ImageDescriptor& descriptor,
                                     const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                                     const ReadRequest& request, const ReadOptions& options,
                                     std::size_t decode_threads);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_READ_PIECES_H_
