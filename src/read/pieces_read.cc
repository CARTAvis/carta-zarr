/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The half of an ordinary read that reads. See pieces.cc for why the two are apart.
//
// It names no Store: pixels arrive through the PixelSource seam, so everything under src/read/
// compiles without one, as src/reduce/ does.

#include "chunk_blocks.h"
#include "pixel_mask.h"
#include "read/pieces.h"
#include "zarr/pixel_selection.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace carta::zarr::internal {
namespace {

// What every pixel read of an image asks before it touches storage.
//
// The three checks are one question: can this image serve this request into this buffer, now? A read
// that answered only two of them would not have checked.
Result<zarr::PixelSelection> CheckRead(const ImageDescriptor& descriptor, const ReadRequest& request,
                                       std::size_t destination_size, const ReadControl& control) {
    auto selection = zarr::BuildSelection(descriptor, request, zarr::DestinationOrder::logical);
    if (!selection) {
        return selection.error();
    }
    if (selection.value().elements() > destination_size) {
        return Error{ErrorCode::buffer_too_small, "Destination buffer is too small for the request", descriptor.id};
    }

    // Before allocating a mask or starting any storage work. A cancelled request must not consume
    // temporary memory, or open an array, just to discover that it cannot proceed.
    if (auto allowed = zarr::CheckReadControl(control, descriptor.id); !allowed) {
        return allowed.error();
    }
    return selection;
}

// Applies `mask`, dense in the selection's own order, to the pixels it was read with, which lie in
// `destination` where the selection says: densely for a piece, and in its place in the piece for a
// segment. Axis 0 is the fastest in both, so the mask is applied one run of axis 0's count at a time.
void ApplyPixelMaskWithin(const zarr::PixelSelection& selection, float* destination, const std::uint8_t* mask) {
    const auto& order = selection.destination_to_stored();
    const auto strides = selection.DestinationStrides();
    const std::size_t rank = order.size();
    std::vector<std::uint64_t> count(rank);
    std::vector<std::uint64_t> stride(rank);
    for (std::size_t axis = 0; axis < rank; ++axis) {
        count[axis] = selection.count().at(order[axis]);
        stride[axis] = strides.at(order[axis]);
    }
    float* const base = destination + selection.destination_offset();
    std::vector<std::uint64_t> at(rank, 0);
    while (true) {
        std::uint64_t into = 0;
        for (std::size_t axis = 1; axis < rank; ++axis) {
            into += at[axis] * stride[axis];
        }
        float* run = base + into;
        for (std::uint64_t i = 0; i < count[0]; ++i) {
            if (mask[i] != 0) {
                run[i * stride[0]] = std::numeric_limits<float>::quiet_NaN();
            }
        }
        mask += count[0];
        // The next run, the faster axes first.
        std::size_t axis = 1;
        for (; axis < rank; ++axis) {
            if (++at[axis] < count[axis]) {
                break;
            }
            at[axis] = 0;
        }
        if (axis >= rank) {
            return;
        }
    }
}

}  // namespace

Result<std::size_t> ReadInPieces(const PixelSource& source, const ImageDescriptor& descriptor,
                                 const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                                 const ReadRequest& request, BufferView<float> destination, const ReadOptions& options,
                                 std::size_t decode_threads, const ProgressCallback& progress) {
    const auto checked = CheckRead(descriptor, request, destination.size, options.control);
    if (!checked) {
        return checked.error();
    }
    // A piece's share is measured against the whole read's element count; the selection each piece
    // is read through is built per piece below.
    const auto elements = checked.value().elements();

    const bool apply_mask = AppliesPixelMask(options, descriptor);
    const auto pieces =
        PlanPieces(descriptor, geometry, flag_geometry, request, options, decode_threads, static_cast<bool>(progress));

    std::vector<std::uint8_t> mask;
    // Reads `part` of `within` into `into` (the destination from where `within` lands), at `part`'s
    // place in it: the flag first, so that a flag that cannot be read leaves `into` as it was, then
    // the pixels, then the flag applied to them.
    const auto read_part = [&](const ReadRequest& part, const ReadRequest& within,
                               BufferView<float> into) -> Result<std::size_t> {
        auto selection = zarr::BuildSelectionWithin(descriptor, part, within);
        if (!selection) {
            return selection.error();
        }
        const auto part_elements = static_cast<std::size_t>(selection.value().elements());
        if (apply_mask) {
            // Dense in a buffer of its own, wherever the pixels land.
            auto dense = zarr::BuildSelection(descriptor, part, zarr::DestinationOrder::logical);
            if (!dense) {
                return dense.error();
            }
            mask.assign(part_elements, 0);
            // The read completes before it returns, so the destination remains valid for the next read.
            if (auto mask_read = source.ReadMask(dense.value(), {mask.data(), mask.size()}, options.control);
                !mask_read) {
                return mask_read.error();
            }
        }
        if (auto read = source.ReadPixels(selection.value(), into, options.control); !read) {
            return read.error();
        }
        if (apply_mask && selection.value().destination_extent() == part_elements) {
            ApplyPixelMask(into.data, mask.data(), part_elements);
        } else if (apply_mask) {
            ApplyPixelMaskWithin(selection.value(), into.data, mask.data());
        }
        return part_elements;
    };

    for (const auto& piece : pieces) {
        // The rest of the caller's buffer from where this piece lands, not the piece's own size: the
        // seam checks the piece against the buffer, and passing the piece's size would compare the
        // piece with itself. first_element is below the whole read's count, which CheckRead checked
        // against the buffer, so this never goes negative.
        const BufferView<float> piece_destination{destination.data + piece.first_element,
                                                  destination.size - static_cast<std::size_t>(piece.first_element)};
        std::size_t piece_elements = 0;
        if (piece.segments.empty()) {
            auto read = read_part(piece.request, piece.request, piece_destination);
            if (!read) {
                return read.error();
            }
            piece_elements = read.value();
        } else {
            for (const auto& segment : piece.segments) {
                auto read = read_part(segment, piece.request, piece_destination);
                if (!read) {
                    return read.error();
                }
                piece_elements += read.value();
            }
        }

        const auto finished = static_cast<std::size_t>(piece.first_element + piece_elements);
        if (progress && !progress(finished, static_cast<std::size_t>(elements))) {
            return Error{ErrorCode::cancelled, "The read was cancelled by its progress callback", descriptor.id};
        }
    }
    return static_cast<std::size_t>(elements);
}

namespace {

std::uint64_t Elements(const ReadRequest& request) {
    std::uint64_t elements = 1;
    for (const auto& range : request.axes) {
        elements *= range.count;
    }
    return elements;
}

// The most elements of a sample read at once when the caller sets no budget: 64 Ki, a quarter of a
// MiB of pixels. A sample is one element per chunk and usually a few thousand at most, but an image
// of small chunks can have millions per plane.
constexpr std::uint64_t kSampleElements = std::uint64_t{1} << 16;

// The most elements of `element_bytes` each that one piece of a sample may hold.
std::uint64_t SampleElements(const ReadOptions& options, std::size_t element_bytes) {
    return options.read_budget_bytes == 0 ? kSampleElements
                                          : std::max<std::uint64_t>(1, options.read_budget_bytes / element_bytes);
}

// Calls `read` with `request` cut into pieces of at most `most` elements, in turn, until one fails.
// The axes after some axis k are kept whole, k is cut into runs, and every axis before it is taken an
// element at a time. k is the first axis whose following axes fit together; when not even the last
// axis alone fits, the last axis is cut.
template <typename Read>
Result<void> ForEachPiece(const ReadRequest& request, std::uint64_t most, Read&& read) {
    const std::size_t rank = request.axes.size();
    std::size_t cut = 0;
    std::uint64_t inner = Elements(request);
    while (cut < rank) {
        inner /= std::max<std::uint64_t>(1, request.axes[cut].count);
        if (inner <= most) {
            break;
        }
        ++cut;
    }
    if (cut == rank) {
        return read(request);
    }
    const std::uint64_t run = std::max<std::uint64_t>(1, most / std::max<std::uint64_t>(1, inner));
    ReadRequest piece = request;
    std::vector<std::uint64_t> at(cut, 0);
    while (true) {
        for (std::size_t axis = 0; axis < cut; ++axis) {
            const auto& range = request.axes[axis];
            piece.axes[axis] = Range{range.start + (at[axis] * range.stride), 1, range.stride};
        }
        const auto& along = request.axes[cut];
        for (std::uint64_t first = 0; first < along.count; first += run) {
            piece.axes[cut] =
                Range{along.start + (first * along.stride), std::min(run, along.count - first), along.stride};
            if (auto done = read(piece); !done) {
                return done;
            }
        }
        // The next element of the axes before the cut, the last of them fastest.
        std::size_t axis = cut;
        while (axis > 0) {
            --axis;
            if (++at[axis] < request.axes[axis].count) {
                break;
            }
            at[axis] = 0;
            if (axis == 0) {
                return {};
            }
        }
        if (cut == 0) {
            return {};
        }
    }
}

}  // namespace

Result<std::uint64_t> PrefetchChunks(const PixelSource& source, const ImageDescriptor& descriptor,
                                     const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                                     const ReadRequest& request, const ReadOptions& options,
                                     std::size_t decode_threads) {
    // Checked as the request the caller made, so that a mistake in it is reported in its own terms
    // rather than in terms of the sample made from it.
    if (auto checked = zarr::BuildSelection(descriptor, request, zarr::DestinationOrder::logical); !checked) {
        return checked.error();
    }
    // Before anything is allocated, as a read checks: a cancelled prefetch holds no memory for it.
    if (auto allowed = zarr::CheckReadControl(options.control, descriptor.id); !allowed) {
        return allowed.error();
    }
    // The pixels and the flag are each sampled one element per chunk by their own layout, so the
    // pixels are read here without the flag, which is sampled below.
    //
    // Each in pieces, through one buffer the budget holds. A sample is one element per chunk, and an
    // image of small chunks has a great many: read whole, a 1 KiB budget would allocate 4 MiB for it.
    const auto sample = OneElementPerChunk(geometry, request);
    const auto chunks = Elements(sample);
    auto pixels_only = options;
    pixels_only.apply_pixel_mask = false;
    std::vector<float> discarded(static_cast<std::size_t>(std::min(chunks, SampleElements(options, sizeof(float)))));
    auto read = ForEachPiece(sample, discarded.size(), [&](const ReadRequest& piece) -> Result<void> {
        auto piece_read = ReadInPieces(source, descriptor, geometry, flag_geometry, piece,
                                       {discarded.data(), discarded.size()}, pixels_only, decode_threads, {});
        if (!piece_read) {
            return piece_read.error();
        }
        return {};
    });
    if (!read) {
        return read.error();
    }
    if (!AppliesPixelMask(options, descriptor)) {
        return chunks;
    }

    const auto flag_sample = OneElementPerChunk(flag_geometry, request);
    const auto flag_chunks = Elements(flag_sample);
    std::vector<std::uint8_t> discarded_flags(
        static_cast<std::size_t>(std::min(flag_chunks, SampleElements(options, sizeof(std::uint8_t)))));
    // Each element of the sample decodes a whole flag chunk, so a piece is limited to the chunks the
    // budget allows, not only by its buffer: limited by the buffer alone, one read of 1024 x 1024
    // flags chunked 32 x 32 would decode a MiB against a budget of 8 KiB. A flag chunk holds what any
    // chunk does for its decoded bytes, one byte per element. At least one chunk, as a read holds one
    // however small the budget. The pixels' pieces go through ReadInPieces, which cuts them by what
    // their chunks hold, as it cuts any read.
    const std::uint64_t flag_chunk_held =
        kHeldBytesPerDecodedByte * ChunkElements(flag_geometry.chunk_shape.empty() ? geometry : flag_geometry);
    const std::uint64_t flag_budget =
        options.read_budget_bytes != 0 ? options.read_budget_bytes : DefaultReadBytes(flag_chunk_held, decode_threads);
    const std::uint64_t flag_piece =
        std::min<std::uint64_t>(discarded_flags.size(), std::max<std::uint64_t>(1, flag_budget / flag_chunk_held));
    auto flags = ForEachPiece(flag_sample, flag_piece, [&](const ReadRequest& piece) -> Result<void> {
        auto selection = zarr::BuildSelection(descriptor, piece, zarr::DestinationOrder::logical);
        if (!selection) {
            return selection.error();
        }
        return source.ReadMask(selection.value(), {discarded_flags.data(), discarded_flags.size()}, options.control);
    });
    if (!flags) {
        return flags.error();
    }
    return chunks + flag_chunks;
}

}  // namespace carta::zarr::internal
