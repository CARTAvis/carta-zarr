/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_READ_H_
#define CARTA_ZARR_READ_H_

// Asking for pixels: which of them, into what, and under what limits.

#include "carta-zarr/descriptor.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace carta::zarr {

struct Range {
    std::uint64_t start = 0;
    std::uint64_t count = 0;
    std::uint64_t stride = 1;
};

struct ReadRequest {
    // One range per ImageDescriptor::axes entry, in the same order.
    std::vector<Range> axes;
};

// What a read should do with the decoded-chunk cache.
//
// A scan over a cube touches every chunk once and reuses none of them, so caching what it decodes
// evicts an interactive working set to no purpose -- and rebuilding that working set costs
// decompression, which is the resource the scan is already saturating. `bypass` runs the read
// against a cache pool of zero bytes, leaving the shared one alone.
//
// Arrays are opened per pool, so the first bypassed read of an array pays to open it again. That is
// once per array, against a scan that reads all of it.
enum class CachePolicy {
    inherit,
    bypass,
};

// Called as a read advances, with the number of destination elements that are final and the number
// the request will produce in total. Returning false cancels the read, which then reports cancelled.
//
// Supplying one splits the read into chunk-aligned pieces along the slowest-varying selected axis,
// so that there is somewhere to report from and somewhere to stop. The destination is dense in
// logical order with axis 0 fastest, which is what makes the finished part a prefix rather than a
// scatter -- a caller can render or forward it as it arrives.
//
// It is not the only thing that splits a read; ReadOptions::temporary_memory_limit_bytes does too,
// and a read with neither is issued in one piece.
//
// An argument of Image::Read rather than a field of ReadOptions, because it is the only operation
// that has anywhere to report from in these terms -- a reduction reports through its sink, and a
// cube histogram through a CubeHistogramProgressCallback of its own. As a field it was a field four
// of the five entry points silently ignored; as an argument it is simply not part of what they take.
//
// A read that nothing interrupts is not made slower by supplying one: the pieces are sized to hold
// enough chunks to decode in parallel, and at that size a split read measures the same as an
// unsplit one.
using ProgressCallback = std::function<bool(std::size_t elements_written, std::size_t elements_total)>;

// What a read is allowed to do while it runs, whatever it is reading for.
//
// These three are the whole of what every path through this library honours: an ordinary read and
// all three reductions reach the same storage operations underneath and check the same things at
// the same boundaries. Said in its own type so that an operation which honours only these can take
// only these -- everything below the pixel source seam does, because that is all any of it ever
// looked at.
struct ReadControl {
    // Cooperative cancellation checked before and after each storage operation. The callback
    // must be safe to invoke from the calling thread.
    std::function<bool()> cancellation_requested;
    // A steady-clock deadline checked at the same storage-operation boundaries. An in-flight
    // TensorStore operation is not interrupted, but a request never starts another operation once
    // this deadline has passed.
    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max();
    // Whether this read may put what it decodes in the shared cache. See CachePolicy.
    CachePolicy cache_policy = CachePolicy::inherit;
};

// Everything a read of pixels takes, on top of what any read takes.
//
// The two fields here are the ones that mean something only when pixels are being read into a
// buffer this library sized: whether a flag is folded in on the way, and how much the library may
// hold at once while doing it. A pixel mask read has neither -- it is the flag, and the destination
// is the caller's -- which is why it takes a ReadControl and this cannot be handed to it.
//
// That is the point of the split. Every field of this type is honoured by every operation that
// takes it, so there is no table of which ones apply where, and setting one where it would have
// been ignored does not compile.
struct ReadOptions {
    ReadControl control;
    // Write NaN wherever the pixel mask is false, so that one call answers what would otherwise be
    // a pixel read plus a mask read. On by default: masking during the read costs one pass over
    // data already in hand, while a caller doing it afterwards pays for a second traversal.
    bool apply_pixel_mask = true;
    // Maximum temporary memory one piece of the read may use. Zero means the library's own budget.
    //
    // This bounds the pixel mask buffer, and it is also what a split read sizes its pieces by --
    // both are "how much this read may hold at once", and splitting to fit is a better answer than
    // refusing. So setting it splits the read, whether or not anyone asked to watch: a ceiling is a
    // statement about memory, not about wanting progress. A read that cannot be split far enough --
    // no axis selects more than one element, or a piece of one chunk is still too large -- reports
    // buffer_too_small rather than allocating past the limit.
    //
    // For ReduceSpectral it is a target rather than a limit. That walk splits along x, along the
    // chunk rows and along the spectrum, and each of the three bottoms out at one chunk, which is
    // the smallest thing that can be decoded: asking for part of a chunk decodes all of it anyway,
    // and asking twice decodes it twice. So an image whose chunk is larger than this exceeds it by
    // the ratio, and refusing to reduce would be the worse answer. ChunkGeometry::chunk_shape says
    // in advance when that will happen.
    std::size_t temporary_memory_limit_bytes = 0;
};

// Somewhere for a read to write, counted in elements rather than in bytes.
//
// Typed because the element type is not the caller's to choose: pixels arrive as float, and a pixel
// mask as one byte per pixel. An untyped view with a byte count could express neither fact, so the
// same field meant "bytes" at one entry point and "elements" at the other, and a buffer of the
// wrong kind was a run-time error at best.
//
// The library never holds one past the call it was passed to.
template <typename T>
struct BufferView {
    T* data = nullptr;
    std::size_t size = 0;  // elements, not bytes
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_READ_H_
