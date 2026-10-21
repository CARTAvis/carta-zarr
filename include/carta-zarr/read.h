/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_READ_H_
#define CARTA_ZARR_READ_H_

// Asking for pixels: which of them, into what, and under what limits.

#include "carta-zarr/descriptor.h"
#include "carta-zarr/export.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
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

namespace internal {
struct CachePoolAccess;
}  // namespace internal

// A decoded-chunk cache of a read's own, apart from the one its Context shares among every read.
//
// The shared cache holds a session's interactive working set: the planes and profiles being looked
// at. Two kinds of read do better with a cache of their own:
//
// - A scan that touches every chunk once reuses none of them. Caching what it decodes would only
//   evict the working set, so it reads through a pool of zero bytes.
// - A walk that comes back to its chunks, such as a moment whose neighboring slabs share chunks,
//   needs room for all of them, but only while it runs. It reads through a pool of that size and
//   lets go of it when it ends.
//
// Made by Context::CreateCachePool. It is a shared handle like Context: copies share the pool, and the
// last copy to go frees what it holds. The pool runs on its Context's threads, so its reads take
// their turn among the session's reads. Arrays are opened once per pool, so the first read of an
// array through a new pool opens it again, which takes a few milliseconds.
class CARTA_ZARR_EXPORT CachePool final {
public:
    CachePool(const CachePool&) = default;
    CachePool& operator=(const CachePool&) = default;
    // A copy, so that the handle moved from still refers to what it did. See ADR 0013.
    CachePool(CachePool&& other) noexcept : CachePool(other) {}
    CachePool& operator=(CachePool&& other) noexcept { return *this = other; }
    // Inline, unlike the other handles' destructors: ReadControl holds a CachePool, and some parts of
    // this library are built and tested on their own with a ReadControl. A shared_ptr's deleter is
    // fixed where it was made, so this does not need Impl to be complete.
    ~CachePool() = default;

    // How much decoded chunk data the pool may hold, in bytes, as it was asked for.
    std::size_t bytes() const noexcept;

private:
    class Impl;
    explicit CachePool(std::shared_ptr<const Impl> impl);

    std::shared_ptr<const Impl> _impl;

    friend class Context;
    friend struct internal::CachePoolAccess;
};

/// Called as a read advances, with the number of destination elements that are final and the number
/// the request will produce in total. Returning false cancels the read, which then reports cancelled.
///
/// A read given one is cut into chunk-aligned pieces along the slowest-varying selected axis: as
/// many as ReadOptions::read_budget_bytes needs, and each decoding at most sixteen chunks in all,
/// unless one chunk along that axis already spans more. This is called as each piece finishes, about
/// as often with many decode threads as with few.
///
/// The destination is dense in logical order with axis 0 fastest, so the finished part is always a
/// prefix, which a caller can render or forward as it arrives. A piece too large to read at once is
/// read in parts and reported when all of them are in. The prefix holds after a read fails as well
/// as while it runs; see Image::Read.
///
/// It is an argument of Image::Read, not a field of ReadOptions, because no other operation reports
/// progress this way: a reduction reports through its sink, and a cube histogram through its own
/// CubeHistogramProgressCallback. As a field, the other operations would silently ignore it.
///
/// Supplying one can cut a read finer. Every read is cut to fit its budget, but a read with a
/// callback is also cut to sixteen chunks, so that progress arrives often. A single plane is read
/// whole only without one.
using ProgressCallback = std::function<bool(std::size_t elements_written, std::size_t elements_total)>;

// What a read is allowed to do while it runs, whatever it is reading for.
//
// Every operation in this library honors these three: an ordinary read and all three reductions
// reach the same storage operations and check them at the same points. They are a type of their own
// so that the code below the pixel source seam, which uses nothing else, can take only these.
struct ReadControl {
    // Cooperative cancellation checked before and after each storage operation. The callback
    // must be safe to invoke from the calling thread.
    std::function<bool()> cancellation_requested;
    // A steady-clock deadline checked at the same storage-operation boundaries. An in-flight
    // TensorStore operation is not interrupted, but a request never starts another operation once
    // this deadline has passed.
    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max();
    // The cache this read keeps what it decodes in, and finds what an earlier read kept: its Context's
    // shared one when there is none. See CachePool.
    //
    // To bypass the cache, pass a pool of zero bytes.
    std::optional<CachePool> cache_pool;
};

// Everything a read of pixels takes, on top of what any read takes.
//
// The two fields here only mean something when pixels are read: whether a flag is applied on the
// way, and how much the library may hold at once. A pixel mask read has neither, since it reads the
// flag itself into the caller's buffer, so it takes only a ReadControl.
//
// Every operation that takes this type honors every field of it, so no field is silently ignored.
struct ReadOptions {
    ReadControl control;
    // Write NaN wherever the image's flag is true, meaning the pixel is bad, so that one call does
    // the work of a pixel read and a flag read. On by default: masking during the read costs one pass
    // over data already in memory, while masking afterwards costs a second pass.
    bool apply_pixel_mask = true;
    /// How much memory one read may hold at once beyond the caller's own destination, in bytes. Zero
    /// means the library's own budget: two chunks for every decode thread of the image's context,
    /// between 256 MiB and 2 GiB.
    ///
    /// A read holds the chunks it is decoding, each about three times its decoded size (the
    /// compressed bytes, the codec's buffer and the decoded copy), and the buffers the library
    /// allocates for them: the pixel mask, one byte per element, and the pixels wherever the library
    /// rather than the caller holds them. Image::Read and Image::Prefetch cut their request into
    /// pieces that fit. ReduceSpectral, ComputeHistogram and ComputeCubeHistogram size each read of
    /// their walk by it.
    ///
    /// One exception: a single plane read with no progress callback is read whole when the budget
    /// allows one chunk per decode thread. Cut to the budget it runs a quarter slower, and the plane
    /// bounds what it holds. See ProgressCallback and ADR 0021.
    ///
    /// A chunk is the smallest thing that can be decoded: asking for part of one decodes all of it.
    /// So a read holds at least one chunk, and under a budget smaller than that it reads one chunk at
    /// a time and holds what that chunk does rather than refusing. Image::DecodedChunkBytes says in
    /// advance what a chunk decodes to.
    ///
    /// The budget is per read. Reads running at once each hold their own.
    std::size_t read_budget_bytes = 0;
};

// A run of elements the caller owns and lends for one call, counted in elements rather than in
// bytes: somewhere for a read to write, or, as BufferView<const T>, something for a reduction to read
// -- its regions and each region's raster -- whose length the library checks rather than assumes.
//
// Typed because the element type is not the caller's to choose: pixels arrive as float, and a pixel
// mask as one byte per pixel. A buffer of the wrong type does not compile.
//
// The library never holds one past the call it was passed to.
template <typename T>
struct BufferView {
    T* data = nullptr;
    std::size_t size = 0;  // elements, not bytes
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_READ_H_
