/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_
#define CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_

#include "store.h"

#include <tensorstore/context.h>
#include <tensorstore/tensorstore.h>

#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>

namespace carta::zarr::internal {

// The resources shared by every read made through one public carta::zarr::Context. Array handles
// are cloned into a per-Store context, so this object retains only shared TensorStore resources.
// Only translation units that talk to TensorStore include this header; store.h forward declares the
// type so that the schema layer never sees TensorStore.
class StoreContext {
public:
    explicit StoreContext(tensorstore::Context context) : context(std::move(context)) {}

    // Make a context for one Dataset/Store. The TensorStore resource handles remain shared, while
    // the array-handle table has the same lifetime as that store rather than the public Context.
    StoreContextPtr CloneForStore() const;

    /**
     * Open an array once and hand back the same handle afterwards.
     *
     * Opening is not cheap: it reads and parses the array's metadata and builds a driver. A pixel
     * read is issued once per casacore cursor step, so paying that per call turns a fixed cost into
     * a per-call one and dominates everything else -- measured at several milliseconds a call,
     * against roughly a hundred milliseconds for a whole 4096-square plane read in one go.
     *
     * The handle is cheap to copy and safe to use from several threads, so callers get a copy and
     * the table is only locked around the lookup.
     */
    Result<tensorstore::TensorStore<>> OpenArray(const std::filesystem::path& array_path,
                                                 std::string_view node) const;

    /**
     * The same resources with a cache pool of `bytes` of its own. See carta::zarr::CachePool.
     *
     * A child context, so the thread pools are the shared ones -- a walk that ran on its own
     * threads would compete with the session rather than take its turn.
     *
     * It keeps its own array table because a handle carries the pool it was opened against, so a
     * read through this pool cannot reuse one opened with the shared pool. A new one each call,
     * and kept by nothing here: whoever asked holds the only reference, so what the pool decoded is
     * freed when they let go. The zero-byte pool used to be built once and kept for as long as the
     * store, which cost nothing at that size and would keep gigabytes at the size a moment asks for.
     */
    Result<StoreContextPtr> WithCachePool(std::size_t bytes) const;

    tensorstore::Context context;

private:
    mutable std::mutex _arrays_mutex;
    mutable std::map<std::string, tensorstore::TensorStore<>> _arrays;
};

// How CachePool's implementation reaches the context it reads through, which is internal and has no
// business in the public header.
struct CachePoolAccess {
    static const StoreContextPtr& StoreContextOf(const CachePool& pool);
};

}  // namespace carta::zarr::internal

namespace carta::zarr {

// The size a pool was asked for, beside the context that holds it; TensorStore keeps the size too,
// but only behind a resource lookup that can fail.
class CachePool::Impl {
public:
    Impl(std::size_t bytes, internal::StoreContextPtr store_context)
        : bytes(bytes), store_context(std::move(store_context)) {}

    std::size_t bytes;
    internal::StoreContextPtr store_context;
};

}  // namespace carta::zarr

namespace carta::zarr::internal {

inline const StoreContextPtr& CachePoolAccess::StoreContextOf(const CachePool& pool) {
    return pool._impl->store_context;
}

/**
 * Open the Zarr array held in a located directory.
 *
 * The one place this library says what its arrays are made of -- zarr3 over a file kvstore -- so a
 * transport that reaches bytes some other way has one function to change rather than a search to
 * run. ADR 0004 records why the transport still hands up a filesystem path for one to be built
 * from, and what that costs.
 *
 * With a context the handle comes from that context's table and is shared. Without one, TensorStore
 * default resources are used and nothing is kept; that is what a probe gets, and a probe opens no
 * arrays.
 */
Result<tensorstore::TensorStore<>> OpenZarrArray(const std::filesystem::path& array_directory,
                                                 const StoreContextPtr& context, std::string_view node);

/**
 * Whether an opened array is the one the store's canonical metadata describes: its rank, its
 * extent, the names and order of its dimensions, and its data type.
 *
 * The store decides what an image is from the document it parsed -- with consolidated metadata, the
 * root's copy -- while TensorStore opens an array from the array's own. The two are meant to be the
 * same document and are not always: a copy left stale by a rewrite, or a store written by hand. An
 * array that disagrees is refused as invalid_metadata before any value of it is read, whether it
 * holds pixels or a coordinate, because a coordinate of another length describes another image.
 */
Result<void> VerifyArrayMatchesMetadata(const tensorstore::TensorStore<>& array,
                                        const zarr::ArrayMetadata& expected, std::string_view node);

// Translate the public options into TensorStore context resources. A concurrency limit of zero
// keeps TensorStore's default; a cache size is written whenever the caller gave one, zero included,
// because zero is the pool that holds nothing rather than the absence of an answer.
Result<StoreContextPtr> MakeStoreContext(const ContextOptions& options);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_
