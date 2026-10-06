/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_
#define CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_

#include "store.h"

#include <tensorstore/context.h>
#include <tensorstore/tensorstore.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

namespace carta::zarr::internal {

// The resources shared by every read made through one public carta::zarr::Context. Array handles
// are kept by a per-Store clone, or by a CachePool's own context for the reads through it, so this
// object retains only shared TensorStore resources.
// Only translation units that talk to TensorStore include this header; store.h forward declares the
// type so that the schema layer never sees TensorStore.
class StoreContext : public std::enable_shared_from_this<StoreContext> {
public:
    explicit StoreContext(tensorstore::Context context);

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
     *
     * `through` is the pool a read has of its own, when it has one. A handle carries the pool it was
     * opened against and keeps what that pool decoded alive, so it is kept in the pool's table and
     * goes when the last holder of the pool lets go -- a pool is made for one task, and an image read
     * through it may stay open through any number of them. A store kept them for a while instead,
     * and five pools let go of left everything they had decoded alive until the image closed.
     *
     * In the pool's table under this store, by its id: a pool outlives the stores read through it,
     * and it once kept its handles by path alone, so a dataset closed, rewritten and opened again
     * through the same pool was read through the handle opened for the one before -- the chunks it
     * had decoded, and the data type it had then. An id is never reused, unlike an address. The
     * handles of a store that has closed are dropped the next time the pool opens one.
     */
    Result<tensorstore::TensorStore<>> OpenArray(const std::filesystem::path& array_path, std::string_view node,
                                                 const StoreContext* through = nullptr) const;

    /**
     * The same resources with a cache pool of `bytes` of its own. See carta::zarr::CachePool.
     *
     * A child context, so the thread pools are the shared ones -- a walk that ran on its own
     * threads would compete with the session rather than take its turn.
     *
     * It keeps its own array table because a handle carries the pool it was opened against, so a
     * read through this pool cannot reuse one opened with the shared pool -- one per store reading
     * through it, see OpenArray. A new one each call,
     * and kept by nothing here: whoever asked holds the only reference, so what the pool decoded is
     * freed when they let go. The zero-byte pool used to be built once and kept for as long as the
     * store, which cost nothing at that size and would keep gigabytes at the size a moment asks for.
     */
    Result<StoreContextPtr> WithCachePool(std::size_t bytes) const;

    tensorstore::Context context;

private:
    // The array at `array_path` of `store`, opened against this context and kept here.
    Result<tensorstore::TensorStore<>> OpenArrayOf(const StoreContext& store, const std::filesystem::path& array_path,
                                                   std::string_view node) const;

    struct Opened {
        // Who it was opened for, so that it can go once they have.
        std::weak_ptr<const StoreContext> store;
        tensorstore::TensorStore<> array;
    };

    // Never reused, unlike an address.
    std::uint64_t _id;
    mutable std::mutex _arrays_mutex;
    // By the id of the store it was opened for, and its path.
    mutable std::map<std::pair<std::uint64_t, std::string>, Opened> _arrays;
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

// Translate the public options into TensorStore context resources. A concurrency limit of zero
// keeps TensorStore's default; a cache size is written whenever the caller gave one, zero included,
// because zero is the pool that holds nothing rather than the absence of an answer.
Result<StoreContextPtr> MakeStoreContext(const ContextOptions& options);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_ZARR_STORE_CONTEXT_H_
