/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a pool of its own is made of, asked of the TensorStore context underneath. A consumer sees a
// CachePool only by reading through one, and every read through one returns the same pixels whatever
// its size, so the three things that make it worth having -- that it is the size asked for, that it
// runs on the session's threads rather than its own, and that nothing keeps it once its holder lets
// go -- are visible only here.

#include "zarr/store_context.h"

#include <tensorstore/internal/cache/cache_pool_resource.h>
#include <tensorstore/internal/data_copy_concurrency_resource.h>
#include <tensorstore/internal/file_io_concurrency_resource.h>

#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include "support/check.h"

namespace {

using carta::zarr::testing::Require;
using carta::zarr::internal::StoreContextPtr;

std::size_t PoolLimit(const StoreContextPtr& store_context) {
    auto pool = store_context->context.GetResource<tensorstore::internal::CachePoolResource>();
    Require(pool.ok(), "the context has no cache pool resource");
    const auto& weak = **pool;
    Require(static_cast<bool>(weak), "the cache pool resource is empty");
    return weak->limits().total_bytes_limit;
}

template <typename Provider>
const void* ResourceOf(const StoreContextPtr& store_context) {
    auto resource = store_context->context.GetResource<Provider>();
    Require(resource.ok(), std::string("the context has no ") + Provider::id + " resource");
    return resource->get();
}

StoreContextPtr Session() {
    carta::zarr::ContextOptions options;
    options.cache_bytes = std::size_t{8} << 20;
    options.io_threads = 3;
    options.decode_threads = 5;
    auto session = carta::zarr::internal::MakeStoreContext(options);
    Require(static_cast<bool>(session), "the session context could not be made");
    return session.value();
}

// The size asked for, zero included: zero is the pool that holds nothing, which is what a scan that
// should leave the session's working set alone asks for.
void TestAPoolIsTheSizeAskedFor() {
    const auto session = Session();
    for (const std::size_t bytes : {std::size_t{0}, std::size_t{1}, std::size_t{3} << 30}) {
        auto own = session->WithCachePool(bytes);
        Require(static_cast<bool>(own), "a pool of " + std::to_string(bytes) + " bytes could not be made");
        Require(PoolLimit(own.value()) == bytes,
                "a pool of " + std::to_string(bytes) + " bytes holds " + std::to_string(PoolLimit(own.value())));
    }
    Require(PoolLimit(session) == (std::size_t{8} << 20), "making a pool of its own changed the session's");
}

// On the session's threads: a walk that ran on threads of its own would compete with the session's
// reads rather than take its turn among them.
void TestAPoolRunsOnTheSessionsThreads() {
    const auto session = Session();
    auto own = session->WithCachePool(std::size_t{1} << 20);
    Require(static_cast<bool>(own), "a pool of its own could not be made");
    Require(ResourceOf<tensorstore::internal::DataCopyConcurrencyResource>(own.value()) ==
                ResourceOf<tensorstore::internal::DataCopyConcurrencyResource>(session),
            "a pool of its own decodes on threads other than the session's");
    Require(ResourceOf<tensorstore::internal::FileIoConcurrencyResource>(own.value()) ==
                ResourceOf<tensorstore::internal::FileIoConcurrencyResource>(session),
            "a pool of its own reads files on threads other than the session's");
    Require(ResourceOf<tensorstore::internal::CachePoolResource>(own.value()) !=
                ResourceOf<tensorstore::internal::CachePoolResource>(session),
            "a pool of its own is the session's pool");
}

// Nothing keeps a pool once its holder lets go, and two asked for are two. The bypass this replaces
// was built once and kept for as long as the store, which cost nothing at zero bytes and would keep
// gigabytes of decoded chunks at the size a moment asks for.
void TestNothingKeepsAPoolItsHolderLetGo() {
    const auto session = Session();
    std::weak_ptr<const carta::zarr::internal::StoreContext> first;
    std::weak_ptr<const carta::zarr::internal::StoreContext> last;
    {
        auto own = session->WithCachePool(std::size_t{1} << 20);
        Require(static_cast<bool>(own), "a pool of its own could not be made");
        first = own.value();
        auto other = session->WithCachePool(std::size_t{1} << 20);
        Require(static_cast<bool>(other) && other.value() != own.value(), "two pools asked for are one");
        last = other.value();
    }
    Require(first.expired() && last.expired(), "a pool outlived the last holder of it");
}

}  // namespace

int main() {
    try {
        TestAPoolIsTheSizeAskedFor();
        TestAPoolRunsOnTheSessionsThreads();
        TestNothingKeepsAPoolItsHolderLetGo();
    } catch (const std::exception& error) {
        std::cerr << "store context test failed: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
