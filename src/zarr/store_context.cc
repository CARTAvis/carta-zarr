/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "store_context.h"

#include <tensorstore/data_type.h>
#include <tensorstore/index.h>
#include <tensorstore/open.h>
#include <tensorstore/open_mode.h>
#include <tensorstore/spec.h>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

namespace carta::zarr::internal {
namespace {

// What a carta-zarr array is, stated once. Every array this library reads -- pixels, flags and
// coordinate values alike -- is opened through here.
Result<tensorstore::TensorStore<>> OpenZarr3File(const std::string& path, const tensorstore::Context& context,
                                                 std::string_view node) {
    auto spec = tensorstore::Spec::FromJson({
        {"driver", "zarr3"},
        {"kvstore", {{"driver", "file"}, {"path", path}}},
        // A chunk cached after the array was opened is used without asking storage whether it has
        // changed. TensorStore's default asks on every read, which on a parallel file system is a
        // metadata round trip per chunk even when the cache holds it. ADR 0015.
        {"recheck_cached_data", "open"},
    });
    if (!spec.ok()) {
        return Error{ErrorCode::io_error, "Failed to create TensorStore spec: " + spec.status().ToString(),
                     std::string(node)};
    }
    auto opened = tensorstore::Open(spec.value(), context, tensorstore::OpenMode::open,
                                    tensorstore::ReadWriteMode::read)
                      .result();
    if (!opened.ok()) {
        return Error{ErrorCode::io_error, "Failed to open TensorStore: " + opened.status().ToString(),
                     std::string(node)};
    }
    return std::move(opened).value();
}

// The one thing about a data type that cannot go in the table beside the others: a
// tensorstore::dtype_v is a template, so the mapping has to be written as code. It goes through the
// shared table for the half that can be shared -- what the name means -- and spells out only the
// half that cannot.
//
// A type this does not answer for is a type no read accepts, which is why complex is absent: nothing
// reaches here holding one, because an image and its coordinates are required to be real long
// before a value of them is asked for.
bool MatchesDataType(std::string_view expected, tensorstore::DataType actual) {
    switch (zarr::ParseDataType(expected)) {
        case DataType::boolean: return actual == tensorstore::dtype_v<bool>;
        case DataType::int8: return actual == tensorstore::dtype_v<std::int8_t>;
        case DataType::uint8: return actual == tensorstore::dtype_v<std::uint8_t>;
        case DataType::int16: return actual == tensorstore::dtype_v<std::int16_t>;
        case DataType::uint16: return actual == tensorstore::dtype_v<std::uint16_t>;
        case DataType::int32: return actual == tensorstore::dtype_v<std::int32_t>;
        case DataType::uint32: return actual == tensorstore::dtype_v<std::uint32_t>;
        case DataType::int64: return actual == tensorstore::dtype_v<std::int64_t>;
        case DataType::uint64: return actual == tensorstore::dtype_v<std::uint64_t>;
        case DataType::float16: return actual == tensorstore::dtype_v<tensorstore::dtypes::float16_t>;
        case DataType::float32: return actual == tensorstore::dtype_v<float>;
        case DataType::float64: return actual == tensorstore::dtype_v<double>;
        default: return false;
    }
}

// A pool of `bytes`, said once for the session's pool and for a read's own.
nlohmann::json CachePoolSpec(std::size_t bytes) {
    return {{"total_bytes_limit", bytes}};
}

}  // namespace

StoreContextPtr StoreContext::CloneForStore() const {
    return std::make_shared<const StoreContext>(context);
}

Result<tensorstore::TensorStore<>> StoreContext::OpenArray(const std::filesystem::path& array_path,
                                                           std::string_view node) const {
    const std::string key = array_path.string();
    {
        const std::scoped_lock lock(_arrays_mutex);
        if (const auto found = _arrays.find(key); found != _arrays.end()) {
            return found->second;
        }
    }

    // Opened outside the lock so that a slow open of one array does not stall reads of another.
    auto opened = OpenZarr3File(key, context, node);
    if (!opened) {
        return opened.error();
    }

    const std::scoped_lock lock(_arrays_mutex);
    // Another thread may have opened the same array first; either handle is equivalent, so keep
    // whichever landed in the table.
    return _arrays.emplace(key, std::move(opened.value())).first->second;
}

Result<tensorstore::TensorStore<>> OpenZarrArray(const std::filesystem::path& array_directory,
                                                 const StoreContextPtr& context, std::string_view node) {
    if (context) {
        return context->OpenArray(array_directory, node);
    }
    return OpenZarr3File(array_directory.string(), tensorstore::Context::Default(), node);
}

Result<void> VerifyArrayMatchesMetadata(const tensorstore::TensorStore<>& array,
                                        const zarr::ArrayMetadata& expected, std::string_view node) {
    const auto rank = expected.shape.size();
    if (static_cast<std::size_t>(array.rank()) != rank) {
        return Error{ErrorCode::invalid_metadata, "Array rank differs between metadata sources",
                     std::string(node)};
    }
    const auto actual_shape = array.domain().shape();
    if (actual_shape.size() != rank) {
        return Error{ErrorCode::invalid_metadata, "Array rank differs between metadata sources",
                     std::string(node)};
    }
    for (std::size_t axis = 0; axis < rank; ++axis) {
        if (actual_shape[axis] != static_cast<tensorstore::Index>(expected.shape.at(axis))) {
            return Error{ErrorCode::invalid_metadata,
                         "Array shape differs between canonical metadata and the array store",
                         std::string(node)};
        }
    }
    // Two axes of one length exchanged leave every extent agreeing, and a read addressed by the
    // store's names then puts one axis's values where the other's should be. A 1-D array whose own
    // document names nothing is let through: it has no order to get wrong.
    const auto actual_dimension_names = array.domain().labels();
    const bool unnamed_vector =
        rank == 1 && actual_dimension_names.size() == 1 && actual_dimension_names[0].empty();
    if (!unnamed_vector) {
        if (actual_dimension_names.size() != rank || expected.dimension_names.size() != rank) {
            return Error{ErrorCode::invalid_metadata, "Array dimension names differ between metadata sources",
                         std::string(node)};
        }
        for (std::size_t axis = 0; axis < rank; ++axis) {
            if (actual_dimension_names[axis] != expected.dimension_names.at(axis)) {
                return Error{ErrorCode::invalid_metadata,
                             "Array dimension names differ between canonical metadata and the array store",
                             std::string(node)};
            }
        }
    }
    if (!MatchesDataType(expected.data_type, array.dtype())) {
        return Error{ErrorCode::invalid_metadata,
                     "Array data type differs between canonical metadata and the array store",
                     std::string(node)};
    }
    return {};
}

Result<StoreContextPtr> StoreContext::WithCachePool(std::size_t bytes) const {
    nlohmann::json spec = nlohmann::json::object();
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
    spec["cache_pool"] = CachePoolSpec(bytes);
    auto child = tensorstore::Context::FromJson(std::move(spec), context);
    if (!child.ok()) {
        return Error{ErrorCode::invalid_argument, "Unable to make a cache pool: " + child.status().ToString(), {}};
    }
    return std::make_shared<const StoreContext>(std::move(child.value()));
}

Result<StoreContextPtr> MakeStoreContext(const ContextOptions& options) {
    nlohmann::json spec = nlohmann::json::object();
    if (options.cache_bytes) {
        // Zero is a size like any other here: it is the pool that holds nothing, which is what a
        // caller declining the cache asks for. No value at all is the only thing that leaves
        // TensorStore's default in place.
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        spec["cache_pool"] = CachePoolSpec(*options.cache_bytes);
    }
    if (options.io_threads > 0) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        spec["file_io_concurrency"] = {{"limit", options.io_threads}};
    }
    if (options.decode_threads > 0) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-avoid-unchecked-container-access)
        spec["data_copy_concurrency"] = {{"limit", options.decode_threads}};
    }

    if (spec.empty()) {
        return std::make_shared<const StoreContext>(tensorstore::Context::Default());
    }

    auto context = tensorstore::Context::FromJson(std::move(spec));
    if (!context.ok()) {
        return Error{ErrorCode::invalid_argument,
                     "Unable to apply the requested resource limits: " + context.status().ToString(), {}};
    }
    return std::make_shared<const StoreContext>(std::move(context.value()));
}

}  // namespace carta::zarr::internal
