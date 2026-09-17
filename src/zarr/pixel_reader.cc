/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "pixel_reader.h"

#include "store_context.h"

#include <tensorstore/array.h>
#include <tensorstore/cast.h>
#include <tensorstore/context.h>
#include <tensorstore/data_type.h>
#include <tensorstore/index.h>
#include <tensorstore/index_space/dim_expression.h>
#include <tensorstore/internal/unowned_to_shared.h>
#include <tensorstore/tensorstore.h>
#include <tensorstore/util/result.h>

#include <algorithm>
#include <chrono>
#include <exception>
#include <limits>
#include <string>
#include <utility>

namespace carta::zarr::internal::zarr {
namespace {

// The one thing about a data type that cannot go in the table beside the others: a
// tensorstore::dtype_v is a template, so the mapping has to be written as code. It goes through the
// shared table for the half that can be shared -- what the name means -- and spells out only the
// half that cannot.
//
// A type this does not answer for is a type a pixel read refuses, which is why complex is absent:
// nothing reaches here holding one, because an image is required to be real long before a pixel of
// it is asked for.
bool MatchesDataType(std::string_view expected, tensorstore::DataType actual) {
    switch (ParseDataType(expected)) {
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

bool SelectionIsWellFormed(const PixelSelection& selection) {
    const auto rank = selection.start.size();
    if (rank == 0 || selection.count.size() != rank || selection.stride.size() != rank || selection.shape.size() != rank ||
        selection.dimension_names.size() != rank ||
        selection.logical_to_stored.size() != rank) {
        return false;
    }
    std::vector<bool> seen(rank, false);
    for (const auto stored : selection.logical_to_stored) {
        if (stored >= rank || seen.at(stored)) {
            return false;
        }
        seen.at(stored) = true;
    }
    return std::all_of(selection.stride.begin(), selection.stride.end(),
                       [](std::uint64_t value) { return value > 0; });
}

// What the array on disk has to agree with the store's canonical metadata about before a single
// pixel of it is read: its rank, its extent, what its dimensions are called, and what it holds.
//
// Its own function because it is the one stretch of ReadInto that decides nothing about the read.
// What is left reads as the seven steps it is: check the request, open the array, verify it, slice,
// transpose, convert, read.
Result<void> VerifyStoreMatchesSelection(const tensorstore::TensorStore<>& store, const PixelSelection& selection,
                                         std::string_view expected_data_type, std::string_view node) {
    const auto rank = selection.start.size();
    if (static_cast<std::size_t>(store.rank()) != rank) {
        return Error{ErrorCode::invalid_argument, "Selection rank does not match the array rank",
                     std::string(node)};
    }
    const auto actual_shape = store.domain().shape();
    if (actual_shape.size() != rank) {
        return Error{ErrorCode::invalid_metadata, "Array rank differs between metadata sources",
                     std::string(node)};
    }
    for (std::size_t axis = 0; axis < rank; ++axis) {
        if (actual_shape[axis] != static_cast<tensorstore::Index>(selection.shape.at(axis))) {
            return Error{ErrorCode::invalid_metadata,
                         "Array shape differs between canonical metadata and the array store",
                         std::string(node)};
        }
    }
    const auto actual_dimension_names = store.domain().labels();
    if (actual_dimension_names.size() != rank) {
        return Error{ErrorCode::invalid_metadata, "Array dimension names differ between metadata sources",
                     std::string(node)};
    }
    for (std::size_t axis = 0; axis < rank; ++axis) {
        if (actual_dimension_names[axis] != selection.dimension_names.at(axis)) {
            return Error{ErrorCode::invalid_metadata,
                         "Array dimension names differ between canonical metadata and the array store",
                         std::string(node)};
        }
    }
    if (!MatchesDataType(expected_data_type, store.dtype())) {
        return Error{ErrorCode::invalid_metadata,
                     "Array data type differs between canonical metadata and the array store",
                     std::string(node)};
    }
    return {};
}

template <typename Element>
Result<void> ReadInto(const std::filesystem::path& array_path, const StoreContextPtr& context,
                      std::string_view node, std::string_view expected_data_type, const PixelSelection& selection,
                      tensorstore::DataType target_dtype,
                      Element* destination, std::size_t destination_elements, const ReadOptions& options) {
    if (destination == nullptr) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is null", std::string(node)};
    }
    if (!SelectionIsWellFormed(selection)) {
        return Error{ErrorCode::invalid_argument, "Malformed pixel selection", std::string(node)};
    }
    const auto elements = SelectionElementCount(selection);
    if (elements == 0) {
        return Error{ErrorCode::invalid_argument, "Pixel selection is empty", std::string(node)};
    }
    if (elements > destination_elements) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is too small", std::string(node)};
    }

    if (!context) {
        return Error{ErrorCode::invalid_argument, "Pixel reads require a context", std::string(node)};
    }

    auto control = CheckReadControl(options, node);
    if (!control) {
        return control.error();
    }

    try {
        // Reused across calls. A slice read happens once per casacore cursor step, so opening here
        // would make a fixed cost a per-call one.
        // A read that says it should not pollute the shared cache runs against a child context
        // whose pool holds nothing. Chosen here rather than by the caller because this is where the
        // array handle is taken, and a handle carries the pool it was opened against.
        const StoreContextPtr pool =
            options.cache_policy == CachePolicy::bypass ? context->WithoutCache() : context;
        auto opened = pool->OpenArray(array_path, node);
        if (!opened) {
            return opened.error();
        }
        control = CheckReadControl(options, node);
        if (!control) {
            return control.error();
        }
        auto const store = std::move(opened).value();
        if (auto agreed = VerifyStoreMatchesSelection(store, selection, expected_data_type, node); !agreed) {
            return agreed.error();
        }
        const auto rank = selection.start.size();
        std::vector<tensorstore::Index> start(rank);
        std::vector<tensorstore::Index> count(rank);
        std::vector<tensorstore::Index> stride(rank);
        std::vector<tensorstore::DimensionIndex> order(rank);
        for (std::size_t i = 0; i < rank; ++i) {
            start.at(i) = static_cast<tensorstore::Index>(selection.start.at(i));
            count.at(i) = static_cast<tensorstore::Index>(selection.count.at(i));
            stride.at(i) = static_cast<tensorstore::Index>(selection.stride.at(i));
            order.at(i) = static_cast<tensorstore::DimensionIndex>(selection.logical_to_stored.at(i));
        }

        // Slice in stored order, then move the stored dimensions into logical order. Both are index
        // transforms, so TensorStore composes them into the one copy the read already performs.
        auto sliced = store | tensorstore::AllDims().TranslateSizedInterval(start, count, stride);
        if (!sliced.ok()) {
            return Error{ErrorCode::invalid_argument,
                         "Requested region is outside the array: " + sliced.status().ToString(),
                         std::string(node)};
        }
        auto transposed = std::move(sliced).value() | tensorstore::Dims(order).Transpose();
        if (!transposed.ok()) {
            return Error{ErrorCode::invalid_argument,
                         "Failed to reorder axes: " + transposed.status().ToString(), std::string(node)};
        }

        // Conversion rides the same copy, so a float64 or integer array is never materialized in
        // its stored type first.
        auto converted = tensorstore::Cast(std::move(transposed).value(), target_dtype);
        if (!converted.ok()) {
            return Error{ErrorCode::unsupported_data_type,
                         "Array cannot be converted to the requested output type: " +
                             converted.status().ToString(),
                         std::string(node)};
        }

        std::vector<tensorstore::Index> shape(rank);
        for (std::size_t i = 0; i < rank; ++i) {
            shape.at(i) = count.at(order.at(i));
        }
        // Axis 0 is the fastest-varying destination dimension, so a logical-order shape over a
        // densely packed buffer is Fortran-ordered. TensorStore requires a shared array here; the
        // caller owns this buffer and the read below is awaited before returning, so a non-owning
        // shared pointer is what the ownership actually is rather than a way around the check.
        auto target = tensorstore::Array(tensorstore::internal::UnownedToShared(destination), shape,
                                         tensorstore::fortran_order);

        auto const read_result = tensorstore::Read(std::move(converted).value(), target).result();
        if (!read_result.ok()) {
            return Error{ErrorCode::io_error, "TensorStore read failed: " + read_result.status().ToString(),
                         std::string(node)};
        }
        return CheckReadControl(options, node);
    } catch (const std::exception& error) {
        return Error{ErrorCode::io_error, error.what(), std::string(node)};
    }
}

}  // namespace

Result<void> ReadFloat32(const std::filesystem::path& array_path, const StoreContextPtr& context,
                         std::string_view node, std::string_view expected_data_type, const PixelSelection& selection, float* destination,
                         std::size_t destination_elements, const ReadOptions& options) {
    return ReadInto(array_path, context, node, expected_data_type, selection, tensorstore::dtype_v<float>, destination,
                    destination_elements, options);
}

Result<void> ReadMaskBytes(const std::filesystem::path& array_path, const StoreContextPtr& context,
                           std::string_view node, std::string_view expected_data_type, const PixelSelection& selection, std::uint8_t* destination,
                           std::size_t destination_elements, const ReadOptions& options) {
    // The caller's buffer holds bytes, so the read converts into bytes. Asking TensorStore for
    // bool and writing it through a reinterpret_cast of that buffer assumed bool and uint8_t are
    // the same object, which C++ does not say they are; the conversion costs nothing here because
    // it rides the copy the read already performs.
    return ReadInto(array_path, context, node, expected_data_type, selection, tensorstore::dtype_v<std::uint8_t>,
                    destination, destination_elements, options);
}

}  // namespace carta::zarr::internal::zarr
