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
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <utility>

namespace carta::zarr::internal::zarr {
namespace {

template <typename Element>
Result<void> ReadInto(const std::filesystem::path& array_path, const StoreContextPtr& context,
                      std::string_view node, const PixelSelection& selection, tensorstore::DataType target_dtype,
                      BufferView<Element> destination, const ReadControl& control) {
    if (destination.data == nullptr) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is null", std::string(node)};
    }
    // The selection is well formed by construction, so nothing about it is checked again here. The
    // destination is not the selection's: it is whatever the caller handed across the seam, and the
    // read writes through it as far as the selection reaches. Holding the one to the other is what
    // stands between a caller's mistake and a write past the end of its buffer.
    const auto elements = selection.elements();
    if (elements > destination.size) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is too small", std::string(node)};
    }

    if (!context) {
        return Error{ErrorCode::invalid_argument, "Pixel reads require a context", std::string(node)};
    }

    auto allowed = CheckReadControl(control, node);
    if (!allowed) {
        return allowed.error();
    }

    try {
        // Reused across calls. A slice read happens once per casacore cursor step, so opening here
        // would make a fixed cost a per-call one.
        // A read with a pool of its own runs against that pool's context rather than the session's.
        // Chosen here rather than by the caller because this is where the array handle is taken,
        // and a handle carries the pool it was opened against.
        const StoreContextPtr& pool =
            control.cache_pool ? CachePoolAccess::StoreContextOf(*control.cache_pool) : context;
        auto opened = pool->OpenArray(array_path, node);
        if (!opened) {
            return opened.error();
        }
        allowed = CheckReadControl(control, node);
        if (!allowed) {
            return allowed.error();
        }
        auto const store = std::move(opened).value();
        // Whether this is the array the store described was settled before the store handed over
        // where it lives. What is the read's own is that the selection addresses an array of this
        // rank.
        const auto rank = selection.start().size();
        if (static_cast<std::size_t>(store.rank()) != rank) {
            return Error{ErrorCode::invalid_argument, "Selection rank does not match the array rank",
                         std::string(node)};
        }
        std::vector<tensorstore::Index> start(rank);
        std::vector<tensorstore::Index> count(rank);
        std::vector<tensorstore::Index> stride(rank);
        std::vector<tensorstore::DimensionIndex> order(rank);
        for (std::size_t i = 0; i < rank; ++i) {
            start.at(i) = static_cast<tensorstore::Index>(selection.start().at(i));
            count.at(i) = static_cast<tensorstore::Index>(selection.count().at(i));
            stride.at(i) = static_cast<tensorstore::Index>(selection.stride().at(i));
            order.at(i) = static_cast<tensorstore::DimensionIndex>(selection.destination_to_stored().at(i));
        }

        // Slice in stored order, then move the stored dimensions into the destination's order. Both
        // are index transforms, so TensorStore composes them into the one copy the read already
        // performs.
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
        auto target = tensorstore::Array(tensorstore::internal::UnownedToShared(destination.data), shape,
                                         tensorstore::fortran_order);

        auto const read_result = tensorstore::Read(std::move(converted).value(), target).result();
        if (!read_result.ok()) {
            return Error{ErrorCode::io_error, "TensorStore read failed: " + read_result.status().ToString(),
                         std::string(node)};
        }
        return CheckReadControl(control, node);
    } catch (const std::exception& error) {
        return Error{ErrorCode::io_error, error.what(), std::string(node)};
    }
}

}  // namespace

Result<void> ReadFloat32(const std::filesystem::path& array_path, const StoreContextPtr& context,
                         std::string_view node, const PixelSelection& selection, BufferView<float> destination,
                         const ReadControl& control) {
    return ReadInto(array_path, context, node, selection, tensorstore::dtype_v<float>, destination, control);
}

Result<void> ReadMaskBytes(const std::filesystem::path& array_path, const StoreContextPtr& context,
                           std::string_view node, const PixelSelection& selection,
                           BufferView<std::uint8_t> destination, const ReadControl& control) {
    // The caller's buffer holds bytes, so the read converts into bytes. Asking TensorStore for
    // bool and writing it through a reinterpret_cast of that buffer assumed bool and uint8_t are
    // the same object, which C++ does not say they are; the conversion costs nothing here because
    // it rides the copy the read already performs.
    return ReadInto(array_path, context, node, selection, tensorstore::dtype_v<std::uint8_t>, destination,
                    control);
}

}  // namespace carta::zarr::internal::zarr
