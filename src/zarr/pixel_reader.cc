/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
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
#include <tensorstore/strided_layout.h>
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
Result<void> ReadInto(const std::filesystem::path& array_path, const StoreContextPtr& context, std::string_view node,
                      const PixelSelection& selection, tensorstore::DataType target_dtype,
                      BufferView<Element> destination, const ReadControl& control) {
    if (destination.data == nullptr) {
        return Error{ErrorCode::invalid_argument, "Destination buffer is null", std::string(node)};
    }
    // The selection is well formed by construction, so it is not checked again. The destination is
    // whatever the caller passed across the seam, and the read writes through it as far as the
    // selection reaches, so this check is what stops a caller's mistake from writing past the end
    // of its buffer.
    if (selection.destination_extent() > destination.size) {
        return Error{ErrorCode::buffer_too_small, "Destination buffer is too small", std::string(node)};
    }

    if (!context) {
        return Error{ErrorCode::invalid_argument, "Pixel reads require a context", std::string(node)};
    }

    auto allowed = CheckReadControl(control, node);
    if (!allowed) {
        return allowed.error();
    }

    try {
        // The handle is opened once and reused: a slice read happens once per casacore cursor step.
        // A read with a pool of its own opens the array against that pool's context, because a
        // handle carries the pool it was opened against. The pool keys the handle by this store, so
        // a pool that outlives the store does not keep it.
        const StoreContext* pool =
            control.cache_pool ? CachePoolAccess::StoreContextOf(*control.cache_pool).get() : nullptr;
        auto opened = context->OpenArray(array_path, node, pool);
        if (!opened) {
            return opened.error();
        }
        allowed = CheckReadControl(control, node);
        if (!allowed) {
            return allowed.error();
        }
        auto const store = std::move(opened).value();
        // The store checked that this is the array it described before handing over where it lives,
        // including its rank, which matches the selection because the selection was built from that
        // array's descriptor.
        const auto rank = selection.start().size();
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
                         "Requested region is outside the array: " + sliced.status().ToString(), std::string(node)};
        }
        auto transposed = std::move(sliced).value() | tensorstore::Dims(order).Transpose();
        if (!transposed.ok()) {
            return Error{ErrorCode::invalid_argument, "Failed to reorder axes: " + transposed.status().ToString(),
                         std::string(node)};
        }

        // The conversion happens in the same copy, so a float64 or integer array is never
        // materialized in its stored type first.
        auto converted = tensorstore::Cast(std::move(transposed).value(), target_dtype);
        if (!converted.ok()) {
            return Error{ErrorCode::unsupported_data_type,
                         "Array cannot be converted to the requested output type: " + converted.status().ToString(),
                         std::string(node)};
        }

        // Axis 0 is the fastest-varying destination dimension: over a densely packed buffer the
        // strides are Fortran order's, and within a larger read they are that read's, from where this
        // part begins. TensorStore requires a shared array here. The caller owns this buffer and the
        // read below completes before returning, so a non-owning shared pointer describes the
        // ownership accurately.
        const auto strides = selection.DestinationStrides();
        std::vector<tensorstore::Index> shape(rank);
        std::vector<tensorstore::Index> byte_strides(rank);
        for (std::size_t i = 0; i < rank; ++i) {
            shape.at(i) = count.at(order.at(i));
            byte_strides.at(i) = static_cast<tensorstore::Index>(strides.at(order.at(i)) * sizeof(Element));
        }
        tensorstore::SharedArray<Element> target(
            tensorstore::internal::UnownedToShared(destination.data + selection.destination_offset()),
            tensorstore::StridedLayout<>(shape, byte_strides));

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

Result<void> ReadFloat32(const std::filesystem::path& array_path, const StoreContextPtr& context, std::string_view node,
                         const PixelSelection& selection, BufferView<float> destination, const ReadControl& control) {
    return ReadInto(array_path, context, node, selection, tensorstore::dtype_v<float>, destination, control);
}

Result<void> ReadMaskBytes(const std::filesystem::path& array_path, const StoreContextPtr& context,
                           std::string_view node, const PixelSelection& selection, BufferView<std::uint8_t> destination,
                           const ReadControl& control) {
    // The caller's buffer holds bytes, so the read converts into bytes rather than writing bool
    // through a reinterpret_cast, which C++ does not allow. The conversion costs nothing, because it
    // happens in the copy the read already performs.
    return ReadInto(array_path, context, node, selection, tensorstore::dtype_v<std::uint8_t>, destination, control);
}

}  // namespace carta::zarr::internal::zarr
