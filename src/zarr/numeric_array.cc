/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "numeric_array.h"

#include <limits>
#include <string>

namespace carta::zarr::internal::zarr {

Result<NumericArray> NumericArray::Make(std::string node, ArrayMetadata metadata, std::vector<double> values) {
    // Counted with a check rather than multiplied, so a shape whose product wraps cannot come out
    // equal to the buffer's length by accident.
    std::uint64_t declared = 1;
    for (const auto extent : metadata.shape) {
        if (extent != 0 && declared > std::numeric_limits<std::uint64_t>::max() / extent) {
            return Error{ErrorCode::invalid_metadata, "Array shape declares more values than can be counted", node};
        }
        declared *= extent;
    }
    if (values.size() != declared) {
        return Error{ErrorCode::invalid_metadata,
                     "Array holds " + std::to_string(values.size()) + " values where its shape declares " +
                         std::to_string(declared),
                     node};
    }
    return NumericArray(std::move(node), std::move(metadata), std::move(values));
}

Result<double> NumericArray::At(const std::vector<NamedIndex>& indices) const {
    const std::size_t rank = _metadata.shape.size();
    if (_metadata.dimension_names.size() != rank) {
        return Error{ErrorCode::invalid_metadata, "Array cannot be addressed by name without dimension names", _node};
    }

    // C order: the last dimension varies fastest.
    std::vector<std::uint64_t> strides(rank, 1);
    for (std::size_t dimension = rank; dimension-- > 1;) {
        strides.at(dimension - 1) = strides.at(dimension) * _metadata.shape.at(dimension);
    }

    std::uint64_t offset = 0;
    for (const auto& [name, index] : indices) {
        const auto dimension = FindDimensionIndex(_metadata, name);
        if (!dimension) {
            return Error{ErrorCode::invalid_slice, "Array has no dimension named " + std::string(name), _node};
        }
        if (index >= _metadata.shape.at(*dimension)) {
            return Error{ErrorCode::invalid_slice,
                         "Index " + std::to_string(index) + " is past the end of dimension " + std::string(name),
                         _node};
        }
        offset += index * strides.at(*dimension);
    }
    // Inside the buffer: every index is inside its dimension, and Make held the buffer to the shape.
    return _values.at(offset);
}

}  // namespace carta::zarr::internal::zarr
