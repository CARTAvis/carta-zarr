/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_ZARR_NUMERIC_ARRAY_H_
#define CARTA_ZARR_SRC_ZARR_NUMERIC_ARRAY_H_

#include "carta-zarr/result.h"

#include "array_metadata.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace carta::zarr::internal::zarr {

/**
 * A numeric Zarr array's values together with the document they were decoded with, addressed by
 * dimension name rather than by offset.
 *
 * Values are read as one flat buffer in C order, so a caller wanting one element of a
 * multi-dimensional array had to derive the strides itself -- which means knowing the array's stored
 * order, the very thing a schema profile should not have to track. And it had to pair the buffer
 * with metadata it fetched separately, which with consolidated metadata was the root's copy rather
 * than the document the values came from. Here the two arrive together and cannot be paired wrongly:
 * Store::ReadNumericArray hands back the array's own document with its values.
 *
 * Holds exactly as many values as its shape declares; Make refuses anything else, so an element
 * the shape promises is always there. A dimension the caller does not name is taken at index 0.
 * Naming a dimension the array does not have, or an index past the end of one it does, is an error
 * rather than a value.
 */
class NumericArray {
public:
    using NamedIndex = std::pair<std::string_view, std::uint64_t>;

    // Reports invalid_metadata, naming `node`, when the values are not as many as the shape declares.
    static Result<NumericArray> Make(std::string node, ArrayMetadata metadata, std::vector<double> values);

    // The node the values were read from, so an error about them can say which variable it is about.
    const std::string& node() const noexcept { return _node; }
    const ArrayMetadata& metadata() const noexcept { return _metadata; }
    // In C order, flattened.
    const std::vector<double>& values() const noexcept { return _values; }

    // Takes a vector so a caller can name a dimension conditionally; a braced list still works.
    Result<double> At(const std::vector<NamedIndex>& indices) const;

private:
    NumericArray(std::string node, ArrayMetadata metadata, std::vector<double> values)
        : _node(std::move(node)), _metadata(std::move(metadata)), _values(std::move(values)) {}

    std::string _node;
    ArrayMetadata _metadata;
    std::vector<double> _values;
};

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_NUMERIC_ARRAY_H_
