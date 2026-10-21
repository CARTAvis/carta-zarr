/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_
#define CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include "zarr/storage_layout.h"

#include <nlohmann/json.hpp>

#include <optional>

namespace carta::zarr::internal::zarr {

struct ArrayMetadata {
    std::vector<std::uint64_t> shape;
    std::vector<std::string> dimension_names;
    std::string data_type;
    nlohmann::json data_type_configuration;
    nlohmann::json attributes;
    std::vector<std::uint64_t> chunk_shape;
    // Copied out of the document, as attributes is, so that no caller reaches back into the JSON and
    // decides for itself what "absent" means.
    //
    // Carried, not interpreted, because each reader wants something different: a storage layout
    // wants to know whether it is sharded and what compressed it, and a string decoder wants the
    // order the codecs run in and refuses those it cannot run.
    //
    // An array with no codecs gets an empty array, and one with no chunk key encoding an empty
    // object, so a reader asks what is in them rather than whether they are there.
    //
    // One check looks inside: ParseArrayMetadata refuses a `sharding_indexed` codec whose inner chunk
    // shape is not positive extents of the array's rank, the same rule it applies to `shape`,
    // `dimension_names` and the outer `chunk_grid`. It is checked here rather than in
    // ParseStorageLayout, which cannot fail, so that such an array never opens reporting every chunk
    // as the whole image.
    nlohmann::json codecs = nlohmann::json::array();
    nlohmann::json chunk_key_encoding = nlohmann::json::object();
    // What a chunk never written reads as, carried as written for the same reason as the codecs: its
    // form depends on the data type, and only a reader of that type knows what it means. TensorStore
    // reads it itself for every array it opens; the string decoder is the one reader here that asks
    // for it. Null when the document has none.
    nlohmann::json fill_value;
};

Result<ArrayMetadata> ParseArrayMetadata(const nlohmann::json& metadata, std::string_view node);

// How an array is laid out on whatever is storing it: its chunk shape, whether the chunks are
// gathered into shards, and what compressed them.
//
// Read from metadata that ParseArrayMetadata has already accepted, so it cannot fail: everything it
// could refuse was refused there. It is here rather than on Store so that what a sharding codec says
// about the real chunk shape can be tested without a store.
StorageLayout ParseStorageLayout(const ArrayMetadata& metadata);

bool IsNonNegativeInteger(const nlohmann::json& value);
bool IsPositiveInteger(const nlohmann::json& value);
bool IsNumericVector(const nlohmann::json& value, std::size_t length);
bool IsNumericMatrix(const nlohmann::json& value, std::size_t rows, std::size_t columns);
bool IsRealDataType(std::string_view data_type);
DataType ParseDataType(std::string_view data_type);
bool IsFixedLengthUtf32(const ArrayMetadata& metadata);
std::optional<std::size_t> FindDimensionIndex(const ArrayMetadata& metadata, std::string_view name);
const char* ErrorCodeName(ErrorCode code) noexcept;

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_
