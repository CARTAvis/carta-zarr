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
    // Copied out of the document rather than left in it, the same way attributes is, so that no
    // caller reaches back into the json and decides for itself what "absent" means.
    //
    // Carried, not interpreted. What a chain means differs by who is asking -- a storage layout
    // wants to know whether it is sharded and what compressed it, a string decoder wants the order
    // the codecs run in and refuses the ones it cannot run -- so the chain is handed on whole and
    // each of them reads it its own way.
    //
    // An array with no codecs gets an empty array and no chunk key encoding gets an empty object,
    // so a reader asks what is in them rather than whether they are there.
    //
    // ParseArrayMetadata is the one exception, and it is not an interpretation: it looks inside a
    // `sharding_indexed` codec far enough to refuse an inner chunk shape that is not a positive
    // extent of the array's rank. That is the same well-formedness rule it applies to `shape`,
    // `dimension_names` and the outer `chunk_grid`, asked at the same gate, and it keeps nothing.
    //
    // Checked here rather than in ParseStorageLayout, which is infallible: an array with a malformed
    // sharding codec must not open and report a chunk geometry synthesised from its own shape -- an
    // image that says every one of its chunks is the whole image. Checking at the gate makes that
    // unreachable rather than guarded against.
    nlohmann::json codecs = nlohmann::json::array();
    nlohmann::json chunk_key_encoding = nlohmann::json::object();
    // What a chunk never written reads as, carried as written for the same reason the codecs are:
    // its spelling depends on the data type, and only a reader of that type can say what it means.
    // TensorStore reads it for itself for every array it opens; the string decoder, which TensorStore
    // cannot stand in for, is the one reader here that asks. Null when the document has none.
    nlohmann::json fill_value;
};

Result<ArrayMetadata> ParseArrayMetadata(const nlohmann::json& metadata, std::string_view node);

// How an array is laid out on whatever is storing it: its chunk shape, whether the chunks are
// gathered into shards, and what compressed them.
//
// A projection of the document rather than a second reading of it, which is why it is here and not
// on Store: what a sharding codec says about the real chunk shape can be tested without a store.
// Nothing in it touches a transport.
//
// Reports invalid_metadata for a sharding codec whose chunk_shape is not positive integers, or does
// not have the rank of the shard it sits in.
// How the array is stored, read out of metadata ParseArrayMetadata has already accepted. Infallible
// for that reason: everything it could have refused is refused at the gate.
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
