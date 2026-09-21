/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_
#define CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

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
    // Copied out of the document rather than left in it, the same way attributes is, because the
    // alternative is what this used to be: the parsed metadata and the raw document travelling
    // together to every caller that needed a field this struct did not carry, each of them reaching
    // back into the json and each deciding for itself what "absent" meant.
    //
    // Carried, not interpreted. What a chain means differs by who is asking -- a storage layout
    // wants to know whether it is sharded and what compressed it, a string decoder wants the order
    // the codecs run in and refuses the ones it cannot run -- so the chain is handed on whole and
    // each of them reads it its own way.
    //
    // An array with no codecs gets an empty array and no chunk key encoding gets an empty object,
    // so a reader asks what is in them rather than whether they are there.
    nlohmann::json codecs = nlohmann::json::array();
    nlohmann::json chunk_key_encoding = nlohmann::json::object();
};

Result<ArrayMetadata> ParseArrayMetadata(const nlohmann::json& metadata, std::string_view node);

// How an array is laid out on whatever is storing it: its chunk shape, whether the chunks are
// gathered into shards, and what compressed them.
//
// A projection of the document rather than a second reading of it, which is why it is here and not
// on Store. It was a Store method reaching into the raw node metadata, so the one thing it decides
// -- what a sharding codec says about the real chunk shape -- could only be exercised by writing a
// store to a directory and opening an image out of it. Nothing in it touches a transport.
//
// Reports invalid_metadata for a sharding codec whose chunk_shape is not positive integers, or does
// not have the rank of the shard it sits in.
Result<StorageLayout> ParseStorageLayout(const ArrayMetadata& metadata, std::string_view node);

bool IsNonNegativeInteger(const nlohmann::json& value);
bool IsPositiveInteger(const nlohmann::json& value);
bool IsNumericVector(const nlohmann::json& value, std::size_t length);
bool IsNumericMatrix(const nlohmann::json& value, std::size_t rows, std::size_t columns);
bool IsRealDataType(std::string_view data_type);
DataType ParseDataType(std::string_view data_type);
bool IsFixedLengthUtf32(const ArrayMetadata& metadata);
std::optional<std::size_t> FindDimensionIndex(const ArrayMetadata& metadata, std::string_view name);
const char* ErrorCodeName(ErrorCode code) noexcept;

} // namespace carta::zarr::internal::zarr

#endif // CARTA_ZARR_SRC_ZARR_ARRAY_METADATA_H_
