/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// How large a dataset is, which is a question about a dataset and not about a node. It was a method
// on Store, which is what made Store look like a dataset facade: it uses nothing private, holds no
// state, and has one caller.

#include "store.h"

#include <limits>
#include <map>
#include <string>
#include <string_view>

namespace carta::zarr::internal {
namespace {

Result<std::uint64_t> ElementSizeBytes(const zarr::ArrayMetadata& metadata, std::string_view node) {
    static const std::map<std::string_view, std::uint64_t> element_sizes{
        {"bool", 1},      {"int8", 1},       {"uint8", 1},      {"int16", 2},
        {"uint16", 2},    {"int32", 4},      {"uint32", 4},     {"int64", 8},
        {"uint64", 8},    {"float16", 2},    {"float32", 4},    {"float64", 8},
        {"complex64", 8}, {"complex128", 16},
    };
    if (const auto found = element_sizes.find(metadata.data_type); found != element_sizes.end()) {
        return found->second;
    }

    // XRADIO coordinate labels use the fixed_length_utf32 extension data type. Other fixed-length
    // extension types can be sized the same way when they declare length_bytes.
    if (metadata.data_type_configuration.is_object() &&
        metadata.data_type_configuration.contains("length_bytes")) {
        const auto& length_bytes = metadata.data_type_configuration.at("length_bytes");
        if (zarr::IsPositiveInteger(length_bytes)) {
            return length_bytes.get<std::uint64_t>();
        }
    }

    return Error{ErrorCode::unsupported_data_type,
                 "Array " + std::string(node) + " has unsupported data_type " + metadata.data_type,
                 std::string(node)};
}

}  // namespace

Result<std::uint64_t> TotalArraySizeBytes(const Store& store) {
    const auto& nodes_result = store.ListNodes();
    if (!nodes_result) {
        return nodes_result.error();
    }

    std::uint64_t total_bytes = 0;
    std::size_t array_count = 0;
    for (const auto& node : nodes_result.value()) {
        const auto& metadata = store.ReadNodeMetadata(node);
        if (!metadata || !metadata.value().is_object() || metadata.value().value("node_type", "") != "array") {
            continue;
        }
        ++array_count;

        const auto& array_metadata_result = store.ReadArrayMetadata(node);
        if (!array_metadata_result) {
            return array_metadata_result.error();
        }
        const auto& array_metadata = array_metadata_result.value();
        auto element_size_result = ElementSizeBytes(array_metadata, node);
        if (!element_size_result) {
            return element_size_result.error();
        }

        std::uint64_t array_bytes = element_size_result.value();
        for (const auto dimension : array_metadata.shape) {
            if (dimension == 0) {
                array_bytes = 0;
                break;
            }
            if (array_bytes > std::numeric_limits<std::uint64_t>::max() / dimension) {
                return Error{ErrorCode::invalid_metadata,
                             "Array " + node + " byte size overflows uint64_t", node};
            }
            array_bytes *= dimension;
        }
        if (total_bytes > std::numeric_limits<std::uint64_t>::max() - array_bytes) {
            return Error{ErrorCode::invalid_metadata, "Total Zarr array byte size overflows uint64_t"};
        }
        total_bytes += array_bytes;
    }

    if (array_count == 0) {
        return Error{ErrorCode::invalid_metadata, "Zarr store contains no arrays"};
    }
    return total_bytes;
}

}  // namespace carta::zarr::internal
