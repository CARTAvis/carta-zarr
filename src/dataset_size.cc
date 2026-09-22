/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// How large a dataset is, which is a question about a dataset and not about a node. It was a method
// on Store, which is what made Store look like a dataset facade: it uses nothing private, holds no
// state, and has one caller.
//
// Both halves of the answer live here now. The facade used to keep the measured half -- it parsed
// the location string a second time and walked the directory itself -- which meant the rule for
// what a location means existed twice, and the copy up there was the weaker one: it never resolved
// the path, so a process that changed directory after opening the store measured somewhere else, or
// nowhere, and silently fell back to the declared size instead.

#include "store.h"

#include "zarr/data_type.h"

#include <chrono>
#include <limits>
#include <string>
#include <string_view>

namespace carta::zarr::internal {
namespace {

Result<std::uint64_t> ElementSizeBytes(const zarr::ArrayMetadata& metadata, std::string_view node) {
    if (const auto* const info = zarr::FindDataType(metadata.data_type); info != nullptr) {
        return info->element_bytes;
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

// Every array in the hierarchy, at its uncompressed size: what the metadata declares the dataset
// holds, and what a caller gets when the store itself cannot be measured.
//
// Not a bound on what the store occupies, in either direction. It counts array data only, so the
// per-node zarr.json documents and any sharding indices are missing from it -- in this repo's own
// fixtures those outweigh the compressed chunks several times over -- while compression pushes the
// other way. Which wins is a property of the store, and this function runs when the store could not
// be read. See ADR 0008.

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

}  // namespace

Result<DatasetSize> DatasetSizeBytes(const Store& store, std::chrono::steady_clock::time_point deadline) {
    if (auto stored = store.StoredSizeBytes(deadline)) {
        return DatasetSize{stored.value(), SizeBasis::measured};
    }
    // Every way of failing to measure leaves the same thing to answer with: the transport has no
    // bytes to count, the deadline passed part-way through, a directory refused to be enumerated,
    // or the total overflowed. The declared size is a more useful answer than an error to all four,
    // and is what this reported before there was an error to swallow -- the walk was a bool and
    // every one of its failure paths returned false.
    //
    // What is lost on the way is which of the four it was, and that is what would be needed to say
    // anything about the number that comes back: a deadline that expired says the store is large,
    // so the declared size is almost certainly above what it occupies, while a directory that
    // refused to be read says nothing at all. The caller is told which question was answered, not
    // how the answers compare.
    auto declared = TotalArraySizeBytes(store);
    if (!declared) {
        return declared.error();
    }
    return DatasetSize{declared.value(), SizeBasis::declared};
}

}  // namespace carta::zarr::internal
