/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// How large a dataset is, which is a question about a dataset and not about a node, so it is not a
// method on Store. Both halves of the answer, measured and declared, live here, and the measured
// one walks the location the store resolved when it opened rather than parsing the string again:
// a process that changed directory since would otherwise measure somewhere else.

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
    if (metadata.data_type_configuration.is_object() && metadata.data_type_configuration.contains("length_bytes")) {
        const auto& length_bytes = metadata.data_type_configuration.at("length_bytes");
        if (zarr::IsPositiveInteger(length_bytes)) {
            return length_bytes.get<std::uint64_t>();
        }
    }

    return Error{ErrorCode::unsupported_data_type,
                 "Array " + std::string(node) + " has unsupported data_type " + metadata.data_type, std::string(node)};
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
    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }

    std::uint64_t total_bytes = 0;
    std::size_t array_count = 0;
    for (const auto& entry : inventory.value()) {
        // A group holds no array data, and neither does a node that does not say it is an array --
        // including one whose document would not parse, which might have been an array and cannot
        // be sized either way. Leaving it out makes the total smaller than it might be, which a
        // declared size already does not claim to rule out (ADR 0008); refusing would leave a
        // dataset that opens with no size at all.
        if (entry.kind != NodeKind::array) {
            continue;
        }
        ++array_count;

        // An array this total knows is there and cannot size makes the total wrong rather than
        // smaller, so it refuses instead of leaving the array out.
        const auto& node = entry.name;
        const auto& array_metadata_result = *entry.array;
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
                return Error{ErrorCode::invalid_metadata, "Array " + node + " byte size overflows uint64_t", node};
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
