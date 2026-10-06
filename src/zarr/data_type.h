/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_ZARR_DATA_TYPE_H_
#define CARTA_ZARR_SRC_ZARR_DATA_TYPE_H_

// What this library knows about a Zarr v3 data type, in one table.
//
// It was in five: the real types a coordinate may use, the name-to-enum map a descriptor is filled
// from, the byte sizes a dataset's size is summed with, the name-to-TensorStore-dtype check a pixel
// read makes, and a second byte-size switch keyed on the enum rather than the name. Each covered a
// different subset -- eleven types, twelve, fourteen, twelve, and "everything else is four bytes"
// -- so adding a type meant finding all five and knowing which of them cared.
//
// Header-only, and it includes no JSON: chunk_blocks.h reads it, and that header is compiled by a
// test target that links nothing at all.

#include "carta-zarr/descriptor.h"

#include <array>
#include <cstdint>
#include <string_view>

namespace carta::zarr::internal::zarr {

struct DataTypeInfo {
    // What a zarr.json spells.
    std::string_view name;
    // What a descriptor reports.
    DataType kind;
    // What one element occupies once decoded.
    std::uint64_t element_bytes;
    // Whether this is a real-valued type, which is what a coordinate array and an openable image
    // are each required to be. bool is not one: a boolean variable is a flag, never an image.
    bool real;
};

// Every data type this library names. `fixed_length_utf32` is deliberately absent: it is an
// extension type whose width is declared per array in its own configuration, so it is sized from
// that rather than from here.
inline constexpr std::array<DataTypeInfo, 14> kDataTypes{{
    {"bool", DataType::boolean, 1, false},
    {"int8", DataType::int8, 1, true},
    {"uint8", DataType::uint8, 1, true},
    {"int16", DataType::int16, 2, true},
    {"uint16", DataType::uint16, 2, true},
    {"int32", DataType::int32, 4, true},
    {"uint32", DataType::uint32, 4, true},
    {"int64", DataType::int64, 8, true},
    {"uint64", DataType::uint64, 8, true},
    {"float16", DataType::float16, 2, true},
    {"float32", DataType::float32, 4, true},
    {"float64", DataType::float64, 8, true},
    {"complex64", DataType::complex64, 8, false},
    {"complex128", DataType::complex128, 16, false},
}};

// The entry for a name as a zarr.json spells it, or nullptr for one this library does not name.
constexpr const DataTypeInfo* FindDataType(std::string_view name) noexcept {
    for (const auto& info : kDataTypes) {
        if (info.name == name) {
            return &info;
        }
    }
    return nullptr;
}

// The same by the enum a descriptor carries. DataType::unknown has no entry, which is what it
// means.
constexpr const DataTypeInfo* FindDataType(DataType kind) noexcept {
    for (const auto& info : kDataTypes) {
        if (info.kind == kind) {
            return &info;
        }
    }
    return nullptr;
}

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_DATA_TYPE_H_
