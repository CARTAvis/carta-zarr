/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_STRING_ARRAY_H_
#define CARTA_ZARR_SRC_ZARR_STRING_ARRAY_H_

#include "carta-zarr/result.h"

#include "array_metadata.h"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::zarr {

/**
 * Decode a 1-D Zarr v3 "fixed_length_utf32" string array by hand.
 *
 * TensorStore's zarr3 driver does not support string data types, so the chunk is read and decoded
 * directly. 1-D arrays of any number of chunks using the default or v2 chunk key encoding are
 * supported, with a codec chain of: an optional identity transpose, bytes, at most one of
 * zstd/gzip/blosc, and any number of crc32c codecs (checksums are verified). This covers the layouts
 * produced by XRADIO for coordinate label arrays. A missing chunk yields the empty fill value for
 * each of its elements.
 */
Result<std::vector<std::string>> ReadFixedLengthUtf32StringArray(const std::filesystem::path& array_directory,
                                                                 const ArrayMetadata& array_metadata,
                                                                 std::string_view node);

/**
 * Whether ReadFixedLengthUtf32StringArray would decode an array with this metadata, asked of the
 * metadata alone: its type, its rank, its chunk key encoding and its codec chain. Fails with the
 * error the read would, so that a probe refuses what an open would rather than listing an image
 * that cannot be opened. A chunk that is missing or corrupt is not a metadata question and is left
 * to the read.
 */
Result<void> CheckFixedLengthUtf32StringArray(const ArrayMetadata& array_metadata, std::string_view node);

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_STRING_ARRAY_H_
