/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_PIXEL_READER_H_
#define CARTA_ZARR_SRC_ZARR_PIXEL_READER_H_

#include "store.h"
#include "zarr/pixel_selection.h"

#include <cstdint>
#include <filesystem>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::zarr {

/**
 * Read pixels as float32, converting from the stored type during the read.
 *
 * `expected` is the array's metadata as the store parsed it. The array on disk is held to it --
 * rank, extent, dimension names and data type -- before anything is read.
 *
 * Missing chunks resolve to the array's fill value, which is what a Zarr reader is required to do
 * and is the only definition of "absent pixel" the format offers.
 */
Result<void> ReadFloat32(const std::filesystem::path& array_directory, const StoreContextPtr& context,
                         std::string_view node, const ArrayMetadata& expected, const PixelSelection& selection, float* destination,
                         std::size_t destination_elements, const ReadControl& control);

// Read a boolean array as one byte per element, true meaning a good pixel.
Result<void> ReadMaskBytes(const std::filesystem::path& array_directory, const StoreContextPtr& context,
                           std::string_view node, const ArrayMetadata& expected, const PixelSelection& selection, std::uint8_t* destination,
                           std::size_t destination_elements, const ReadControl& control);

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_PIXEL_READER_H_
