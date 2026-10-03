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
 * The array is opened from its own document, which the store has already held to the one it
 * described the image from (Store::VerifyArray); nothing here asks again.
 *
 * Missing chunks resolve to the array's fill value, which is what a Zarr reader is required to do
 * and is the only definition of "absent pixel" the format offers.
 */
Result<void> ReadFloat32(const std::filesystem::path& array_directory, const StoreContextPtr& context,
                         std::string_view node, const PixelSelection& selection, BufferView<float> destination,
                         const ReadControl& control);

// Read a boolean array as one byte per element, true meaning a good pixel.
Result<void> ReadMaskBytes(const std::filesystem::path& array_directory, const StoreContextPtr& context,
                           std::string_view node, const PixelSelection& selection,
                           BufferView<std::uint8_t> destination, const ReadControl& control);

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_PIXEL_READER_H_
