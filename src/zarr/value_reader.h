/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_VALUE_READER_H_
#define CARTA_ZARR_SRC_ZARR_VALUE_READER_H_

#include "store.h"

#include <filesystem>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::zarr {

/**
 * Read a numeric Zarr array as doubles through TensorStore.
 *
 * One of the two translation units that read array data through TensorStore -- pixels are the
 * other -- so a consumer that needs metadata but never values, a schema profile and its tests, can
 * link without either. A null context uses TensorStore's own default resources.
 *
 * Nothing here asks whether the array is the one the store described: the store has, before it
 * hands over the location. See Store::VerifyArray.
 */
Result<std::vector<double>> ReadNumericValues(const std::filesystem::path& array_directory,
                                              const StoreContextPtr& context, std::string_view node);

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_VALUE_READER_H_
