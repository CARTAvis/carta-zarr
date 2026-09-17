/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_

#include "../profile.h"

namespace carta::zarr::internal::xradio {

// One enumeration of the store, answering both what it is and what is in it. See SchemaInspection.
Result<SchemaInspection> InspectImages(const Store& store);
Result<ImageDescriptor> DescribeImage(const Store& store, std::string_view image_id);
Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_
