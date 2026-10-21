/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_

#include "../profile.h"
#include "coordinates.h"

#include <optional>
#include <string>
#include <vector>

namespace carta::zarr::internal::xradio {

// One enumeration of the store, answering both what it is and what is in it. See SchemaInspection.
Result<SchemaInspection> InspectImages(const Store& store);

// Reads the coordinate values, then describes the image from them. See DescribeImageFrom.
Result<DescribedImage> DescribeImage(const Store& store, std::string_view image_id);

// Every rule that turns an image's metadata and its coordinates' values into a descriptor, with the
// values handed in rather than read.
//
// Split from the reads by ADR 0006's rule: a module takes the narrowest input its tests can build.
// Reading a value goes to the filesystem through TensorStore, and the profile tests build without
// either. Given the values, they can test the paths that describe an image, not only those that
// refuse one.
Result<DescribedImage> DescribeImageFrom(const Store& store, std::string_view image_id, const CoordinateValues& values);
Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_IMAGE_H_
