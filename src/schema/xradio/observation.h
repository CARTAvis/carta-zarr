/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_OBSERVATION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_OBSERVATION_H_

#include "carta-zarr/descriptor.h"

#include "../../zarr/array_metadata.h"

namespace carta::zarr::internal::xradio {

/**
 * What an image variable's own attributes say about the observation that produced it.
 *
 * Optional throughout, and deliberately silent: a field written in a type this cannot read is
 * skipped and the rest is kept. Nothing here can close an image, which is why it returns a value
 * rather than a Result and reports no diagnostic. A telescope position that is half unreadable is
 * left empty rather than half converted -- every element is checked before it is converted, because
 * the alternative is an exception out of nlohmann on a path whose caller is holding a Result.
 *
 * Takes one array's metadata and returns a descriptor field, so the geodetic-to-cartesian
 * conversion and the obsdate scale/value split are checkable without a store, a transport or a
 * directory tree.
 */
ObservationInfo DescribeObservation(const zarr::ArrayMetadata& image);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_OBSERVATION_H_
