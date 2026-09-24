/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include <nlohmann/json.hpp>

#include <optional>
#include <vector>

namespace carta::zarr::internal::xradio {

/**
 * What the root's coordinate_system_info says: everything in a DirectionCoordinate but the two
 * things the l and m samples decide, the increment and the reference pixel.
 *
 * The one reading of that attribute. The probe asks it whether to accept a store and the description
 * asks it what to report, and each used to read the attribute its own way: the probe refused a
 * projection that was missing, a pole that was not two numbers, a matrix that was not two by two,
 * while the description filled in a default for each -- a leniency nothing could reach, because the
 * probe had already refused. And the probe looked only when the attribute was there, so a store
 * without one opened with a projection of "" and a reference of (0, 0), which read as answers.
 *
 * Every XRADIO writes it, so a store without it, or with one missing any of those four, is
 * malformed: invalid_metadata, naming the attribute. What may be absent is the reference frame, the
 * equinox and the projection parameters, which leave their fields empty.
 */
Result<DirectionCoordinate> ReadCoordinateSystem(const nlohmann::json& root_attributes);

/**
 * The direction coordinate an image dataset's root attributes and its l and m samples describe.
 *
 * Everything casacore's DirectionCoordinate constructor takes, which is why it is more than a
 * projection and a reference value: the projection parameters and the native pole become longPole
 * and latPole, and omitting them shifts coordinates silently rather than raising an error. See
 * ADR 0002.
 *
 * Takes JSON and two vectors of samples and returns a descriptor, so everything it decides -- the
 * radians the attributes are written in against the degrees the descriptor reports, which slot of
 * the transformation matrix each element belongs in, how an equinox spelled "J2000" becomes a
 * number -- is checkable without a store, a transport or a directory tree. It was file-local inside
 * a 779-line translation unit before, and none of those was.
 */
Result<DirectionCoordinate> DescribeDirection(const nlohmann::json& root_attributes,
                                              const std::vector<double>& l_values,
                                              const std::vector<double>& m_values,
                                              std::vector<Diagnostic>& diagnostics);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
