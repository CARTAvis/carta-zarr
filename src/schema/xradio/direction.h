/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include <nlohmann/json.hpp>

#include <optional>
#include <vector>

namespace carta::zarr::internal::xradio {

// A sky direction measure as XRADIO writes one (SkyDirectionArray): two angles in radians under
// `data`, its frame and, for an equatorial frame, its equinox under `attrs`. Nothing when it has no
// two numeric angles. Both the reference direction and the pointing center are read with it, so the
// two name a frame the same way.
std::optional<SkyDirection> ReadSkyDirection(const nlohmann::json* measure);

/**
 * What the root's coordinate_system_info says: everything in a DirectionCoordinate but the two
 * things the l and m samples decide, the increment and the reference pixel.
 *
 * The one reader of that attribute: the probe uses it to decide whether to accept a store, and the
 * description to decide what to report, so the two cannot disagree.
 *
 * Every XRADIO writes it, so a store without it, or with one missing the projection, the reference
 * direction, the native pole or the matrix, is malformed: invalid_metadata, naming the attribute. The
 * reference frame, the equinox and the projection parameters may be absent, and leave their fields
 * empty.
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
 * Takes JSON and two vectors of samples and returns a descriptor, so everything it decides can be
 * tested without a store, a transport or a directory tree: converting the attributes' radians to the
 * descriptor's degrees, which slot of the transformation matrix each element goes in, and how an
 * equinox written "J2000" becomes a number.
 */
Result<DirectionCoordinate> DescribeDirection(const nlohmann::json& root_attributes,
                                              const std::vector<double>& l_values, const std::vector<double>& m_values,
                                              std::vector<Diagnostic>& diagnostics);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
