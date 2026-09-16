/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_

#include "carta-zarr/carta_zarr.h"

#include <nlohmann/json.hpp>

#include <optional>
#include <vector>

namespace carta::zarr::internal::xradio {

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
std::optional<DirectionCoordinate> DescribeDirection(const nlohmann::json& root_attributes,
                                                     const std::vector<double>& l_values,
                                                     const std::vector<double>& m_values,
                                                     std::vector<Diagnostic>& diagnostics);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_DIRECTION_H_
