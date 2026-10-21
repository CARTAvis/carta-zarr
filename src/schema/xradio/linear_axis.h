/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_LINEAR_AXIS_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_LINEAR_AXIS_H_

#include "carta-zarr/descriptor.h"

#include <cmath>
#include <optional>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

// The linear description of a sampled coordinate: the reference pixel, reference value, and
// increment a consumer needs to build a linear axis. Absent fields mean the samples do not support
// one; the consumer must build a tabular axis from the values instead.
struct LinearAxisFit {
    std::optional<double> reference_pixel;  // CRPIX, 1-based
    std::optional<double> reference_value;  // CRVAL, in the coordinate's own unit
    std::optional<double> increment;        // CDELT, in the coordinate's own unit
    bool uniform = false;
    std::vector<Diagnostic> diagnostics;
};

/**
 * Fit a linear description to a coordinate's samples.
 *
 * `reference` is the world value the reference pixel should land on; nullopt asks for the tangent
 * point at 0.0, which is what a direction cosine axis is measured from. `axis_name` names the axis
 * in any diagnostic raised.
 *
 * The fit reports rather than decides. A direction axis is linear by construction, so its caller
 * uses the increment even when `uniform` is false; a spectral axis need not be evenly spaced, so its
 * caller withholds the whole linear description in that case. Both surface the diagnostics.
 *
 * Values are returned in the unit they arrived in. A caller reporting a different unit scales the
 * increment afterwards, which is safe because the reference pixel is derived from the raw samples.
 *
 * A degenerate axis -- fewer than two samples, or a zero increment -- yields the increment where one
 * exists, no reference pixel, and no diagnostic. It is not evenly spaced, but neither is it the
 * unevenly sampled axis the diagnostic is about.
 */
LinearAxisFit FitLinearAxis(const std::vector<double>& values, std::optional<double> reference,
                            std::string_view axis_name);

inline constexpr double kRadToDeg = 180.0 / M_PI;

// What a direction axis makes of a fit.
struct DirectionAxisFit {
    std::optional<double> increment;  // CDELT, in degrees
    std::optional<double> reference_pixel;
    std::vector<Diagnostic> diagnostics;
};

/**
 * The linear description of one direction axis, from its direction cosines.
 *
 * A direction axis is linear by construction, so it keeps the increment whether or not the samples
 * are evenly spaced: uneven samples are rounding, not a sign that the axis is tabular. The diagnostic
 * saying so is kept either way.
 *
 * The samples are direction cosines and the descriptor reports degrees, so the increment is scaled
 * here, where tests can reach it. Without the scaling, an l/m pixel scale would be wrong by a factor
 * of 57.3.
 */
DirectionAxisFit FitDirectionAxis(const std::vector<double>& cosines, std::string_view axis_name);

/**
 * The linear description of one aperture-plane axis, u or v, from its samples in wavelengths.
 *
 * Linear by construction, as a direction axis is (it is a Fourier transform's grid), so the increment
 * is kept whether or not the samples are evenly spaced, and a degenerate axis is reported as one.
 * `reference` is the value the reference pixel lands on, XRADIO's `crval`; nullopt measures from 0.0,
 * the center of the plane. Nothing is scaled: the unit is the coordinate's own.
 */
DirectionAxisFit FitApertureAxis(const std::vector<double>& values, std::optional<double> reference,
                                 std::string_view axis_name);

// What a spectral axis makes of a fit. Every field is absent unless the channels are evenly spaced.
struct SpectralAxisFit {
    std::optional<double> reference_pixel;
    std::optional<double> reference_value;
    std::optional<double> increment;
    std::vector<Diagnostic> diagnostics;
};

/**
 * The linear description of a spectral coordinate, from its channel frequencies.
 *
 * The opposite rule to a direction axis, and for the opposite reason: XRADIO warns that neighboring
 * channels need not be evenly spaced, so unevenly spaced ones get no linear description at all and
 * the consumer builds a tabular axis instead. See ADR 0002.
 *
 * When the description is withheld, a diagnostic about its reference pixel would only confuse the
 * consumer, so it is dropped. The nonuniform_axis diagnostic is kept, because it tells the consumer
 * to build a tabular axis.
 */
SpectralAxisFit FitSpectralAxis(const std::vector<double>& channels, double reference, std::string_view axis_name);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_LINEAR_AXIS_H_
