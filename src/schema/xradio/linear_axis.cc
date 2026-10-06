/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "linear_axis.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

// How close a sample must sit to the reference world value to be treated as landing on it. Scaled
// by magnitude so that it means the same thing for direction cosines and for hertz.
double ExactnessTolerance(const std::vector<double>& values) {
    double maximum = 1.0;
    for (const double value : values) {
        maximum = std::max(maximum, std::abs(value));
    }
    return 1.0e-12 * maximum;
}

std::size_t ClosestIndex(const std::vector<double>& values, double reference) {
    return static_cast<std::size_t>(std::distance(
        values.begin(), std::min_element(values.begin(), values.end(), [reference](double left, double right) {
            return std::abs(left - reference) < std::abs(right - reference);
        })));
}

// Whether the samples are evenly spaced. This is a looser test than exactness above: it asks whether
// consecutive spacings agree, not whether a value lands on a target.
bool IsEvenlySpaced(const std::vector<double>& values, double increment) {
    for (std::size_t index = 2; index < values.size(); ++index) {
        const double spacing = values.at(index) - values.at(index - 1);
        const double tolerance = 1.0e-9 * std::max({1.0, std::abs(increment), std::abs(spacing)});
        if (std::abs(spacing - increment) > tolerance) {
            return false;
        }
    }
    return true;
}

Diagnostic MakeDiagnostic(DiagnosticCode code, std::string message, std::string_view axis_name) {
    return Diagnostic{code, std::move(message), std::string(axis_name)};
}

}  // namespace

LinearAxisFit FitLinearAxis(const std::vector<double>& values, std::optional<double> reference,
                            std::string_view axis_name) {
    LinearAxisFit fit;
    if (values.size() < 2) {
        return fit;
    }

    // Checked before anything is fitted: a sample that is not a number fails every comparison, so
    // the spacing test below would pass it and describe the axis around it as linear.
    if (!std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); })) {
        fit.diagnostics.push_back(MakeDiagnostic(DiagnosticCode::nonuniform_axis,
                                                 "The " + std::string(axis_name) +
                                                     " coordinate has a sample that is not a finite number, so it "
                                                     "has no linear description",
                                                 axis_name));
        return fit;
    }

    double increment = values.at(1) - values.at(0);
    fit.increment = increment;
    if (increment == 0.0) {
        return fit;
    }

    fit.uniform = IsEvenlySpaced(values, increment);
    if (fit.uniform && values.size() > 2) {
        // Adjacent samples of a direction cosine axis differ in their last few digits, so one pair
        // carries several digits less precision than the axis itself does. Spreading the whole span
        // over the samples recovers them, and two axes cut from the same grid then agree exactly --
        // which casacore's square-pixel test requires.
        increment = (values.back() - values.front()) / static_cast<double>(values.size() - 1);
        fit.increment = increment;
    }
    if (!fit.uniform) {
        fit.diagnostics.push_back(
            MakeDiagnostic(DiagnosticCode::nonuniform_axis,
                           "The " + std::string(axis_name) +
                               " coordinate is not evenly spaced; any linear description of it is an approximation",
                           axis_name));
    }

    const double reference_value = reference.value_or(0.0);
    fit.reference_value = reference_value;

    const std::size_t closest = ClosestIndex(values, reference_value);
    if (std::abs(values.at(closest) - reference_value) <= ExactnessTolerance(values)) {
        fit.reference_pixel = static_cast<double>(closest + 1);
    } else {
        fit.reference_pixel = ((reference_value - values.front()) / increment) + 1.0;
        fit.diagnostics.push_back(MakeDiagnostic(DiagnosticCode::inexact_reference_pixel,
                                                 "The " + std::string(axis_name) +
                                                     " coordinate has no sample at its reference world "
                                                     "value; CRPIX was linearly extrapolated",
                                                 axis_name));
    }
    return fit;
}

DirectionAxisFit FitDirectionAxis(const std::vector<double>& cosines, std::string_view axis_name) {
    auto fit = FitLinearAxis(cosines, std::nullopt, axis_name);
    DirectionAxisFit direction;
    if (fit.increment) {
        direction.increment = *fit.increment * kRadToDeg;
    }
    direction.reference_pixel = fit.reference_pixel;
    direction.diagnostics = std::move(fit.diagnostics);
    // A direction axis is linear by construction, so an axis that cannot be described linearly is a
    // store this library cannot make sense of rather than a coordinate it reports tabularly -- which
    // is what a spectral axis in the same position is. Three ways to reach it: fewer than two
    // samples, two samples that are the same, and a sample that is not a finite number. Each leaves the
    // caller's reference pixel and increment at whatever they were, and says so.
    //
    // Diagnosed here rather than in FitLinearAxis because a spectral axis with one channel is not
    // degenerate, it is a continuum image, and it takes the tabular path by design.
    if (!fit.reference_pixel) {
        direction.diagnostics.push_back(
            MakeDiagnostic(DiagnosticCode::degenerate_axis,
                           "The " + std::string(axis_name) +
                               " coordinate has fewer than two distinct finite samples, so it has no linear "
                               "description; any reference pixel or increment reported for it is not usable",
                           axis_name));
    }
    return direction;
}

SpectralAxisFit FitSpectralAxis(const std::vector<double>& channels, double reference, std::string_view axis_name) {
    auto fit = FitLinearAxis(channels, reference, axis_name);
    SpectralAxisFit spectral;
    if (fit.uniform) {
        spectral.reference_pixel = fit.reference_pixel;
        spectral.reference_value = fit.reference_value;
        spectral.increment = fit.increment;
        spectral.diagnostics = std::move(fit.diagnostics);
        return spectral;
    }
    for (auto& diagnostic : fit.diagnostics) {
        if (diagnostic.code == DiagnosticCode::nonuniform_axis) {
            spectral.diagnostics.push_back(std::move(diagnostic));
        }
    }
    return spectral;
}

}  // namespace carta::zarr::internal::xradio
