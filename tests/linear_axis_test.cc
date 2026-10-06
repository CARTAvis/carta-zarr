/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The linear axis fit takes a coordinate's samples and reports the linear description they support.
// It needs no store, so the cases that matter most -- unevenly spaced samples, a reference value no
// sample lands on, degenerate axes -- cost a vector literal each.

#include "schema/xradio/linear_axis.h"

#include "support/check.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

using carta::zarr::internal::xradio::FitDirectionAxis;
using carta::zarr::internal::xradio::FitLinearAxis;
using carta::zarr::internal::xradio::FitSpectralAxis;
using carta::zarr::internal::xradio::kRadToDeg;
using carta::zarr::internal::xradio::LinearAxisFit;

using carta::zarr::testing::Require;

bool Near(double left, double right) {
    return std::abs(left - right) <= 1.0e-9 * std::max({1.0, std::abs(left), std::abs(right)});
}

bool HasDiagnostic(const LinearAxisFit& fit, carta::zarr::DiagnosticCode code) {
    for (const auto& diagnostic : fit.diagnostics) {
        if (diagnostic.code == code) {
            return true;
        }
    }
    return false;
}

// Evenly spaced samples with one sitting exactly on the reference value: the reference pixel is that
// sample's 1-based index, and nothing is diagnosed.
void TestExactReferencePixel() {
    const auto fit = FitLinearAxis({100.0, 102.0, 104.0}, 100.0, "frequency");
    Require(fit.uniform, "evenly spaced samples were not reported as uniform");
    Require(fit.increment && Near(*fit.increment, 2.0), "increment was not the sample spacing");
    Require(fit.reference_value && Near(*fit.reference_value, 100.0), "reference value was not carried through");
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 1.0), "reference pixel was not the exact sample index");
    Require(fit.diagnostics.empty(), "an exact fit raised a diagnostic");

    const auto middle = FitLinearAxis({100.0, 102.0, 104.0}, 102.0, "frequency");
    Require(middle.reference_pixel && Near(*middle.reference_pixel, 2.0), "reference pixel was not 1-based");
}

// A direction cosine axis is measured from the tangent point, so a nullopt reference means 0.0.
void TestTangentPointReference() {
    const auto fit = FitLinearAxis({-0.003, -0.002, -0.001, 0.0}, std::nullopt, "l");
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 4.0), "the tangent point was not located");
    Require(fit.diagnostics.empty(), "locating the tangent point exactly raised a diagnostic");
}

// No sample lands on the reference value, so the reference pixel is extrapolated and said to be.
void TestInexactReferencePixel() {
    const auto fit = FitLinearAxis({-0.003, -0.002, -0.001, -0.0005}, std::nullopt, "l");
    Require(HasDiagnostic(fit, carta::zarr::DiagnosticCode::inexact_reference_pixel),
            "an extrapolated reference pixel was not diagnosed");
    // Spacing is 0.001 from the first pair, and 0.0 lies three increments past -0.003.
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 4.0), "reference pixel was not linearly extrapolated");
    Require(fit.uniform == false, "samples with a changing spacing were reported as uniform");
}

// Unevenly spaced samples still yield an increment -- a direction axis needs one -- but say so.
void TestNonUniformIsReported() {
    const auto fit = FitLinearAxis({1.4e9, 1.401e9, 1.403e9}, 1.4e9, "frequency");
    Require(!fit.uniform, "unevenly spaced samples were reported as uniform");
    Require(HasDiagnostic(fit, carta::zarr::DiagnosticCode::nonuniform_axis),
            "unevenly spaced samples were not diagnosed");
    Require(fit.increment && Near(*fit.increment, 1.0e6), "increment was not taken from the first pair");
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 1.0), "reference pixel was not located");
}

// Two samples are evenly spaced by definition; there is no second spacing to disagree with.
void TestTwoSamplesAreUniform() {
    const auto fit = FitLinearAxis({10.0, 20.0}, 10.0, "frequency");
    Require(fit.uniform, "a two-sample axis was not reported as uniform");
    Require(!HasDiagnostic(fit, carta::zarr::DiagnosticCode::nonuniform_axis),
            "a two-sample axis was diagnosed as uneven");
}

// Degenerate axes are not the unevenly sampled axis the diagnostic is about, so they stay silent.
void TestDegenerateAxesAreSilent() {
    const auto single = FitLinearAxis({0.5}, std::nullopt, "l");
    Require(!single.increment, "a single-sample axis reported an increment");
    Require(!single.reference_pixel, "a single-sample axis reported a reference pixel");
    Require(single.diagnostics.empty(), "a single-sample axis raised a diagnostic");
    Require(!single.uniform, "a single-sample axis claimed to be uniform");

    const auto empty = FitLinearAxis({}, std::nullopt, "l");
    Require(!empty.increment && empty.diagnostics.empty(), "an empty axis was not silent");

    // Repeated samples have a zero increment: an increment exists but no pixel can be located from
    // it, and dividing by it would be the only way to try.
    const auto repeated = FitLinearAxis({7.0, 7.0, 7.0}, 9.0, "frequency");
    Require(repeated.increment && Near(*repeated.increment, 0.0), "a zero increment was not reported");
    Require(!repeated.reference_pixel, "a zero increment produced a reference pixel");
    Require(repeated.diagnostics.empty(), "a zero-increment axis raised a diagnostic");
}

// A descending axis is ordinary: increments are signed.
void TestDescendingAxis() {
    const auto fit = FitLinearAxis({104.0, 102.0, 100.0}, 100.0, "frequency");
    Require(fit.uniform, "a descending axis was not reported as uniform");
    Require(fit.increment && Near(*fit.increment, -2.0), "a descending axis lost the sign of its increment");
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 3.0), "a descending axis mislocated its reference pixel");
}

// ---------------------------------------------------------------------------------------------
// What the two callers make of a fit.
//
// The fit reports; these decide, and they decide opposite things for opposite reasons. Each is
// reached here without a directory tree, and DirectionCoordinate::increment is asserted, so l and m
// cannot reach the consumer wrong by the factor of 57.3 between radians and degrees.

template <typename Fit>
bool Diagnosed(const Fit& fit, carta::zarr::DiagnosticCode code) {
    for (const auto& diagnostic : fit.diagnostics) {
        if (diagnostic.code == code) {
            return true;
        }
    }
    return false;
}

// The samples are direction cosines and the descriptor reports degrees.
void TestADirectionAxisReportsDegrees() {
    // A tenth of a milliradian per pixel, which is about 0.0057 degrees.
    const std::vector<double> cosines{-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4};
    const auto fit = FitDirectionAxis(cosines, "l");
    Require(fit.increment.has_value(), "a direction axis always has an increment");
    Require(Near(*fit.increment, 1.0e-4 * kRadToDeg),
            "the increment was not converted from radians to degrees: expected " + std::to_string(1.0e-4 * kRadToDeg) +
                ", got " + std::to_string(*fit.increment));
    Require(!Near(*fit.increment, 1.0e-4), "the increment was left in radians");
    Require(fit.reference_pixel.has_value() && Near(*fit.reference_pixel, 3.0),
            "the tangent point is the third sample, 1-based");
    Require(fit.diagnostics.empty(), "an evenly spaced direction axis is not worth a diagnostic");
}

// A direction axis is linear by construction, so it keeps its increment even when the samples come
// back uneven. It says so, and the consumer uses it anyway.
void TestADirectionAxisKeepsAnUnevenIncrement() {
    const std::vector<double> cosines{0.0, 1.0e-4, 2.5e-4, 3.0e-4};
    const auto fit = FitDirectionAxis(cosines, "m");
    Require(fit.increment.has_value(), "an uneven direction axis still reports an increment");
    Require(Near(*fit.increment, 1.0e-4 * kRadToDeg), "and it is the first pair's spacing, in degrees");
    Require(Diagnosed(fit, carta::zarr::DiagnosticCode::nonuniform_axis),
            "and it says the samples were not evenly spaced");
}

// A direction axis is linear by construction, so one that cannot be described linearly is a store
// this library cannot make sense of -- unlike a spectral axis in the same position, which is a
// continuum image taking the tabular path. DescribeDirection leaves the caller's reference pixel
// and increment at their defaults, and says so: a DirectionCoordinate reading 0 for both is not to
// be taken for an answer.
void TestADegenerateDirectionAxisSaysSo() {
    const auto one = FitDirectionAxis({1.0e-4}, "l");
    Require(!one.increment && !one.reference_pixel, "a single sample describes no axis");
    Require(Diagnosed(one, carta::zarr::DiagnosticCode::degenerate_axis), "and it should say so");

    Require(!FitDirectionAxis({}, "l").increment, "and neither does none");
    Require(Diagnosed(FitDirectionAxis({}, "l"), carta::zarr::DiagnosticCode::degenerate_axis),
            "which is the same answer");

    // Two samples that are the same: the increment is reported as the zero it is, but there is no
    // reference pixel to go with it, so this is the same non-answer wearing a number.
    const auto flat = FitDirectionAxis({1.0e-4, 1.0e-4}, "m");
    Require(!flat.reference_pixel, "a zero increment locates no reference pixel");
    Require(Diagnosed(flat, carta::zarr::DiagnosticCode::degenerate_axis), "and that is worth saying too");

    // An axis that is merely uneven is not degenerate: it has a linear description and says it is an
    // approximation. The two must not be confused.
    const auto uneven = FitDirectionAxis({0.0, 1.0e-4, 2.5e-4, 3.0e-4}, "l");
    Require(!Diagnosed(uneven, carta::zarr::DiagnosticCode::degenerate_axis),
            "an uneven axis is described, not refused");
}

// The opposite rule: unevenly spaced channels get no linear description at all, because a consumer
// that received one could not tell it was an approximation. See ADR 0002.
void TestASpectralAxisWithholdsWhatItCannotDescribe() {
    const std::vector<double> channels{1.0e9, 1.001e9, 1.0035e9, 1.004e9};
    const auto fit = FitSpectralAxis(channels, 1.0e9, "frequency");
    Require(!fit.increment && !fit.reference_pixel && !fit.reference_value,
            "an unevenly spaced spectral axis reports no linear description");
    Require(Diagnosed(fit, carta::zarr::DiagnosticCode::nonuniform_axis),
            "it keeps the reason, which says to build a tabular axis");
    Require(!Diagnosed(fit, carta::zarr::DiagnosticCode::inexact_reference_pixel),
            "and drops the one describing a reference pixel the consumer never receives");
}

// The diagnostic that is dropped above is genuinely raised underneath, so the dropping is a decision
// rather than a coincidence of this sample set.
void TestTheDroppedDiagnosticWasReallyThere() {
    const std::vector<double> channels{1.0e9, 1.001e9, 1.0035e9, 1.004e9};
    const auto underneath = FitLinearAxis(channels, 1.0005e9, "frequency");
    Require(HasDiagnostic(underneath, carta::zarr::DiagnosticCode::nonuniform_axis) &&
                HasDiagnostic(underneath, carta::zarr::DiagnosticCode::inexact_reference_pixel),
            "the fit underneath raises both, which is what makes the filter above a decision");
    const auto fit = FitSpectralAxis(channels, 1.0005e9, "frequency");
    Require(fit.diagnostics.size() == 1 && Diagnosed(fit, carta::zarr::DiagnosticCode::nonuniform_axis),
            "only one of them survives");
}

void TestAnEvenSpectralAxisKeepsItsDescription() {
    const std::vector<double> channels{1.0e9, 1.001e9, 1.002e9, 1.003e9};
    const auto fit = FitSpectralAxis(channels, 1.002e9, "frequency");
    Require(fit.increment && Near(*fit.increment, 1.0e6), "an even spectral axis reports its increment");
    Require(fit.reference_value && Near(*fit.reference_value, 1.002e9), "and the value it was given");
    Require(fit.reference_pixel && Near(*fit.reference_pixel, 3.0), "and the pixel that value lands on");
    Require(fit.diagnostics.empty(), "with nothing to report");
}

// A sample that is not a number fails every comparison, so a spacing test written as "refuse when
// it differs" passes it: [1e9, NaN, 1.002e9] came out evenly spaced at 1e6. Such an axis has no
// linear description, whichever kind it is, and a consumer told it had one would draw it.
void TestANonFiniteSampleDescribesNoLinearAxis() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto& channels : {std::vector<double>{1.0e9, nan, 1.002e9}, std::vector<double>{nan, 1.0e9, 1.001e9},
                                 std::vector<double>{1.0e9, 1.001e9, inf}}) {
        const auto underneath = FitLinearAxis(channels, 1.0e9, "frequency");
        Require(!underneath.uniform && !underneath.increment && !underneath.reference_pixel,
                "an axis with a non-finite sample was given a linear description");
        const auto spectral = FitSpectralAxis(channels, 1.0e9, "frequency");
        Require(!spectral.increment && !spectral.reference_pixel && !spectral.reference_value,
                "a spectral axis with a non-finite sample reported a linear description");
        Require(Diagnosed(spectral, carta::zarr::DiagnosticCode::nonuniform_axis), "and did not say why");
        const auto direction = FitDirectionAxis({0.0, nan, 2.0e-4}, "l");
        Require(!direction.increment && !direction.reference_pixel,
                "a direction axis with a non-finite sample reported a linear description");
        Require(Diagnosed(direction, carta::zarr::DiagnosticCode::degenerate_axis), "and did not say it is unusable");
    }
}

}  // namespace

int main() {
    try {
        TestExactReferencePixel();
        TestTangentPointReference();
        TestInexactReferencePixel();
        TestNonUniformIsReported();
        TestTwoSamplesAreUniform();
        TestDegenerateAxesAreSilent();
        TestDescendingAxis();
        TestADirectionAxisReportsDegrees();
        TestADirectionAxisKeepsAnUnevenIncrement();
        TestADegenerateDirectionAxisSaysSo();
        TestASpectralAxisWithholdsWhatItCannotDescribe();
        TestTheDroppedDiagnosticWasReallyThere();
        TestAnEvenSpectralAxisKeepsItsDescription();
        TestANonFiniteSampleDescribesNoLinearAxis();
        std::cout << "carta-zarr linear axis tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr linear axis tests failed: " << error.what() << '\n';
        return 1;
    }
}
