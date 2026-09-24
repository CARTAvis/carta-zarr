/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What an image dataset's root attributes and its l and m samples say about where the image points.
//
// ADR 0002 rests on all of it reaching casacore intact, and until this file existed almost none of
// it was checked. Two gaps in particular: DirectionCoordinate::increment was asserted nowhere, and
// the transformation matrix was asserted nowhere while every fixture in the repository writes the
// identity -- so the four element copies could be transposed and nothing would notice.
//
// None of it needs a store: this takes JSON and two vectors.

#include "schema/xradio/direction.h"

#include "schema/xradio/linear_axis.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::Diagnostic;
using carta::zarr::DirectionCoordinate;
using carta::zarr::internal::xradio::DescribeDirection;
using carta::zarr::internal::xradio::kRadToDeg;

using carta::zarr::testing::Require;

bool Near(double left, double right) {
    return std::abs(left - right) <= 1.0e-9 * std::max({1.0, std::abs(left), std::abs(right)});
}

// A tenth of a milliradian per pixel, with the tangent point on the third sample.
const std::vector<double> kSamples{-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4};

nlohmann::json Root(nlohmann::json coordinate_system) {
    return nlohmann::json{{"type", "image_dataset"}, {"coordinate_system_info", std::move(coordinate_system)}};
}

DirectionCoordinate Describe(const nlohmann::json& coordinate_system,
                             const std::vector<double>& l = kSamples,
                             const std::vector<double>& m = kSamples) {
    std::vector<Diagnostic> diagnostics;
    auto direction = DescribeDirection(Root(coordinate_system), l, m, diagnostics);
    Require(direction.has_value(), "a direction coordinate was expected");
    return *direction;
}

// The matrix that was never read back. Every fixture writes the identity, so a transposed copy is
// the same array; only an asymmetric one can tell.
void TestTheTransformationMatrixKeepsItsShape() {
    const auto direction = Describe({{"pixel_coordinate_transformation_matrix",
                                      nlohmann::json{{0.36, 0.48}, {-0.8, 0.6}}}});
    Require(Near(direction.transformation_matrix.at(0).at(0), 0.36), "element (0, 0)");
    Require(Near(direction.transformation_matrix.at(0).at(1), 0.48),
            "element (0, 1) -- a transposed copy puts -0.8 here");
    Require(Near(direction.transformation_matrix.at(1).at(0), -0.8),
            "element (1, 0) -- a transposed copy puts 0.48 here");
    Require(Near(direction.transformation_matrix.at(1).at(1), 0.6), "element (1, 1)");
}

// Absent, or the wrong shape, leaves the identity rather than a half-filled matrix.
void TestAMalformedMatrixLeavesTheIdentity() {
    for (const auto& matrix : {nlohmann::json{{1.0, 2.0, 3.0}}, nlohmann::json{"not a matrix"},
                               nlohmann::json{{"a", "b"}, {"c", "d"}}}) {
        const auto direction = Describe({{"pixel_coordinate_transformation_matrix", matrix}});
        Require(Near(direction.transformation_matrix.at(0).at(0), 1.0) &&
                    Near(direction.transformation_matrix.at(0).at(1), 0.0) &&
                    Near(direction.transformation_matrix.at(1).at(0), 0.0) &&
                    Near(direction.transformation_matrix.at(1).at(1), 1.0),
                "a matrix that is not two by two of numbers should leave the identity alone");
    }
}

// XRADIO writes radians; the descriptor reports degrees. Three separate places convert, and all
// three used to be unreachable.
void TestRadiansBecomeDegrees() {
    const double ra = 1.0;
    const double dec = 0.5;
    const auto direction = Describe({{"reference_direction", {{"data", {ra, dec}}}},
                                     {"native_pole_direction", {{"data", {0.0, M_PI / 2.0}}}}});
    Require(Near(direction.reference_value.at(0), ra * kRadToDeg), "the reference right ascension");
    Require(Near(direction.reference_value.at(1), dec * kRadToDeg), "the reference declination");
    Require(Near(direction.native_pole_direction.at(1), 90.0),
            "a native pole of pi over two is ninety degrees, which is what casacore calls latPole");
    Require(Near(direction.increment.at(0), 1.0e-4 * kRadToDeg), "the l increment");
    Require(Near(direction.increment.at(1), 1.0e-4 * kRadToDeg), "the m increment");
    Require(Near(direction.reference_pixel.at(0), 3.0), "the tangent point is the third sample, 1-based");
}

// An equinox arrives as a number or as a string with a leading letter naming the system.
void TestEquinoxIsReadInEveryFormXradioWrites() {
    const auto with = [](const nlohmann::json& equinox) {
        return Describe({{"reference_direction", {{"data", {0.0, 0.0}}, {"attrs", {{"frame", "fk5"}, {"equinox", equinox}}}}}});
    };
    Require(with(2000.0).equinox == 2000.0, "a numeric equinox");
    Require(with("J2000").equinox == 2000.0, "J2000");
    Require(with("B1950").equinox == 1950.0, "B1950");
    Require(with("j2000").equinox == 2000.0, "a lowercase prefix");
    Require(with("1950").equinox == 1950.0, "a bare year");
    Require(!with("not a year").equinox.has_value(), "something unparseable leaves it absent rather than zero");
}

void TestTheFrameAndProjectionAreUpperCased() {
    const auto direction =
        Describe({{"projection", "sin"},
                  {"reference_direction", {{"data", {0.0, 0.0}}, {"attrs", {{"frame", "icrs"}}}}}});
    Require(direction.projection == "SIN", "the projection is reported upper case");
    Require(direction.reference_frame == "ICRS", "and so is the frame");
}

void TestProjectionParametersAreCarried() {
    const auto direction = Describe({{"projection", "SIN"}, {"projection_parameters", {0.25, -0.5}}});
    Require(direction.projection_parameters.size() == 2, "both parameters");
    Require(Near(direction.projection_parameters.at(0), 0.25) &&
                Near(direction.projection_parameters.at(1), -0.5),
            "in the order they were written -- they become longPole and latPole");
}

// An unevenly sampled direction axis still reports an increment, and says so.
void TestAnUnevenAxisReportsAnIncrementAndADiagnostic() {
    const std::vector<double> uneven{0.0, 1.0e-4, 2.5e-4, 3.0e-4};
    std::vector<Diagnostic> diagnostics;
    const auto direction = DescribeDirection(Root({{"projection", "SIN"}}), uneven, kSamples, diagnostics);
    Require(direction.has_value(), "an uneven axis still describes a direction");
    Require(Near(direction->increment.at(0), 1.0e-4 * kRadToDeg), "with an increment in degrees");
    bool said_so = false;
    for (const auto& diagnostic : diagnostics) {
        said_so = said_so || diagnostic.code == carta::zarr::DiagnosticCode::nonuniform_axis;
    }
    Require(said_so, "and a diagnostic saying the samples were not evenly spaced");
}

// Nothing to describe at all.
void TestNothingToDescribe() {
    std::vector<Diagnostic> diagnostics;
    const nlohmann::json bare{{"type", "image_dataset"}};
    Require(!DescribeDirection(bare, {}, {}, diagnostics).has_value(),
            "no coordinate system and no samples describes no direction");
    Require(DescribeDirection(bare, kSamples, kSamples, diagnostics).has_value(),
            "samples alone still describe the pixel grid, even with no coordinate system");
}

}  // namespace

int main() {
    try {
        TestTheTransformationMatrixKeepsItsShape();
        TestAMalformedMatrixLeavesTheIdentity();
        TestRadiansBecomeDegrees();
        TestEquinoxIsReadInEveryFormXradioWrites();
        TestTheFrameAndProjectionAreUpperCased();
        TestProjectionParametersAreCarried();
        TestAnUnevenAxisReportsAnIncrementAndADiagnostic();
        TestNothingToDescribe();
        std::cout << "carta-zarr direction tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr direction tests failed: " << error.what() << '\n';
        return 1;
    }
}
