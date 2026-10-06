/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The beams a beam fit parameter table holds.
//
// The table is a three- or four-dimensional array read through a flat buffer, and what decides
// which element is which beam is dimension names against parameter labels. DescribeBeams takes the
// arrays already read, so a case is a table and its labels, with no store to write.
//
// Every plane and parameter gets a value that identifies it -- major = 100*frequency +
// 10*polarization, minor and the position angle one and two above it -- so a pair of transposed
// strides disagrees with the oracle rather than producing the same numbers in a different order.

#include "schema/xradio/beam_table.h"

#include "support/check.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using carta::zarr::Beam;
using carta::zarr::ErrorCode;
using carta::zarr::internal::xradio::DescribeBeams;
using carta::zarr::internal::zarr::ArrayMetadata;
using carta::zarr::internal::zarr::NumericArray;

using carta::zarr::testing::Require;

// The value a plane's parameters are built from, so that every beam identifies where it came from.
double Base(std::uint64_t time, std::uint64_t frequency, std::uint64_t polarization) {
    return (1000.0 * static_cast<double>(time)) + (100.0 * static_cast<double>(frequency)) +
           (10.0 * static_cast<double>(polarization));
}

// C order, parameters fastest, laid out as the labels are: minor, major, pa.
std::vector<double> Values(std::uint64_t times, std::uint64_t frequencies, std::uint64_t polarizations) {
    std::vector<double> values;
    for (std::uint64_t t = 0; t < times; ++t) {
        for (std::uint64_t f = 0; f < frequencies; ++f) {
            for (std::uint64_t p = 0; p < polarizations; ++p) {
                const double base = Base(t, f, p);
                values.push_back(base);
                values.push_back(base + 1.0);
                values.push_back(base + 2.0);
            }
        }
    }
    return values;
}

// A table as Store::ReadNumericArray hands it back: its values bound to its own document.
NumericArray Table(std::vector<std::uint64_t> shape, std::vector<std::string> dimensions, std::vector<double> values,
                   const std::string& unit = "rad") {
    ArrayMetadata metadata;
    metadata.shape = std::move(shape);
    metadata.dimension_names = std::move(dimensions);
    metadata.data_type = "float64";
    metadata.attributes = nlohmann::json{{"units", unit}};
    auto table = NumericArray::Make("BEAM", std::move(metadata), std::move(values));
    Require(static_cast<bool>(table), "the table's values are not as many as its shape declares");
    return std::move(table.value());
}

// Deliberately not in major/minor/pa order: a parameter is located by its label, not by position.
const std::vector<std::string> kLabels{"minor", "major", "pa"};

void RequireBeamMatchesItsPlane(const Beam& beam, const std::string& what) {
    const double base = Base(static_cast<std::uint64_t>(beam.time), static_cast<std::uint64_t>(beam.channel),
                             static_cast<std::uint64_t>(beam.polarization));
    Require(beam.major == base + 1.0, "beam major was read from the wrong element " + what);
    Require(beam.minor == base, "beam minor was read from the wrong element " + what);
    Require(beam.position_angle == base + 2.0, "beam position angle was read from the wrong element " + what);
}

void TestParametersAreLocatedByLabel() {
    const auto beams = DescribeBeams(
        Table({1, 3, 2, 3}, {"time", "frequency", "polarization", "beam_params_label"}, Values(1, 3, 2)), kLabels);
    Require(static_cast<bool>(beams), "a well-formed beam table was not read");
    Require(beams.value().size() == 6, "the beam table did not decode one beam per frequency and polarization");
    for (const auto& beam : beams.value()) {
        RequireBeamMatchesItsPlane(beam, "for a single-plane table");
        Require(beam.unit == "rad", "beam unit was not read from the beam array attributes");
    }
}

// A table with more than one time plane is read whole, not as its first plane only. Time varies
// slowest, so a single-plane table reads back in the same order either way.
void TestEveryTimePlaneIsReported() {
    const auto beams = DescribeBeams(
        Table({2, 3, 2, 3}, {"time", "frequency", "polarization", "beam_params_label"}, Values(2, 3, 2)), kLabels);
    Require(static_cast<bool>(beams), "a multi-plane beam table was not read");
    Require(beams.value().size() == 12, "the multi-plane beam table did not report every plane");
    for (const auto& beam : beams.value()) {
        RequireBeamMatchesItsPlane(beam, "for a multi-plane table");
    }
    const auto& first = beams.value().front();
    Require(first.time == 0 && first.channel == 0 && first.polarization == 0,
            "time did not vary slowest, so single-plane callers would see a different order");
}

// A table need not carry a time dimension; an absent one is a single implicit plane. Addressing the
// array must not insist on naming a dimension the array lacks.
void TestAnAbsentTimeDimensionIsOnePlane() {
    const auto beams =
        DescribeBeams(Table({3, 2, 3}, {"frequency", "polarization", "beam_params_label"}, Values(1, 3, 2)), kLabels);
    Require(static_cast<bool>(beams), "a beam table without a time dimension was not readable");
    Require(beams.value().size() == 6, "the time-less beam table did not report one beam per plane");
    for (const auto& beam : beams.value()) {
        Require(beam.time == 0, "a beam table without a time dimension reported a nonzero time index");
        RequireBeamMatchesItsPlane(beam, "without a time dimension");
    }
}

// The image named a beam table, so a table that cannot be addressed is an error rather than an
// empty list -- an empty list is the answer for an image with no beam at all, which this is not.
void TestATableMissingADimensionIsAnError() {
    const auto beams = DescribeBeams(Table({3, 3}, {"frequency", "beam_params_label"}, Values(1, 3, 1)), kLabels);
    Require(!beams && beams.error().code == ErrorCode::invalid_metadata,
            "a beam table with no polarization dimension was read anyway");
    Require(beams.error().node_path == "BEAM", "the error did not name the beam table");
}

// Same reasoning for the labels: parameters that cannot be located are a failure to read a beam,
// not an image without one.
void TestLabelsMustNameEveryParameter() {
    const auto beams =
        DescribeBeams(Table({1, 3, 2, 3}, {"time", "frequency", "polarization", "beam_params_label"}, Values(1, 3, 2)),
                      {"minor", "major", "angle"});
    Require(!beams && beams.error().code == ErrorCode::invalid_metadata,
            "labels that do not name a position angle produced beams anyway");
    Require(beams.error().node_path == "beam_params_label", "the error did not name the label array");
}

// The labels are the coordinate of the table's parameter dimension, one label a parameter. More
// labels than parameters, or fewer, are another table's, and refused as a malformed store rather
// than read as far as they happen to reach.
void TestLabelsAreOneAParameter() {
    const auto more =
        DescribeBeams(Table({1, 3, 2, 3}, {"time", "frequency", "polarization", "beam_params_label"}, Values(1, 3, 2)),
                      {"minor", "major", "pa", "extra"});
    Require(!more && more.error().code == ErrorCode::invalid_metadata,
            "more labels than the table has parameters produced beams anyway");
    Require(more.error().node_path == "beam_params_label", "the error did not name the label array");

    const std::vector<double> four_parameters(static_cast<std::size_t>(1 * 3 * 2 * 4), 1.0);
    const auto fewer = DescribeBeams(
        Table({1, 3, 2, 4}, {"time", "frequency", "polarization", "beam_params_label"}, four_parameters), kLabels);
    Require(!fewer && fewer.error().code == ErrorCode::invalid_metadata,
            "fewer labels than the table has parameters produced beams anyway");
    Require(fewer.error().node_path == "beam_params_label", "the error did not name the label array");
}

}  // namespace

int main() {
    try {
        TestParametersAreLocatedByLabel();
        TestEveryTimePlaneIsReported();
        TestAnAbsentTimeDimensionIsOnePlane();
        TestATableMissingADimensionIsAnError();
        TestLabelsMustNameEveryParameter();
        TestLabelsAreOneAParameter();
        std::cout << "carta-zarr beam table tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr beam table tests failed: " << error.what() << '\n';
        return 1;
    }
}
