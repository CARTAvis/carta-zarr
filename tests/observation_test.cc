/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// What an image variable's attributes say about the observation that produced it.
//
// Optional metadata comes from a file, so its shape is whatever was written rather than whatever
// the schema describes. The contract is that a field it cannot read is skipped and the rest is kept
// -- nothing here can close an image. This takes JSON, so no store is written.
//
// The geodetic-to-cartesian conversion is asserted here for the first time. It was reachable only
// through a store, and no fixture wrote a position that would show a swapped axis.

#include "schema/xradio/observation.h"

#include "support/check.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <optional>
#include <string>

namespace {

using carta::zarr::ObservationInfo;
using carta::zarr::internal::xradio::DescribeObservation;
using ArrayMetadata = carta::zarr::internal::zarr::ArrayMetadata;

using carta::zarr::testing::Require;

bool Near(double left, double right) {
    return std::abs(left - right) <= 1.0e-6 * std::max({1.0, std::abs(left), std::abs(right)});
}

ObservationInfo Describe(nlohmann::json attributes) {
    ArrayMetadata image;
    image.data_type = "float32";
    image.attributes = std::move(attributes);
    return DescribeObservation(image);
}

void TestTheReadableFieldsAreReported() {
    const auto observation = Describe(nlohmann::json{
        {"object_name", "Zarr test source"}, {"observer", "A. Observer"}, {"telescope", {{"name", "Test scope"}}}});
    Require(observation.object_name == "Zarr test source", "object_name was not read");
    Require(observation.observer == "A. Observer", "observer was not read");
    Require(observation.telescope_name == "Test scope", "the telescope name was not read");
}

// longitude, latitude and radius become OBSGEO-X, Y, Z. Nothing in the repository wrote a position
// where swapping two of the three would change the answer, so this one is deliberately asymmetric.
void TestTheTelescopePositionBecomesCartesian() {
    const double lon = 0.5;
    const double lat = 0.25;
    const double radius = 6371000.0;
    const auto observation = Describe(nlohmann::json{
        {"telescope",
         {{"name", "Test scope"}, {"direction", {{"data", {lon, lat}}}}, {"distance", {{"data", {radius}}}}}}});
    Require(observation.observatory_position.has_value(), "a well-formed telescope position was not converted");
    const auto& position = *observation.observatory_position;
    Require(Near(position.at(0), radius * std::cos(lat) * std::cos(lon)), "OBSGEO-X was not the expected element");
    Require(Near(position.at(1), radius * std::cos(lat) * std::sin(lon)), "OBSGEO-Y was not the expected element");
    Require(Near(position.at(2), radius * std::sin(lat)), "OBSGEO-Z was not the expected element");
}

// A telescope position holding a string where a number belongs does not throw out of nlohmann past
// the Result the caller is holding. It is a value the image can do without: the readable metadata
// beside it survives and only the position is missing.
void TestAMalformedPositionIsSkippedRatherThanConverted() {
    const auto observation = Describe(nlohmann::json{
        {"object_name", "Zarr test source"},
        {"telescope",
         {{"name", "Test scope"}, {"direction", {{"data", {"north", 0.0}}}}, {"distance", {{"data", {6371000.0}}}}}}});
    Require(observation.telescope_name == "Test scope", "the readable telescope metadata was dropped as well");
    Require(observation.object_name == "Zarr test source", "an unrelated attribute was dropped as well");
    Require(!observation.observatory_position.has_value(),
            "a telescope position holding a string was converted rather than skipped");
}

// Half a position is not a position: a distance with no direction, or a direction of one element,
// leaves it empty rather than filling what it can.
void TestAnIncompletePositionIsNotHalfConverted() {
    const auto no_direction =
        Describe(nlohmann::json{{"telescope", {{"name", "T"}, {"distance", {{"data", {6371000.0}}}}}}});
    Require(!no_direction.observatory_position.has_value(), "a distance with no direction produced a position");

    const auto short_direction = Describe(nlohmann::json{
        {"telescope", {{"name", "T"}, {"direction", {{"data", {0.5}}}}, {"distance", {{"data", {6371000.0}}}}}}});
    Require(!short_direction.observatory_position.has_value(), "a one-element direction produced a position");
}

// obsdate carries the scale beside the value. The scale is reported upper-cased, because that is what
// a FITS TIMESYS is, and the value as the MJD it names in that scale.
void TestObsdateSplitsIntoScaleAndValue() {
    const auto numeric =
        Describe(nlohmann::json{{"obsdate", {{"attrs", {{"scale", "tai"}, {"format", "mjd"}}}, {"data", 61000.5}}}});
    Require(numeric.timesys == "TAI", "the obsdate scale was not upper-cased");
    Require(numeric.mjd_obs.has_value() && Near(*numeric.mjd_obs, 61000.5),
            "a numeric obsdate was not reported as an MJD");
}

std::optional<double> MjdOf(nlohmann::json obsdate) {
    return Describe(nlohmann::json{{"obsdate", std::move(obsdate)}}).mjd_obs;
}

// A number is read by the format written beside it: XRADIO writes "mjd", and seconds since 1970 under
// "unix" are converted, while a format this does not know yields no MJD rather than a wrong one. A
// number with no format is nothing, as schema 0.0.2 requires the format. A string is no XRADIO's date,
// and is not read.
void TestObsdateSpellings() {
    Require(MjdOf({{"attrs", {{"format", "MJD"}}}, {"data", 59000.5}}) == 59000.5,
            "an mjd-format number was not an MJD");
    const auto unix_seconds = MjdOf({{"attrs", {{"format", "unix"}}}, {"data", 1590926400.0}});
    Require(unix_seconds && Near(*unix_seconds, 59000.5), "a unix-format number was not converted to an MJD");
    Require(!MjdOf({{"attrs", {{"format", "jyear"}}}, {"data", 2020.4}}), "a number of an unknown format was an MJD");

    Require(!MjdOf({{"data", 59000.5}}), "a number with no format was read as an MJD");
    Require(!MjdOf({{"attrs", {{"format", "mjd"}}}, {"data", "2020-05-31T12:00:00"}}), "a string date was read");
}

// An image carrying none of it is described, not refused.
void TestNothingToDescribe() {
    const auto observation = Describe(nlohmann::json::object());
    Require(observation.object_name.empty() && observation.telescope_name.empty() && observation.timesys.empty(),
            "an image with no observation attributes reported some");
    Require(!observation.observatory_position.has_value() && !observation.mjd_obs.has_value(),
            "an image with no observation attributes reported a position or a date");
}

}  // namespace

int main() {
    try {
        TestTheReadableFieldsAreReported();
        TestTheTelescopePositionBecomesCartesian();
        TestAMalformedPositionIsSkippedRatherThanConverted();
        TestAnIncompletePositionIsNotHalfConverted();
        TestObsdateSplitsIntoScaleAndValue();
        TestObsdateSpellings();
        TestNothingToDescribe();
        std::cout << "carta-zarr observation tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr observation tests failed: " << error.what() << '\n';
        return 1;
    }
}
