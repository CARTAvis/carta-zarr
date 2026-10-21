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
#include <cstdint>
#include <exception>
#include <iostream>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace {

using carta::zarr::ObservationInfo;
using carta::zarr::internal::xradio::DescribeObservation;
using carta::zarr::internal::xradio::DescribeUserKeywords;
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
    Require(numeric.time_scale == "TAI", "the obsdate scale was not upper-cased");
    Require(numeric.observation_mjd.has_value() && Near(*numeric.observation_mjd, 61000.5),
            "a numeric obsdate was not reported as an MJD");
}

std::optional<double> MjdOf(nlohmann::json obsdate) {
    return Describe(nlohmann::json{{"obsdate", std::move(obsdate)}}).observation_mjd;
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
    Require(observation.object_name.empty() && observation.telescope_name.empty() && observation.time_scale.empty(),
            "an image with no observation attributes reported some");
    Require(!observation.observatory_position.has_value() && !observation.observation_mjd.has_value(),
            "an image with no observation attributes reported a position or a date");
    Require(!observation.pointing_center.has_value(), "an image with no pointing center reported one");
}

// Where the telescope pointed, in degrees and in the frame named as the reference direction's is,
// equinox included. XRADIO writes it in radians, from OBSRA and OBSDEC or from the reference
// direction; either way it is reported.
void TestThePointingCenterIsReportedInDegrees() {
    const double pi = std::acos(-1.0);
    const auto observation = Describe(nlohmann::json{
        {"pointing_center",
         {{"data", {pi / 4.0, pi / 6.0}},
          {"attrs", {{"frame", "fk5"}, {"type", "sky_coord"}, {"units", "rad"}, {"equinox", "j2000.0"}}}}}});
    Require(observation.pointing_center.has_value(), "the pointing center was not read");
    const auto& pointing = *observation.pointing_center;
    Require(Near(pointing.value.at(0), 45.0) && Near(pointing.value.at(1), 30.0),
            "the pointing center was not converted to degrees");
    Require(pointing.reference_frame == "FK5", "the pointing center's frame was not named as a reference frame is");
    Require(pointing.equinox && Near(*pointing.equinox, 2000.0), "the pointing center's equinox was not read");

    Require(!Describe(nlohmann::json{{"pointing_center", {{"data", {1.0}}}}}).pointing_center,
            "a pointing center of one angle was reported");
    Require(!Describe(nlohmann::json{{"pointing_center", "somewhere"}}).pointing_center,
            "a pointing center that is not a measure was reported");
}

// The dataset's keywords and the image's, the image's winning; each kind of value XRADIO can carry,
// and those it cannot report left out by name, without failing anything else.
void TestUserKeywordsAreMergedAndTyped() {
    ArrayMetadata image;
    image.data_type = "float32";
    image.attributes = nlohmann::json::parse(R"({"user": {
        "instrume": "ALMA", "distance": 101.5, "nchunk": 3, "offset": -2, "flagged": true,
        "shared": "image", "miscinfo": {"a": 1}, "history": [1, 2], "nothing": null,
        "huge": 18446744073709551615}})");
    const auto root = nlohmann::json::parse(R"({"user": {"shared": "dataset", "origin": "root"}})");
    std::vector<carta::zarr::Diagnostic> diagnostics;
    const auto keywords = DescribeUserKeywords(root, image, "SKY", diagnostics);

    const auto find = [&](const std::string& name) -> const carta::zarr::UserKeyword* {
        const auto found = std::find_if(keywords.begin(), keywords.end(),
                                        [&](const carta::zarr::UserKeyword& keyword) { return keyword.name == name; });
        return found == keywords.end() ? nullptr : &*found;
    };
    const auto* const instrume = find("instrume");
    Require(instrume && std::get<std::string>(instrume->value) == "ALMA", "a string keyword was not reported");
    const auto* const distance = find("distance");
    Require(distance && std::holds_alternative<double>(distance->value) && std::get<double>(distance->value) == 101.5,
            "a number keyword was not reported as a number");
    const auto* const nchunk = find("nchunk");
    Require(nchunk && std::holds_alternative<std::int64_t>(nchunk->value) && std::get<std::int64_t>(nchunk->value) == 3,
            "an integer keyword was not reported as an integer");
    const auto* const offset = find("offset");
    Require(offset && std::get<std::int64_t>(offset->value) == -2, "a negative integer keyword was not reported");
    const auto* const flagged = find("flagged");
    Require(flagged && std::get<bool>(flagged->value), "a boolean keyword was not reported as a boolean");
    const auto* const shared = find("shared");
    Require(shared && std::get<std::string>(shared->value) == "image", "the image's keyword did not win");
    const auto* const origin = find("origin");
    Require(origin && std::get<std::string>(origin->value) == "root", "the dataset's own keyword was not reported");
    Require(std::is_sorted(keywords.begin(), keywords.end(),
                           [](const auto& left, const auto& right) { return left.name < right.name; }),
            "the keywords were not sorted by name");

    for (const std::string left_out : {"miscinfo", "history", "nothing", "huge"}) {
        Require(find(left_out) == nullptr, "keyword '" + left_out + "' was reported though it cannot be");
        const auto said = std::any_of(diagnostics.begin(), diagnostics.end(), [&](const auto& diagnostic) {
            return diagnostic.code == carta::zarr::DiagnosticCode::unsupported_data_type &&
                   diagnostic.node_path == "SKY" && diagnostic.message.find("'" + left_out + "'") != std::string::npos;
        });
        Require(said, "keyword '" + left_out + "' was left out without a diagnostic naming it");
    }
    Require(diagnostics.size() == 4, "a keyword that was reported was also said to be left out");

    // One left out of the dataset's own is said to come from the dataset.
    std::vector<carta::zarr::Diagnostic> root_said;
    ArrayMetadata bare;
    bare.data_type = "float32";
    DescribeUserKeywords(nlohmann::json::parse(R"({"user": {"miscinfo": {}}})"), bare, "SKY", root_said);
    Require(root_said.size() == 1 && root_said.front().node_path == "/attributes/user",
            "a keyword left out of the dataset's own was not said to come from the dataset");
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
        TestThePointingCenterIsReportedInDegrees();
        TestUserKeywordsAreMergedAndTyped();
        std::cout << "carta-zarr observation tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr observation tests failed: " << error.what() << '\n';
        return 1;
    }
}
