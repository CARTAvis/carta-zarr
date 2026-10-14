/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "observation.h"

#include "attributes.h"

#include <array>
#include <cmath>
#include <optional>
#include <string>

namespace carta::zarr::internal::xradio {

namespace {

constexpr double kSecondsPerDay = 86400.0;
// MJD of 1970-01-01, where unix time starts.
constexpr double kUnixEpochMjd = 40587.0;

// The MJD a number names, by the format written beside it: XRADIO writes "mjd", and a number with no
// format was read as one for XRADIO 1.2, whose stores did not always write it. Seconds since 1970 are
// what XRADIO's time coordinate holds.
std::optional<double> MjdOfNumber(double value, const std::string& format, const SchemaRow& row) {
    if (format == "MJD" || (format.empty() && row.unformatted_date_is_mjd)) {
        return value;
    }
    if (format == "UNIX") {
        return kUnixEpochMjd + (value / kSecondsPerDay);
    }
    return std::nullopt;
}

}  // namespace

ObservationInfo DescribeObservation(const zarr::ArrayMetadata& image, const SchemaRow& row) {
    ObservationInfo observation;
    observation.object_name = AttributeString(image.attributes, "object_name");
    observation.observer = AttributeString(image.attributes, "observer");
    if (const auto* const telescope = MemberObject(image.attributes, "telescope")) {
        observation.telescope_name = AttributeString(*telescope, "name");

        const auto* const direction = MemberObject(*telescope, "direction");
        const auto* const distance = MemberObject(*telescope, "distance");
        const auto* const direction_data = direction == nullptr ? nullptr : MemberArray(*direction, "data");
        const auto* const distance_data = distance == nullptr ? nullptr : MemberArray(*distance, "data");
        // Every element is checked before it is converted, not just the array around it: a value of
        // the wrong type throws out of nlohmann, and this is optional metadata reached while
        // describing an image, where the caller is holding a Result and expecting a diagnostic at
        // worst rather than an exception.
        if (direction_data != nullptr && direction_data->size() >= 2 && direction_data->at(0).is_number() &&
            direction_data->at(1).is_number() && distance_data != nullptr && !distance_data->empty() &&
            distance_data->at(0).is_number()) {
            const double lon = direction_data->at(0).get<double>();
            const double lat = direction_data->at(1).get<double>();
            const double radius = distance_data->at(0).get<double>();
            observation.observatory_position = std::array<double, 3>{
                radius * std::cos(lat) * std::cos(lon), radius * std::cos(lat) * std::sin(lon), radius * std::sin(lat)};
        }
    }
    if (const auto* const obsdate = MemberObject(image.attributes, "obsdate")) {
        std::string format;
        if (const auto* const attributes = MemberObject(*obsdate, "attrs")) {
            observation.timesys = Upper(AttributeString(*attributes, "scale"));
            format = Upper(AttributeString(*attributes, "format"));
        }
        // A number, as every XRADIO writes it. A string date is no XRADIO's, and is not read.
        if (const auto value = MemberNumber(*obsdate, "data")) {
            observation.mjd_obs = MjdOfNumber(*value, format, row);
        }
    }
    return observation;
}

}  // namespace carta::zarr::internal::xradio
