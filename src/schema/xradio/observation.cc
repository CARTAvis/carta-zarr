/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "observation.h"

#include "attributes.h"
#include "direction.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>

namespace carta::zarr::internal::xradio {

namespace {

constexpr double kSecondsPerDay = 86400.0;
// MJD of 1970-01-01, where unix time starts.
constexpr double kUnixEpochMjd = 40587.0;

// The MJD a number names, by the format written beside it: XRADIO writes "mjd", and schema 0.0.2
// requires the format. Seconds since 1970 are what XRADIO's time coordinate holds.
std::optional<double> MjdOfNumber(double value, const std::string& format) {
    if (format == "MJD") {
        return value;
    }
    if (format == "UNIX") {
        return kUnixEpochMjd + (value / kSecondsPerDay);
    }
    return std::nullopt;
}

}  // namespace

ObservationInfo DescribeObservation(const zarr::ArrayMetadata& image) {
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
            observation.time_scale = Upper(AttributeString(*attributes, "scale"));
            format = Upper(AttributeString(*attributes, "format"));
        }
        // A number, as every XRADIO writes it. A string date is no XRADIO's, and is not read.
        if (const auto value = MemberNumber(*obsdate, "data")) {
            observation.observation_mjd = MjdOfNumber(*value, format);
        }
    }
    observation.pointing_center = ReadSkyDirection(MemberObject(image.attributes, "pointing_center"));
    return observation;
}

std::vector<UserKeyword> DescribeUserKeywords(const nlohmann::json& root_attributes, const zarr::ArrayMetadata& image,
                                              std::string_view image_id, std::vector<Diagnostic>& diagnostics) {
    // The dataset's first, then the image's over them, as XRADIO's writer merges the two. Where each
    // came from is kept for the diagnostic, which names the node the value was read from.
    nlohmann::json merged = nlohmann::json::object();
    std::map<std::string, std::string> source;
    const auto merge = [&](const nlohmann::json& attributes, const std::string& node) {
        if (const auto* const user = MemberObject(attributes, "user")) {
            for (const auto& [name, value] : user->items()) {
                merged[name] = value;
                source[name] = node;
            }
        }
    };
    merge(root_attributes, "/attributes/user");
    merge(image.attributes, std::string(image_id));

    std::vector<UserKeyword> keywords;
    for (const auto& [name, value] : merged.items()) {
        const auto left_out = [&](const std::string& why) {
            diagnostics.push_back(Diagnostic{DiagnosticCode::unsupported_data_type,
                                             "User keyword '" + name + "' " + why + ", and is not reported",
                                             source.at(name)});
        };
        if (value.is_boolean()) {
            keywords.push_back(UserKeyword{name, value.get<bool>()});
        } else if (value.is_number_unsigned()) {
            // A positive integer, which is how JSON is parsed whatever it was written as; one past what
            // a signed 64-bit integer holds is not given as a number it is not.
            const auto unsigned_value = value.get<std::uint64_t>();
            if (unsigned_value > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
                left_out("holds an integer too large for a signed 64-bit integer");
                continue;
            }
            keywords.push_back(UserKeyword{name, static_cast<std::int64_t>(unsigned_value)});
        } else if (value.is_number_integer()) {
            keywords.push_back(UserKeyword{name, value.get<std::int64_t>()});
        } else if (value.is_number_float()) {
            keywords.push_back(UserKeyword{name, value.get<double>()});
        } else if (value.is_string()) {
            keywords.push_back(UserKeyword{name, value.get<std::string>()});
        } else {
            left_out(std::string("holds ") + (value.is_object()  ? "a dictionary"
                                              : value.is_array() ? "a list"
                                                                 : "no value"));
        }
    }
    return keywords;
}

}  // namespace carta::zarr::internal::xradio
