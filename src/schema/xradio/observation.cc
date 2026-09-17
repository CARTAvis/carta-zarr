/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "observation.h"

#include "attributes.h"

#include <array>
#include <cmath>
#include <string>

namespace carta::zarr::internal::xradio {

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
        if (const auto* const attributes = MemberObject(*obsdate, "attrs")) {
            observation.timesys = Upper(AttributeString(*attributes, "scale"));
        }
        if (const auto* const data = Member(*obsdate, "data"); data != nullptr) {
            if (data->is_string()) {
                observation.date_obs = data->get<std::string>();
            } else if (data->is_number()) {
                observation.mjd_obs = data->get<double>();
            }
        }
    }
    return observation;
}

}  // namespace carta::zarr::internal::xradio
