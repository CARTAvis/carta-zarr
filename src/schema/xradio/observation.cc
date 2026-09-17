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
    if (const auto* object_name = ObjectMember(image.attributes, "object_name");
        object_name != nullptr && object_name->is_string()) {
        observation.object_name = object_name->get<std::string>();
    }
    if (const auto* observer = ObjectMember(image.attributes, "observer");
        observer != nullptr && observer->is_string()) {
        observation.observer = observer->get<std::string>();
    }
    if (const auto* telescope = ObjectMember(image.attributes, "telescope");
        telescope != nullptr && telescope->is_object()) {
        if (const auto* name = ObjectMember(*telescope, "name"); name != nullptr && name->is_string()) {
            observation.telescope_name = name->get<std::string>();
        }
        const auto* direction = ObjectMember(*telescope, "direction");
        const auto* distance = ObjectMember(*telescope, "distance");
        const nlohmann::json* direction_data = nullptr;
        if (direction != nullptr) {
            direction_data = ObjectMember(*direction, "data");
        }
        const nlohmann::json* distance_data = nullptr;
        if (distance != nullptr) {
            distance_data = ObjectMember(*distance, "data");
        }
        // Every element is checked before it is converted, not just the array around it: a value of
        // the wrong type throws out of nlohmann, and this is optional metadata reached while
        // describing an image, where the caller is holding a Result and expecting a diagnostic at
        // worst rather than an exception.
        if (direction != nullptr && direction->is_object() && direction_data != nullptr && direction_data->is_array() &&
            direction_data->size() >= 2 && direction_data->at(0).is_number() && direction_data->at(1).is_number() &&
            distance != nullptr && distance->is_object() && distance_data != nullptr && distance_data->is_array() &&
            !distance_data->empty() && distance_data->at(0).is_number()) {
            const double lon = direction_data->at(0).get<double>();
            const double lat = direction_data->at(1).get<double>();
            const double radius = distance_data->at(0).get<double>();
            observation.observatory_position = std::array<double, 3>{
                radius * std::cos(lat) * std::cos(lon), radius * std::cos(lat) * std::sin(lon), radius * std::sin(lat)};
        }
    }
    if (const auto* obsdate = ObjectMember(image.attributes, "obsdate"); obsdate != nullptr && obsdate->is_object()) {
        if (const auto* attributes = ObjectMember(*obsdate, "attrs");
            attributes != nullptr && attributes->is_object()) {
            observation.timesys = Upper(AttributeString(*attributes, "scale"));
        }
        if (const auto* data = ObjectMember(*obsdate, "data"); data != nullptr) {
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
