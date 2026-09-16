/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "direction.h"

#include "attributes.h"
#include "linear_axis.h"

#include "../../zarr/array_metadata.h"

#include <iterator>
#include <string>
#include <utility>

namespace carta::zarr::internal::xradio {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

std::optional<DirectionCoordinate> DescribeDirection(const nlohmann::json& root_attributes,
                                                     const std::vector<double>& l_values,
                                                     const std::vector<double>& m_values,
                                                     std::vector<Diagnostic>& diagnostics) {
    const auto* coordinate_system = ObjectMember(root_attributes, "coordinate_system_info");
    const bool has_coordinate_system = coordinate_system != nullptr && coordinate_system->is_object();
    if (!has_coordinate_system && (l_values.empty() || m_values.empty())) {
        return std::nullopt;
    }

    DirectionCoordinate direction;
    if (has_coordinate_system) {
        const auto& cs_info = *coordinate_system;
        if (const auto* projection = ObjectMember(cs_info, "projection");
            projection != nullptr && projection->is_string()) {
            direction.projection = Upper(projection->get<std::string>());
        }
        if (const auto* reference_direction = ObjectMember(cs_info, "reference_direction");
            reference_direction != nullptr && reference_direction->is_object()) {
            if (const auto* data = ObjectMember(*reference_direction, "data");
                data != nullptr && zarr_metadata::IsNumericVector(*data, 2)) {
                direction.reference_value.at(0) = data->at(0).get<double>() * kRadToDeg;
                direction.reference_value.at(1) = data->at(1).get<double>() * kRadToDeg;
            }
            if (const auto* attributes = ObjectMember(*reference_direction, "attrs");
                attributes != nullptr && attributes->is_object()) {
                direction.reference_frame = Upper(AttributeString(*attributes, "frame"));
                if (const auto* equinox = ObjectMember(*attributes, "equinox"); equinox != nullptr) {
                    if (equinox->is_number()) {
                        direction.equinox = equinox->get<double>();
                    } else if (equinox->is_string()) {
                        const std::string value = equinox->get<std::string>();
                        const std::size_t position = (value.size() > 1 && (value.at(0) == 'J' || value.at(0) == 'B' ||
                                                                           value.at(0) == 'j' || value.at(0) == 'b'))
                                                       ? 1
                                                       : 0;
                        try {
                            direction.equinox = std::stod(value.substr(position));
                        } catch (...) {
                        }
                    }
                }
            }
        }
        if (const auto* parameters = ObjectMember(cs_info, "projection_parameters");
            parameters != nullptr && parameters->is_array()) {
            for (const auto& value : *parameters) {
                if (value.is_number()) {
                    direction.projection_parameters.push_back(value.get<double>());
                }
            }
        }
        if (const auto* native_pole = ObjectMember(cs_info, "native_pole_direction");
            native_pole != nullptr && native_pole->is_object()) {
            const auto* data = ObjectMember(*native_pole, "data");
            if (data != nullptr && zarr_metadata::IsNumericVector(*data, 2)) {
                direction.native_pole_direction.at(0) = data->at(0).get<double>() * kRadToDeg;
                direction.native_pole_direction.at(1) = data->at(1).get<double>() * kRadToDeg;
            }
        }
        if (const auto* matrix = ObjectMember(cs_info, "pixel_coordinate_transformation_matrix");
            matrix != nullptr && zarr_metadata::IsNumericMatrix(*matrix, 2, 2)) {
            direction.transformation_matrix.at(0).at(0) = matrix->at(0).at(0).get<double>();
            direction.transformation_matrix.at(0).at(1) = matrix->at(0).at(1).get<double>();
            direction.transformation_matrix.at(1).at(0) = matrix->at(1).at(0).get<double>();
            direction.transformation_matrix.at(1).at(1) = matrix->at(1).at(1).get<double>();
        }
    }

    const auto set_direction_axis = [&](const std::vector<double>& values, double& increment, double& reference_pixel,
                                        std::string_view name) {
        auto fit = FitDirectionAxis(values, name);
        if (fit.increment) {
            increment = *fit.increment;
        }
        if (fit.reference_pixel) {
            reference_pixel = *fit.reference_pixel;
        }
        diagnostics.insert(diagnostics.end(), std::make_move_iterator(fit.diagnostics.begin()),
                           std::make_move_iterator(fit.diagnostics.end()));
    };
    set_direction_axis(l_values, direction.increment.at(0), direction.reference_pixel.at(0), "l");
    set_direction_axis(m_values, direction.increment.at(1), direction.reference_pixel.at(1), "m");
    return direction;
}

}  // namespace carta::zarr::internal::xradio
