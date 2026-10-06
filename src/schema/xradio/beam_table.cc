/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "beam_table.h"

#include "attributes.h"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

struct BeamParameterIndices {
    std::optional<std::size_t> major;
    std::optional<std::size_t> minor;
    std::optional<std::size_t> position_angle;
};

BeamParameterIndices FindBeamParameterIndices(const std::vector<std::string>& labels) {
    BeamParameterIndices indices;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        const std::string label = Upper(labels.at(index));
        if (label == "MAJOR") {
            indices.major = index;
        } else if (label == "MINOR") {
            indices.minor = index;
        } else if (label == "PA") {
            indices.position_angle = index;
        }
    }
    return indices;
}

Result<double> ReadBeamValue(const zarr::NumericArray& table, std::uint64_t channel, std::uint64_t polarization,
                             std::uint64_t parameter, bool has_time_dimension, std::uint64_t time) {
    std::vector<zarr::NumericArray::NamedIndex> indices{
        {"frequency", channel}, {"polarization", polarization}, {"beam_params_label", parameter}};
    if (has_time_dimension) {
        indices.emplace_back("time", time);
    }
    return table.At(indices);
}

}  // namespace

Result<std::vector<Beam>> DescribeBeams(const zarr::NumericArray& table,
                                        const std::vector<std::string>& parameter_labels) {
    const auto& beam_metadata = table.metadata();
    const std::string beam_unit = AttributeString(beam_metadata.attributes, "units");
    const auto parameter_indices = FindBeamParameterIndices(parameter_labels);

    const auto freq_dim = zarr::FindDimensionIndex(beam_metadata, "frequency");
    const auto pol_dim = zarr::FindDimensionIndex(beam_metadata, "polarization");
    const auto param_dim = zarr::FindDimensionIndex(beam_metadata, "beam_params_label");
    const auto time_dim = zarr::FindDimensionIndex(beam_metadata, "time");

    if (!freq_dim || !pol_dim || !param_dim) {
        return Error{ErrorCode::invalid_metadata,
                     "Beam table does not carry the frequency, polarization and parameter dimensions", table.node()};
    }
    // The labels are the coordinate of the parameter dimension, so there is one for each parameter.
    // A count of either other than the other's is labels written for another table.
    if (parameter_labels.size() != beam_metadata.shape.at(*param_dim)) {
        return Error{ErrorCode::invalid_metadata,
                     "Beam parameter labels are " + std::to_string(parameter_labels.size()) +
                         " where the beam table has " + std::to_string(beam_metadata.shape.at(*param_dim)) +
                         " parameters",
                     "beam_params_label"};
    }
    if (!parameter_indices.major || !parameter_indices.minor || !parameter_indices.position_angle) {
        return Error{ErrorCode::invalid_metadata,
                     "Beam parameter labels do not name a major axis, a minor axis and a position angle",
                     "beam_params_label"};
    }

    const std::uint64_t n_time = time_dim ? beam_metadata.shape.at(*time_dim) : 1;
    const std::uint64_t n_chan = beam_metadata.shape.at(*freq_dim);
    const std::uint64_t n_pol = beam_metadata.shape.at(*pol_dim);

    // Time varies slowest so that a single-plane beam table reads back in the order it always has.
    std::vector<Beam> beams;
    beams.reserve(static_cast<std::size_t>(n_time * n_chan * n_pol));
    for (std::uint64_t t = 0; t < n_time; ++t) {
        for (std::uint64_t c = 0; c < n_chan; ++c) {
            for (std::uint64_t p = 0; p < n_pol; ++p) {
                Beam beam;
                beam.time = static_cast<std::size_t>(t);
                beam.channel = static_cast<std::size_t>(c);
                beam.polarization = static_cast<std::size_t>(p);
                beam.unit = beam_unit;

                for (const auto& [parameter, field] :
                     {std::pair{*parameter_indices.major, &beam.major},
                      std::pair{*parameter_indices.minor, &beam.minor},
                      std::pair{*parameter_indices.position_angle, &beam.position_angle}}) {
                    auto value = ReadBeamValue(table, c, p, parameter, time_dim.has_value(), t);
                    if (!value) {
                        return value.error();
                    }
                    *field = value.value();
                }
                beams.push_back(std::move(beam));
            }
        }
    }

    return beams;
}

}  // namespace carta::zarr::internal::xradio
