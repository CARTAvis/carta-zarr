/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "beam_table.h"

#include "../../zarr/array_view.h"
#include "attributes.h"

#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

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

Result<double> ReadBeamValue(const zarr::ArrayView& values, std::uint64_t channel, std::uint64_t polarization,
                             std::uint64_t parameter, bool has_time_dimension, std::uint64_t time) {
    std::vector<zarr::ArrayView::NamedIndex> indices{
        {"frequency", channel}, {"polarization", polarization}, {"beam_params_label", parameter}};
    if (has_time_dimension) {
        indices.emplace_back("time", time);
    }
    return values.At(indices);
}

}  // namespace

Result<std::vector<Beam>> DescribeBeams(const zarr::ArrayMetadata& beam_metadata, std::string_view beam_node,
                                        const std::vector<std::string>& parameter_labels,
                                        const std::vector<double>& values) {
    const std::string beam_unit = AttributeString(beam_metadata.attributes, "units");
    const auto parameter_indices = FindBeamParameterIndices(parameter_labels);

    const auto freq_dim = zarr::FindDimensionIndex(beam_metadata, "frequency");
    const auto pol_dim = zarr::FindDimensionIndex(beam_metadata, "polarization");
    const auto param_dim = zarr::FindDimensionIndex(beam_metadata, "beam_params_label");
    const auto time_dim = zarr::FindDimensionIndex(beam_metadata, "time");

    if (!freq_dim || !pol_dim || !param_dim) {
        return MakeError(ErrorCode::invalid_metadata,
                         "Beam table does not carry the frequency, polarization and parameter dimensions",
                         std::string(beam_node));
    }
    if (!parameter_indices.major || !parameter_indices.minor || !parameter_indices.position_angle) {
        return MakeError(ErrorCode::invalid_metadata,
                         "Beam parameter labels do not name a major axis, a minor axis and a position angle",
                         "beam_params_label");
    }

    const std::uint64_t n_time = time_dim ? beam_metadata.shape.at(*time_dim) : 1;
    const std::uint64_t n_chan = beam_metadata.shape.at(*freq_dim);
    const std::uint64_t n_pol = beam_metadata.shape.at(*pol_dim);

    // Time varies slowest so that a single-plane beam table reads back in the order it always has.
    const zarr::ArrayView beam_values(beam_metadata, values);
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
                    auto value = ReadBeamValue(beam_values, c, p, parameter, time_dim.has_value(), t);
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
