/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "image.h"

#include "../../zarr/array_metadata.h"
#include "attributes.h"
#include "beam_table.h"
#include "direction.h"
#include "flag.h"
#include "linear_axis.h"
#include "observation.h"
#include "probe_report.h"

#include <carta-zarr/descriptor.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

constexpr std::string_view kVersion = "1.2";
constexpr std::array<std::string_view, 5> kSkyAxes{"time", "frequency", "polarization", "l", "m"};

int KnownImageRank(std::string_view image_id) {
    static constexpr std::array<std::string_view, 6> known{
        "SKY", "MODEL", "RESIDUAL", "POINT_SPREAD_FUNCTION", "PRIMARY_BEAM", "MASK_DECONVOLVE"};
    const auto* const found = std::find(known.begin(), known.end(), image_id);
    return found == known.end() ? static_cast<int>(known.size()) : static_cast<int>(found - known.begin());
}

// Every image carries a coordinate array for each axis it uses. Both XRADIO readers write all five
// unconditionally (xds_from_casacore.py and xds_from_fits.py both assign coords["time"]), so a
// missing coordinate means the store is malformed.
void RequirePresentCoordinates(ProbeReport& report, const zarr_metadata::ArrayMetadata& image) {
    for (const auto axis : kSkyAxes) {
        report.RequireCoordinateOf(image, axis,
                                   axis == "polarization" ? CoordinateKind::labels : CoordinateKind::numeric);
    }
}

void AppendDiagnostics(ImageDescriptor& descriptor, std::vector<Diagnostic> diagnostics) {
    descriptor.diagnostics.insert(descriptor.diagnostics.end(), std::make_move_iterator(diagnostics.begin()),
                                  std::make_move_iterator(diagnostics.end()));
}

constexpr std::array<std::string_view, 5> kApertureAxes{"time", "frequency", "polarization", "u", "v"};

// An image carries every axis of its plane. XRADIO writes optional coordinate arrays that share the
// image's spatial axes without being images at all: right_ascension and declination are float64 over
// (l, m) and carry no type attribute, so a rule keyed only on "has l and m" mistakes them for
// openable images. Matching the whole axis set also drops the (time, frequency, polarization)
// normalization variables and the beam fit parameters, and it means the optional coordinates are
// never read.
bool HasAllAxes(const zarr_metadata::ArrayMetadata& metadata, const std::array<std::string_view, 5>& axes) {
    return std::all_of(axes.begin(), axes.end(), [&metadata](const auto axis) {
        return zarr_metadata::FindDimensionIndex(metadata, axis).has_value();
    });
}

std::vector<std::string> FindDataGroups(const nlohmann::json& root_attributes, std::string_view image_id) {
    std::vector<std::string> data_groups;
    const auto* groups = ObjectMember(root_attributes, "data_groups");
    if (groups == nullptr || !groups->is_object()) {
        return data_groups;
    }

    for (const auto& [group_name, group] : groups->items()) {
        if (!group.is_object()) {
            continue;
        }
        const auto references_image = std::any_of(group.begin(), group.end(), [&](const nlohmann::json& value) {
            return value.is_string() && value.get<std::string>() == image_id;
        });
        if (references_image) {
            data_groups.push_back(group_name);
        }
    }
    return data_groups;
}

std::vector<AxisDescriptor> DescribeAxes(const Store& store, const zarr_metadata::ArrayMetadata& image) {
    constexpr std::array<std::string_view, kXradioImageAxisOrder.size()> logical_axis_names{"l", "m", "frequency",
                                                                                            "polarization", "time"};
    std::vector<AxisDescriptor> axes;
    axes.reserve(logical_axis_names.size());
    for (std::size_t logical = 0; logical < logical_axis_names.size(); ++logical) {
        const auto name = logical_axis_names.at(logical);
        const auto index = zarr_metadata::FindDimensionIndex(image, name);
        if (!index) {
            continue;
        }
        std::string unit;
        const auto& coordinate_metadata = store.ReadArrayMetadata(name);
        if (coordinate_metadata) {
            unit = AttributeString(coordinate_metadata.value().attributes, "units");
        }
        axes.push_back(AxisDescriptor{std::string(name), kXradioImageAxisOrder.at(logical), image.shape.at(*index),
                                      std::move(unit), *index});
    }
    return axes;
}

// A dataset stores each coordinate once and every image references it by dimension name, so an
// image whose own extent disagrees with the coordinate it names cannot be described with it. The
// probe checks this against the default image, which is the only one it looks at; a dataset may
// hold several images, and this is where the rest of them are held to the same requirement rather
// than being described with a coordinate vector of another image's length.
Result<void> RequireMatchingCoordinates(const Store& store, const zarr_metadata::ArrayMetadata& image,
                                        std::string_view image_id) {
    const auto rank = std::min(image.dimension_names.size(), image.shape.size());
    for (std::size_t axis = 0; axis < rank; ++axis) {
        const auto& name = image.dimension_names.at(axis);
        const auto& coordinate = store.ReadArrayMetadata(name);
        // A coordinate the dataset does not carry at all is the probe's business, and reading one
        // reports its own absence. This is only about the two disagreeing.
        if (!coordinate) {
            continue;
        }
        if (coordinate.value().shape.size() != 1 || coordinate.value().shape.front() != image.shape.at(axis)) {
            return MakeError(ErrorCode::invalid_metadata,
                             "Image dimension '" + name + "' is not the length of the coordinate of that name",
                             std::string(image_id));
        }
    }
    return {};
}

Result<std::vector<double>> ReadNumericCoordinate(const Store& store, std::string_view name) {
    const auto& metadata = store.ReadNodeMetadata(name);
    if (!metadata) {
        if (metadata.error().code == ErrorCode::not_found) {
            return std::vector<double>{};
        }
        return metadata.error();
    }
    return store.ReadNumericArray(name);
}

std::optional<SpectralCoordinate> DescribeSpectralCoordinate(const Store& store,
                                                             const std::vector<double>& frequency_values,
                                                             ImageDescriptor& descriptor) {
    if (frequency_values.empty()) {
        return std::nullopt;
    }

    SpectralCoordinate spectral;
    spectral.channel_frequencies = frequency_values;
    nlohmann::json frequency_attributes = nlohmann::json::object();
    if (const auto& frequency_metadata = store.ReadNodeMetadata("frequency"); frequency_metadata) {
        if (const auto* attributes = ObjectMember(frequency_metadata.value(), "attributes");
            attributes != nullptr && attributes->is_object()) {
            frequency_attributes = *attributes;
        }
    }

    spectral.unit = AttributeString(frequency_attributes, "units");
    if (const auto* reference_frequency = ObjectMember(frequency_attributes, "reference_frequency");
        reference_frequency != nullptr && reference_frequency->is_object()) {
        if (const auto* attributes = ObjectMember(*reference_frequency, "attrs");
            attributes != nullptr && attributes->is_object()) {
            if (spectral.unit.empty()) {
                spectral.unit = AttributeString(*attributes, "units");
            }
            spectral.system = Upper(AttributeString(*attributes, "observer"));
        }
    }

    double reference_value = spectral.channel_frequencies.front();
    if (const auto* reference_frequency = ObjectMember(frequency_attributes, "reference_frequency");
        reference_frequency != nullptr && reference_frequency->is_object()) {
        if (const auto* data = ObjectMember(*reference_frequency, "data"); data != nullptr) {
            if (const auto value = AttributeNumber(*data)) {
                reference_value = *value;
            }
        }
    }
    if (const auto* rest_frequency = ObjectMember(frequency_attributes, "rest_frequency");
        rest_frequency != nullptr && rest_frequency->is_object()) {
        if (const auto* data = ObjectMember(*rest_frequency, "data"); data != nullptr && data->is_number()) {
            spectral.rest_frequency = data->get<double>();
        }
    }

    auto fit = FitSpectralAxis(spectral.channel_frequencies, reference_value, "frequency");
    spectral.reference_pixel = fit.reference_pixel;
    spectral.reference_value = fit.reference_value;
    spectral.increment = fit.increment;
    AppendDiagnostics(descriptor, std::move(fit.diagnostics));
    return spectral;
}

std::optional<TemporalCoordinate> DescribeTemporalCoordinate(const Store& store, std::vector<double> values) {
    if (values.empty()) {
        return std::nullopt;
    }

    TemporalCoordinate temporal;
    temporal.values = std::move(values);
    if (const auto& metadata = store.ReadNodeMetadata("time"); metadata) {
        if (const auto* attributes = ObjectMember(metadata.value(), "attributes");
            attributes != nullptr && attributes->is_object()) {
            temporal.unit = AttributeString(*attributes, "units");
            temporal.scale = Upper(AttributeString(*attributes, "scale"));
            temporal.format = Upper(AttributeString(*attributes, "format"));
        }
    }
    return temporal;
}

// Which variables of this store are images this profile will open. Only InspectImages calls it:
// deciding whether the store matches means enumerating it, and the answer to both questions comes
// back together so that nothing has to enumerate twice or cache the result.
Result<ImageDiscovery> DiscoverImages(const Store& store) {
    const auto& nodes = store.ListNodes();
    if (!nodes) {
        return nodes.error();
    }

    ImageDiscovery result;
    for (const auto& node : nodes.value()) {
        const auto& array_result = store.ReadArrayMetadata(node);
        if (!array_result) {
            continue;
        }
        const auto& array = array_result.value();
        if (IsFlag(array)) {
            continue;
        }

        if (HasAllAxes(array, kSkyAxes)) {
            if (zarr_metadata::IsRealDataType(array.data_type)) {
                result.images.push_back(ImageEntry{node, true, {}});
            } else {
                const Diagnostic diagnostic{"unsupported_data_type", "Complex sky-plane variables are not openable",
                                            node};
                result.images.push_back(ImageEntry{node, false, {diagnostic}});
                result.diagnostics.push_back(diagnostic);
            }
        } else if (HasAllAxes(array, kApertureAxes)) {
            const Diagnostic diagnostic{"unsupported_coordinate_plane", "Aperture-plane variables are not openable",
                                        node};
            result.images.push_back(ImageEntry{node, false, {diagnostic}});
            result.diagnostics.push_back(diagnostic);
        }
    }

    const auto sort_images = [](std::vector<ImageEntry>& images) {
        std::sort(images.begin(), images.end(), [](const auto& left, const auto& right) {
            const int left_rank = KnownImageRank(left.id);
            const int right_rank = KnownImageRank(right.id);
            return left_rank == right_rank ? left.id < right.id : left_rank < right_rank;
        });
    };
    sort_images(result.images);
    const auto readable = std::find_if(result.images.begin(), result.images.end(),
                                       [](const ImageEntry& image) { return image.readable; });
    if (readable != result.images.end()) {
        result.default_image_id = readable->id;
    }
    return result;
}

}  // namespace

Result<SchemaInspection> InspectImages(const Store& store) {
    ProbeReport report(store, "image dataset");

    const auto& root_attributes = store.RootAttributes();
    auto discovery = DiscoverImages(store);
    if (!discovery) {
        return discovery.error();
    }
    report.SetDiagnostics(discovery.value().diagnostics);

    // What was found goes back with what was decided from it, whichever way the decision went.
    const auto finish = [&](SchemaMatchKind kind) -> Result<SchemaInspection> {
        auto probe = report.Finish(kind, std::string(kVersion));
        if (!probe) {
            return probe.error();
        }
        return SchemaInspection{std::move(probe.value()), std::move(discovery.value())};
    };

    if (!discovery.value().default_image_id) {
        // A valid Zarr store without an image that this profile can open is a non-match. The
        // discovery diagnostics still explain why variables such as complex or aperture-plane
        // arrays were not openable.
        return finish(SchemaMatchKind::no_match);
    }

    // Once discovery found an openable image, validate the metadata needed by the image reader.
    const auto& first_image = *discovery.value().default_image_id;
    const auto& array_result = store.ReadArrayMetadata(first_image);
    if (report.RequireArrayMetadata(array_result, first_image)) {
        RequirePresentCoordinates(report, array_result.value());
        if (report.ok() && HasAttribute(root_attributes, "coordinate_system_info")) {
            report.RequireCoordinateSystem(root_attributes);
        }
    }
    return finish(report.ok() ? SchemaMatchKind::match : SchemaMatchKind::invalid);
}

Result<ImageDescriptor> DescribeImage(const Store& store, std::string_view image_id) {
    const auto& array_result = store.ReadArrayMetadata(image_id);
    if (!array_result) {
        return array_result.error();
    }
    const auto& image = array_result.value();
    if (IsFlag(image) || !zarr_metadata::FindDimensionIndex(image, "l") ||
        !zarr_metadata::FindDimensionIndex(image, "m") || !zarr_metadata::IsRealDataType(image.data_type)) {
        return MakeError(ErrorCode::unsupported_data_type, "Image variable is not an openable sky-plane image",
                         std::string(image_id));
    }

    if (auto matching = RequireMatchingCoordinates(store, image, image_id); !matching) {
        return matching.error();
    }

    ImageDescriptor descriptor;
    descriptor.id = std::string(image_id);
    descriptor.stored_type = zarr_metadata::ParseDataType(image.data_type);
    descriptor.unit = AttributeString(image.attributes, "units");
    // XRADIO writes the role on the variable's own "type" attribute, lowercased ("sky", "model",
    // "residual"); the same attribute spells "flag" for pixel masks, which are never images. The v2
    // schema proposes a separate "image_type" attribute that no released XRADIO writes yet, so it is
    // only consulted as a forward-compatible fallback.
    const std::string declared_role = AttributeString(image.attributes, "type");
    descriptor.image_role = declared_role == "flag" ? std::string{} : declared_role;
    if (descriptor.image_role.empty()) {
        descriptor.image_role = AttributeString(image.attributes, "image_type");
    }

    const auto& root_attrs = store.RootAttributes();
    descriptor.data_groups = FindDataGroups(root_attrs, image_id);
    descriptor.axes = DescribeAxes(store, image);

    std::vector<double> l_values;
    std::vector<double> m_values;
    std::vector<double> frequency_values;
    std::vector<double> time_values;
    auto l_result = ReadNumericCoordinate(store, "l");
    if (!l_result) {
        return l_result.error();
    }
    l_values = std::move(l_result.value());
    auto m_result = ReadNumericCoordinate(store, "m");
    if (!m_result) {
        return m_result.error();
    }
    m_values = std::move(m_result.value());

    // A direction axis is linear by construction, so its increment is reported even when the samples
    // are not evenly spaced; the fit says so in a diagnostic rather than withholding the value.
    if (auto direction = DescribeDirection(root_attrs, l_values, m_values, descriptor.diagnostics); direction) {
        descriptor.direction = std::move(direction);
    }

    auto frequency_result = ReadNumericCoordinate(store, "frequency");
    if (!frequency_result) {
        return frequency_result.error();
    }
    frequency_values = std::move(frequency_result.value());
    if (auto spectral = DescribeSpectralCoordinate(store, frequency_values, descriptor); spectral) {
        descriptor.spectral = std::move(spectral);
    }

    const auto& polarization_metadata = store.ReadNodeMetadata("polarization");
    if (polarization_metadata) {
        auto pol_labels = store.ReadStringArray1D("polarization");
        if (!pol_labels) {
            return pol_labels.error();
        }
        PolarizationCoordinate pol;
        pol.labels = std::move(pol_labels.value());
        descriptor.polarization = std::move(pol);
    }

    auto time_result = ReadNumericCoordinate(store, "time");
    if (!time_result) {
        return time_result.error();
    }
    time_values = std::move(time_result.value());
    if (auto temporal = DescribeTemporalCoordinate(store, std::move(time_values)); temporal) {
        descriptor.temporal = std::move(temporal);
    }

    // Parse Observation & Telescope Metadata from the selected image's attributes.
    descriptor.observation = DescribeObservation(image);

    auto layout_res = store.ReadStorageLayout(image_id);
    if (layout_res) {
        descriptor.storage = std::move(layout_res.value());
    }

    auto pixel_mask = DetermineFlag(store, image, image_id, descriptor.diagnostics);
    if (!pixel_mask) {
        return pixel_mask.error();
    }
    descriptor.pixel_mask_id = pixel_mask.value();
    descriptor.has_pixel_mask = !descriptor.pixel_mask_id.empty();

    return descriptor;
}

Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id) {
    const auto& sky_meta = store.ReadNodeMetadata(image_id);
    if (!sky_meta) {
        return sky_meta.error();
    }
    std::string beam_array_name;
    if (const auto* attributes = ObjectMember(sky_meta.value(), "attributes"); attributes != nullptr) {
        beam_array_name = AttributeString(*attributes, "beam_fit_params");
    }
    if (beam_array_name.empty()) {
        // No beam array associated
        return std::vector<Beam>{};
    }

    const auto& beam_metadata = store.ReadArrayMetadata(beam_array_name);
    if (!beam_metadata) {
        return beam_metadata.error();
    }

    // The image named a beam table, so what the table needs to be read is required rather than
    // optional. A label array that cannot be read -- absent, malformed, or written with a codec
    // this build does not carry -- used to be discarded here, and the parameter indices it did not
    // yield then returned an empty beam list: the answer for an image that has no beam at all,
    // which this one is not.
    auto parameter_labels = store.ReadStringArray1D("beam_params_label");
    if (!parameter_labels) {
        return parameter_labels.error();
    }

    auto values = store.ReadNumericArray(beam_array_name);
    if (!values) {
        return values.error();
    }

    return DescribeBeams(beam_metadata.value(), beam_array_name, parameter_labels.value(), values.value());
}

}  // namespace carta::zarr::internal::xradio
