/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "coordinates.h"

#include "../../zarr/string_array.h"
#include "attributes.h"

#include <algorithm>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// Where a numeric coordinate's samples go, by the axis it is on its plane.
std::vector<double>* NumericValuesOf(CoordinateValues& values, const Coordinate& coordinate) {
    const bool sky = coordinate.plane == Plane::sky;
    switch (coordinate.role) {
        case AxisRole::spatial_x:
            return sky ? &values.l : &values.u;
        case AxisRole::spatial_y:
            return sky ? &values.m : &values.v;
        case AxisRole::spectral:
            return &values.frequency;
        case AxisRole::time:
            return &values.time;
        case AxisRole::polarization:
        case AxisRole::other:
            return nullptr;
    }
    return nullptr;
}

}  // namespace

std::size_t AxisCount(Plane plane) noexcept {
    return static_cast<std::size_t>(
        std::count_if(kCoordinates.begin(), kCoordinates.end(),
                      [plane](const Coordinate& coordinate) { return OnPlane(coordinate, plane); }));
}

bool CarriesPlane(const zarr_metadata::ArrayMetadata& variable, Plane plane) {
    return std::all_of(kCoordinates.begin(), kCoordinates.end(), [&](const Coordinate& coordinate) {
        return !OnPlane(coordinate, plane) || zarr_metadata::FindDimensionIndex(variable, coordinate.name).has_value();
    });
}

std::optional<Plane> ImagePlane(std::string_view type) {
    // Schema 0.0.2's image arrays, by the dimensions each is defined over.
    static constexpr std::array<std::string_view, 4> kSky{"sky", "point_spread_function", "primary_beam", "mask"};
    static constexpr std::array<std::string_view, 3> kAperture{"aperture", "visibility", "uv_sampling"};
    if (std::find(kSky.begin(), kSky.end(), type) != kSky.end()) {
        return Plane::sky;
    }
    if (std::find(kAperture.begin(), kAperture.end(), type) != kAperture.end()) {
        return Plane::aperture;
    }
    return std::nullopt;
}

std::optional<Plane> PlaneOf(const zarr_metadata::ArrayMetadata& variable) {
    const auto plane = ImagePlane(AttributeString(variable.attributes, "type"));
    if (!plane || !CarriesPlane(variable, *plane)) {
        return std::nullopt;
    }
    return plane;
}

Result<void> CheckCoordinate(const zarr_metadata::ArrayMetadata& array, const Coordinate& coordinate) {
    const std::string node(coordinate.name);
    // One name, so one dimension: ParseArrayMetadata holds dimension_names to the rank.
    if (array.dimension_names.size() != 1 || array.dimension_names.front() != coordinate.name) {
        return Error{ErrorCode::invalid_metadata, "Coordinate shape or dimension name does not match image dataset",
                     node};
    }
    const bool typed = coordinate.kind == CoordinateKind::labels ? zarr_metadata::IsFixedLengthUtf32(array)
                                                                 : zarr_metadata::IsRealDataType(array.data_type);
    if (!typed) {
        return Error{ErrorCode::unsupported_data_type, "Coordinate array has an unsupported data type", node};
    }
    if (coordinate.kind == CoordinateKind::labels) {
        return zarr_metadata::CheckFixedLengthUtf32StringArray(array, node);
    }
    return {};
}

std::optional<Diagnostic> ExtentDisagreement(const Store& store, const zarr_metadata::ArrayMetadata& image,
                                             std::string_view node) {
    const auto rank = std::min(image.dimension_names.size(), image.shape.size());
    for (std::size_t axis = 0; axis < rank; ++axis) {
        const auto& name = image.dimension_names.at(axis);
        const auto& coordinate = store.ReadArrayMetadata(name);
        if (!coordinate) {
            continue;
        }
        if (coordinate.value().shape.size() != 1 || coordinate.value().shape.front() != image.shape.at(axis)) {
            return Diagnostic{DiagnosticCode::invalid_metadata,
                              "Image dimension '" + name + "' is not the length of the coordinate of that name",
                              std::string(node)};
        }
    }
    return std::nullopt;
}

std::string CoordinateUnit(const Store& store, const Coordinate& coordinate) {
    const auto& metadata = store.ReadArrayMetadata(coordinate.name);
    if (!metadata) {
        return {};
    }
    return AttributeString(metadata.value().attributes, "units");
}

std::vector<AxisDescriptor> DescribeAxes(const zarr_metadata::ArrayMetadata& image) {
    const auto plane = PlaneOf(image).value_or(Plane::sky);
    std::vector<AxisDescriptor> axes;
    axes.reserve(AxisCount(plane));
    for (const auto& coordinate : kCoordinates) {
        if (!OnPlane(coordinate, plane)) {
            continue;
        }
        const auto index = zarr_metadata::FindDimensionIndex(image, coordinate.name);
        if (!index) {
            continue;
        }
        axes.push_back(AxisDescriptor{std::string(coordinate.name), coordinate.role, image.shape.at(*index), *index});
    }
    return axes;
}

Result<CoordinateValues> ReadCoordinateValues(const Store& store, Plane plane) {
    // Every one of the plane's, because every one is required: the probe refuses a dataset missing
    // one before an image of it can be described. A coordinate that cannot be read, the polarization
    // labels included, fails the description rather than leaving a hole.
    CoordinateValues values;
    for (const auto& coordinate : kCoordinates) {
        if (!OnPlane(coordinate, plane)) {
            continue;
        }
        if (coordinate.kind == CoordinateKind::labels) {
            auto labels = store.ReadStringArray1D(coordinate.name);
            if (!labels) {
                return labels.error();
            }
            values.polarization = std::move(labels.value());
            continue;
        }
        const auto& samples = store.ReadNumericArray(coordinate.name);
        if (!samples) {
            return samples.error();
        }
        *NumericValuesOf(values, coordinate) = samples.value().values();
    }
    return values;
}

}  // namespace carta::zarr::internal::xradio
