/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "coordinates.h"

#include "../../zarr/string_array.h"
#include "attributes.h"

#include <algorithm>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// Where a numeric coordinate's samples go, by the axis it is.
std::vector<double>* NumericValuesOf(CoordinateValues& values, AxisRole role) {
    switch (role) {
        case AxisRole::spatial_x:
            return &values.l;
        case AxisRole::spatial_y:
            return &values.m;
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
    const auto& attributes = metadata.value().attributes;
    auto unit = AttributeString(attributes, "units");
    if (unit.empty() && coordinate.role == AxisRole::spectral) {
        if (const auto* const reference = MemberObject(attributes, "reference_frequency")) {
            if (const auto* const measure = MemberObject(*reference, "attrs")) {
                unit = AttributeString(*measure, "units");
            }
        }
    }
    return unit;
}

std::vector<AxisDescriptor> DescribeAxes(const Store& store, const zarr_metadata::ArrayMetadata& image) {
    std::vector<AxisDescriptor> axes;
    axes.reserve(AxisCount(Plane::sky));
    for (const auto& coordinate : kCoordinates) {
        if (!OnPlane(coordinate, Plane::sky)) {
            continue;
        }
        const auto index = zarr_metadata::FindDimensionIndex(image, coordinate.name);
        if (!index) {
            continue;
        }
        axes.push_back(AxisDescriptor{std::string(coordinate.name), coordinate.role, image.shape.at(*index),
                                      CoordinateUnit(store, coordinate), *index});
    }
    return axes;
}

Result<CoordinateValues> ReadCoordinateValues(const Store& store) {
    // Every one of them, because every one is required: the probe refused a dataset missing one
    // before an image of it could be described. A coordinate that cannot be read fails describing
    // rather than leaving a hole -- the polarization labels used to be skipped when their metadata
    // failed to read for any reason, an I/O error included.
    CoordinateValues values;
    for (const auto& coordinate : kCoordinates) {
        if (!OnPlane(coordinate, Plane::sky)) {
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
        *NumericValuesOf(values, coordinate.role) = samples.value().values();
    }
    return values;
}

}  // namespace carta::zarr::internal::xradio
