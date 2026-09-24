/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "probe_report.h"

#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

}  // namespace

namespace {

// A metadata read that failed, as the diagnostic the probe reports it by. A data type this library
// does not read is that; anything else a read can fail with -- a codec or a Zarr version it does not
// read, or the store not answering -- is metadata the probe could not use, and the message says which.
//
// Every code is listed, with no default, so that an ErrorCode added later is a warning here until
// someone says what a probe should call it.
DiagnosticCode DiagnosticCodeFor(ErrorCode code) {
    switch (code) {
        case ErrorCode::unsupported_data_type:
            return DiagnosticCode::unsupported_data_type;
        case ErrorCode::ambiguous_schema:
            return DiagnosticCode::ambiguous_schema;
        case ErrorCode::invalid_metadata:
        case ErrorCode::not_found:
        case ErrorCode::not_zarr:
        case ErrorCode::unsupported_transport:
        case ErrorCode::unsupported_zarr_version:
        case ErrorCode::unsupported_schema:
        case ErrorCode::unsupported_schema_version:
        case ErrorCode::invalid_argument:
        case ErrorCode::unsupported_codec:
        case ErrorCode::invalid_slice:
        case ErrorCode::buffer_too_small:
        case ErrorCode::io_error:
        case ErrorCode::decode_error:
        case ErrorCode::cancelled:
        case ErrorCode::not_implemented:
            return DiagnosticCode::invalid_metadata;
    }
    return DiagnosticCode::invalid_metadata;
}

}  // namespace

ProbeReport::ProbeReport(const Store& store, std::string profile_name)
    : _store(&store), _profile_name(std::move(profile_name)) {}

bool ProbeReport::ok() const noexcept {
    return _met && !_error.has_value();
}

void ProbeReport::AddDiagnostic(DiagnosticCode code, std::string message, std::string node_path) {
    _diagnostics.push_back(Diagnostic{code, std::move(message), std::move(node_path)});
}

void ProbeReport::SetDiagnostics(std::vector<Diagnostic> diagnostics) {
    _diagnostics = std::move(diagnostics);
}

const std::vector<Diagnostic>& ProbeReport::diagnostics() const noexcept {
    return _diagnostics;
}

bool ProbeReport::Fail(DiagnosticCode code, std::string message, std::string node_path) {
    AddDiagnostic(code, std::move(message), std::move(node_path));
    _met = false;
    return false;
}

bool ProbeReport::RequireArrayMetadata(const Result<zarr::ArrayMetadata>& metadata, std::string_view node) {
    if (!ok()) {
        return false;
    }
    if (metadata) {
        return true;
    }
    return Fail(DiagnosticCodeFor(metadata.error().code), metadata.error().message, std::string(node));
}

bool ProbeReport::RequireCoordinateOf(const zarr::ArrayMetadata& image, std::string_view axis, CoordinateKind kind) {
    if (!ok()) {
        return false;
    }
    const auto index = zarr_metadata::FindDimensionIndex(image, axis);
    if (!index) {
        // An axis the image does not carry has nothing to validate against.
        return true;
    }

    const std::string node(axis);
    const auto& metadata_result = _store->ReadNodeMetadata(axis);
    if (!metadata_result) {
        if (metadata_result.error().code == ErrorCode::not_found) {
            return Fail(DiagnosticCode::invalid_metadata, "Missing required coordinate array", node);
        }
        _error = metadata_result.error();
        return false;
    }

    // ReadArrayMetadata reuses both the raw metadata and parsed array metadata caches. Keeping the
    // raw read above preserves the distinction between missing metadata and other I/O errors.
    const auto& array_result = _store->ReadArrayMetadata(axis);
    if (!RequireArrayMetadata(array_result, node)) {
        return false;
    }

    const auto& coordinate = array_result.value();
    // The length the image expects is not checked here. Agreeing with it is part of being an image
    // this profile opens, so the default image -- which is the first one that qualified -- agrees
    // with every coordinate it can read before the probe ever sees it. What is left is whether the
    // coordinate is well formed in itself, which is a fact about the dataset.
    if (coordinate.shape.size() != 1 || coordinate.dimension_names.size() != 1 ||
        coordinate.dimension_names.front() != axis) {
        return Fail(DiagnosticCode::invalid_metadata,
                    "Coordinate shape or dimension name does not match " + _profile_name, node);
    }
    const bool typed = kind == CoordinateKind::labels ? zarr_metadata::IsFixedLengthUtf32(coordinate)
                                                      : zarr_metadata::IsRealDataType(coordinate.data_type);
    if (!typed) {
        return Fail(DiagnosticCode::unsupported_data_type, "Coordinate array has an unsupported data type", node);
    }
    return true;
}

bool ProbeReport::RequireCoordinateSystem(const nlohmann::json& root_attributes) {
    if (!ok()) {
        return false;
    }
    const std::string node = "/attributes/coordinate_system_info";
    if (!root_attributes.is_object() || !root_attributes.contains("coordinate_system_info") ||
        !root_attributes.at("coordinate_system_info").is_object()) {
        return Fail(DiagnosticCode::invalid_metadata, "XRADIO requires coordinate_system_info metadata", node);
    }
    const auto& coordinate = root_attributes.at("coordinate_system_info");
    if (!coordinate.contains("projection") || !coordinate.at("projection").is_string() ||
        coordinate.at("projection").get<std::string>().empty()) {
        return Fail(DiagnosticCode::invalid_metadata, "coordinate_system_info requires a projection", node);
    }
    if (!coordinate.contains("reference_direction") || !coordinate.at("reference_direction").is_object() ||
        !coordinate.at("reference_direction").contains("data") ||
        !zarr_metadata::IsNumericVector(coordinate.at("reference_direction").at("data"), 2)) {
        return Fail(DiagnosticCode::invalid_metadata, "coordinate_system_info requires a two-value reference direction",
                    node);
    }
    if (!coordinate.contains("native_pole_direction") || !coordinate.at("native_pole_direction").is_object() ||
        !coordinate.at("native_pole_direction").contains("data") ||
        !zarr_metadata::IsNumericVector(coordinate.at("native_pole_direction").at("data"), 2)) {
        return Fail(DiagnosticCode::invalid_metadata,
                    "coordinate_system_info requires a two-value native pole direction", node);
    }
    if (!coordinate.contains("pixel_coordinate_transformation_matrix") ||
        !zarr_metadata::IsNumericMatrix(coordinate.at("pixel_coordinate_transformation_matrix"), 2, 2)) {
        return Fail(DiagnosticCode::invalid_metadata,
                    "coordinate_system_info requires a 2x2 pixel transformation matrix", node);
    }
    return true;
}

std::string ProbeReport::InvalidMessage() const {
    std::string message = "XRADIO " + _profile_name + " metadata is invalid";
    if (!_diagnostics.empty()) {
        message += ": " + _diagnostics.front().message;
    }
    return message;
}

Result<SchemaProbeResult> ProbeReport::Finish(SchemaMatchKind kind, std::string schema_version) const {
    if (_error) {
        return *_error;
    }
    return SchemaProbeResult{kind, std::string(kXradioImageSchema), std::move(schema_version), _diagnostics};
}

}  // namespace carta::zarr::internal::xradio
