/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "profile.h"

#include "xradio/image.h"

#include "zarr/array_metadata.h"

#include <algorithm>
#include <utility>

namespace carta::zarr::internal {
namespace {

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

}  // namespace

Result<void> RequireOpenable(const std::vector<ImageEntry>& images, std::string_view image_id) {
    const auto found = std::find_if(images.begin(), images.end(),
                                    [&](const ImageEntry& image) { return image.id == image_id; });
    if (found == images.end()) {
        return MakeError(ErrorCode::not_found, "Image variable was not found", std::string(image_id));
    }
    if (found->readable) {
        return {};
    }
    // Whatever the profile said about why it would not open this one is more use than a generic
    // refusal, and it is already in hand.
    const auto message =
        found->diagnostics.empty() ? "Image variable is not openable by this profile" : found->diagnostics.front().message;
    return MakeError(ErrorCode::unsupported_data_type, message, std::string(image_id));
}

std::string RejectionMessage(const std::vector<Diagnostic>& diagnostics, std::string_view fallback) {
    return diagnostics.empty() ? std::string(fallback) : diagnostics.front().message;
}

Result<void> RequireOpenableDataset(const ProbeResult& probe, std::string_view location) {
    if (probe.kind != ProbeKind::supported_dataset) {
        // Malformed and unrecognised are different answers, and a consumer acts on the difference:
        // one is a file to complain about, the other is a file this library is not for.
        const ErrorCode code = probe.kind == ProbeKind::invalid_dataset ? ErrorCode::invalid_metadata
                                                                       : ErrorCode::unsupported_schema;
        return MakeError(code, RejectionMessage(probe.diagnostics, "No built-in schema profile matched the Zarr store"),
                         std::string(location));
    }
    // A store the profile recognised and found nothing openable in. Not a rejection -- the profile
    // matched -- so there are no diagnostics to report, and a consumer opening this would get a
    // dataset it can do nothing with.
    if (probe.images.empty()) {
        return MakeError(ErrorCode::invalid_metadata, "Supported schema has no image variables",
                         std::string(location));
    }
    return {};
}

const std::vector<SchemaProfile::Entry>& SchemaProfile::BuiltIn() {
    static const std::vector<SchemaProfile::Entry> profiles{
        {SchemaProfile::Entry{SchemaId(kXradioImageSchema), &xradio::InspectImages, &xradio::DescribeImage,
                              &xradio::ReadBeams}}};
    return profiles;
}

Result<SchemaProfile> SchemaProfile::For(std::string_view schema_id) {
    const auto& profiles = BuiltIn();
    const auto found = std::find_if(profiles.begin(), profiles.end(),
                                    [&](const Entry& candidate) { return candidate.id == schema_id; });
    if (found == profiles.end()) {
        return MakeError(ErrorCode::unsupported_schema,
                         "No built-in profile exists for schema " + std::string(schema_id));
    }
    return SchemaProfile{*found};
}

Result<SchemaInspection> SchemaProfile::Inspect(const Store& store) const {
    return _entry->inspect(store);
}

Result<SchemaProbeResult> SchemaProfile::Probe(const Store& store) const {
    auto inspection = Inspect(store);
    if (!inspection) {
        return inspection.error();
    }
    return std::move(inspection.value().probe);
}

Result<ImageDiscovery> SchemaProfile::Discover(const Store& store) const {
    auto inspection = Inspect(store);
    if (!inspection) {
        return inspection.error();
    }
    return std::move(inspection.value().discovery);
}

Result<ImageDescriptor> SchemaProfile::Describe(const Store& store, std::string_view image_id) const {
    auto discovery = Discover(store);
    if (!discovery) {
        return discovery.error();
    }
    if (auto openable = RequireOpenable(discovery.value().images, image_id); !openable) {
        return openable.error();
    }
    return _entry->describe(store, image_id);
}

Result<ImageDescriptor> SchemaProfile::DescribeVerified(const Store& store, std::string_view image_id) const {
    return _entry->describe(store, image_id);
}

Result<std::vector<Beam>> SchemaProfile::ReadBeams(const Store& store, std::string_view image_id) const {
    return _entry->read_beams(store, image_id);
}

Result<ProbeResult> ProbeStore(const Store& store) {
    ProbeResult result;
    // Each profile is asked once, and answers with what it found as well as what it decided. The
    // store used to be enumerated twice for a match -- once to probe and once to discover -- which
    // is what a cache inside Store was there to hide.
    std::vector<SchemaInspection> matches;
    std::vector<SchemaProbeResult> invalid;
    for (const auto& entry : SchemaProfile::BuiltIn()) {
        auto inspection = entry.inspect(store);
        if (!inspection) {
            return inspection.error();
        }
        if (inspection.value().probe.kind == SchemaMatchKind::match) {
            matches.push_back(std::move(inspection.value()));
        } else if (inspection.value().probe.kind == SchemaMatchKind::invalid) {
            invalid.push_back(std::move(inspection.value().probe));
        }
    }

    if (matches.size() > 1) {
        // Unreachable while BuiltIn holds one entry, and the only thing here that a second profile
        // turns on. The code name comes from the enumerator rather than from a literal beside it so
        // that the two cannot drift apart in the meantime -- ErrorCode::ambiguous_schema is not
        // emitted anywhere yet, because this answers with a ProbeResult rather than an Error.
        result.kind = ProbeKind::invalid_dataset;
        result.diagnostics.push_back(Diagnostic{zarr::ErrorCodeName(ErrorCode::ambiguous_schema),
                                                "More than one built-in schema profile matched the Zarr store",
                                                {}});
    } else if (matches.size() == 1) {
        auto& match = matches.front();
        result.kind = ProbeKind::supported_dataset;
        result.schema_id = std::move(match.probe.schema_id);
        result.schema_version = std::move(match.probe.schema_version);
        result.images = std::move(match.discovery.images);
        result.default_image_id = std::move(match.discovery.default_image_id);
        result.diagnostics = std::move(match.probe.diagnostics);
    } else if (!invalid.empty()) {
        const auto& invalid_result = invalid.front();
        result.kind = ProbeKind::invalid_dataset;
        result.schema_id = invalid_result.schema_id;
        result.schema_version = invalid_result.schema_version;
        result.diagnostics = invalid_result.diagnostics;
    } else {
        result.kind = ProbeKind::zarr_without_supported_schema;
    }
    return result;
}

}  // namespace carta::zarr::internal
