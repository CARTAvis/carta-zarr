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

std::string RejectionMessage(const std::vector<Diagnostic>& diagnostics, std::string_view fallback) {
    return diagnostics.empty() ? std::string(fallback) : diagnostics.front().message;
}

Result<void> RequireOpenableDataset(const ProbeResult& probe, std::string_view location) {
    if (probe.kind != ProbeKind::supported_dataset) {
        // Malformed and unrecognised are different answers, and a consumer acts on the difference:
        // one is a file to complain about, the other is a file this library is not for.
        const ErrorCode code = probe.kind == ProbeKind::invalid_dataset ? ErrorCode::invalid_metadata
                                                                       : ErrorCode::unsupported_schema;
        return Error{code, RejectionMessage(probe.diagnostics, "No built-in schema profile matched the Zarr store"),
                     std::string(location)};
    }
    // A store the profile recognised and found nothing openable in. Not a rejection -- the profile
    // matched -- so there are no diagnostics to report, and a consumer opening this would get a
    // dataset it can do nothing with.
    if (probe.images.empty()) {
        return Error{ErrorCode::invalid_metadata, "Supported schema has no image variables",
                     std::string(location)};
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
        return Error{ErrorCode::unsupported_schema,
                     "No built-in profile exists for schema " + std::string(schema_id)};
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

Result<DescribedImage> SchemaProfile::Describe(const Store& store, std::string_view image_id) const {
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
    // What a profile that did not match had to say about the store anyway: a store of nothing but
    // complex or aperture-plane variables is one the profile recognised the images of and will not
    // open, and that is the reason to give, not "nothing matched". A profile that found nothing it
    // knew says nothing, so a store this library is simply not for still carries no diagnostic.
    // The first profile to say anything is the one heard, as the first invalid one is below.
    std::vector<Diagnostic> unmatched;
    for (const auto& entry : SchemaProfile::BuiltIn()) {
        auto inspection = entry.inspect(store);
        if (!inspection) {
            return inspection.error();
        }
        auto& probe = inspection.value().probe;
        if (probe.kind == SchemaMatchKind::match) {
            matches.push_back(std::move(inspection.value()));
        } else if (probe.kind == SchemaMatchKind::invalid) {
            invalid.push_back(std::move(probe));
        } else if (unmatched.empty()) {
            unmatched = std::move(probe.diagnostics);
        }
    }

    if (matches.size() > 1) {
        // Unreachable while BuiltIn holds one entry, and the only thing here that a second profile
        // turns on. ErrorCode::ambiguous_schema is not emitted anywhere yet, because this answers
        // with a ProbeResult rather than an Error.
        result.kind = ProbeKind::invalid_dataset;
        result.diagnostics.push_back(Diagnostic{
            DiagnosticCode::ambiguous_schema, "More than one built-in schema profile matched the Zarr store", {}});
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
        // Diagnostics alone, as ProbeResult promises for this kind: naming the schema would say the
        // store is one, which is what not matching denies.
        result.kind = ProbeKind::zarr_without_supported_schema;
        result.diagnostics = std::move(unmatched);
    }
    return result;
}

}  // namespace carta::zarr::internal
