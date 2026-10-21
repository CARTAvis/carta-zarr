/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "schema_version.h"

#include "attributes.h"

namespace carta::zarr::internal::xradio {

namespace {

// Why a store naming 0.0.2 does not conform, or nothing when it does. A store's version is only a
// claim, so this checks the three things schema 0.0.2 requires that reading the store depends on,
// rather than the whole schema:
//
// - It has data groups.
// - The sky its data groups name, if they name one, is an array. Whether the sky is typed "sky" is
//   checked by qualification, as for every image a data group names (ADR 0001). The schema does not
//   require a sky: a store of aperture-plane images, or of a point spread function alone, is well
//   formed without one.
// - Its frequency coordinate, if it has one, gives a frame and units. Schema 0.0.2 requires these of
//   every store, with or without a sky.
std::optional<Diagnostic> WhyNot002(const Store& store) {
    const auto not_002 = [](const std::string& why, std::string node_path) {
        return Diagnostic{DiagnosticCode::invalid_metadata,
                          "The store names XRADIO image schema " + std::string(kSchemaVersion) + ", but " + why,
                          std::move(node_path)};
    };
    const auto* const groups = MemberObject(store.RootAttributes(), "data_groups");
    if (groups == nullptr) {
        return not_002("it has no data_groups", "/attributes/data_groups");
    }
    // The base group's sky, which is the one XRADIO writes first; another group's when it has none.
    std::string sky;
    if (const auto* const base = MemberObject(*groups, "base")) {
        sky = AttributeString(*base, "sky");
    }
    for (auto group = groups->begin(); sky.empty() && group != groups->end(); ++group) {
        sky = AttributeString(group.value(), "sky");
    }
    if (!sky.empty()) {
        const auto& image = store.ReadArrayMetadata(sky);
        if (!image) {
            return not_002("its data groups name " + sky + " as the sky, which is not an array", sky);
        }
    }
    if (const auto& frequency = store.ReadArrayMetadata("frequency"); frequency) {
        const auto& attributes = frequency.value().attributes;
        if (AttributeString(attributes, "frame").empty() || AttributeString(attributes, "units").empty()) {
            return not_002("its frequency coordinate has no frame or no units", "frequency");
        }
    }
    return std::nullopt;
}

}  // namespace

SchemaVersionCheck CheckSchemaVersion(const Store& store) {
    SchemaVersionCheck checked;
    const auto* const named = Member(store.RootAttributes(), "schema_version");
    if (named == nullptr) {
        checked.refusal = Diagnostic{DiagnosticCode::unsupported_schema_version,
                                     "The store names no XRADIO image schema version, as XRADIO 1.2.4 and earlier "
                                     "wrote; this library reads schema " +
                                         std::string(kSchemaVersion) + ". Write it again with XRADIO 1.2.5 or later",
                                     "/attributes/schema_version"};
        return checked;
    }
    if (named->is_string()) {
        checked.version = named->get<std::string>();
    }
    if (checked.version != kSchemaVersion) {
        const std::string named_as = named->is_string() ? checked.version : named->dump();
        checked.refusal = Diagnostic{DiagnosticCode::unsupported_schema_version,
                                     "XRADIO image schema " + named_as + " is not one this library reads; it reads " +
                                         std::string(kSchemaVersion),
                                     "/attributes/schema_version"};
        return checked;
    }
    checked.refusal = WhyNot002(store);
    return checked;
}

}  // namespace carta::zarr::internal::xradio
