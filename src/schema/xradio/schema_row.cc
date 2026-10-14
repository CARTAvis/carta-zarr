/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "schema_row.h"

#include "attributes.h"

#include <algorithm>
#include <optional>

namespace carta::zarr::internal::xradio {

namespace {

// Why a store naming 0.0.2 is not one, or nothing when it is. The three things schema 0.0.2 requires
// that the rows read differently, checked rather than the whole schema: a store XRADIO 1.2.5 upgraded
// on read is stamped 0.0.2 whether or not the upgrade could reach them.
std::optional<std::string> WhyNot002(const Store& store) {
    const auto* const groups = MemberObject(store.RootAttributes(), "data_groups");
    if (groups == nullptr) {
        return "it has no data_groups";
    }
    // The base group's sky, which is the one XRADIO writes first; another group's when it has none.
    std::string sky;
    if (const auto* const base = MemberObject(*groups, "base")) {
        sky = AttributeString(*base, "sky");
    }
    for (auto group = groups->begin(); sky.empty() && group != groups->end(); ++group) {
        sky = AttributeString(group.value(), "sky");
    }
    if (sky.empty()) {
        return "no data group names a sky image";
    }
    const auto& image = store.ReadArrayMetadata(sky);
    if (!image) {
        return "its data groups name " + sky + " as the sky, which is not an array";
    }
    if (AttributeString(image.value().attributes, "type") != "sky") {
        return "its sky image " + sky + " is not typed \"sky\"";
    }
    if (const auto& frequency = store.ReadArrayMetadata("frequency"); frequency) {
        const auto& attributes = frequency.value().attributes;
        if (AttributeString(attributes, "frame").empty() || AttributeString(attributes, "units").empty()) {
            return "its frequency coordinate has no frame or no units";
        }
    }
    return std::nullopt;
}

// The row reading stores that name `version`, or none.
const SchemaRow* RowFor(std::string_view version) {
    const auto& rows = SchemaRows();
    const auto found =
        std::find_if(rows.begin(), rows.end(), [&](const SchemaRow& row) { return row.version == version; });
    return found == rows.end() ? nullptr : &*found;
}

}  // namespace

ChosenRow ChooseSchemaRow(const Store& store) {
    ChosenRow chosen;
    const auto* const named = Member(store.RootAttributes(), "schema_version");
    if (named == nullptr) {
        chosen.row = RowFor("");
        chosen.diagnostics.push_back(Diagnostic{DiagnosticCode::unversioned_schema,
                                                "The store names no schema version, so it is read as XRADIO 1.2 "
                                                "wrote image stores",
                                                "/attributes/schema_version"});
        return chosen;
    }
    if (named->is_string()) {
        chosen.version = named->get<std::string>();
    }
    const auto* const row = chosen.version.empty() ? nullptr : RowFor(chosen.version);
    if (row == nullptr) {
        std::string readable;
        for (const auto& candidate : SchemaRows()) {
            if (!candidate.version.empty()) {
                readable += (readable.empty() ? "" : ", ") + std::string(candidate.version);
            }
        }
        const std::string named_as = named->is_string() ? chosen.version : named->dump();
        chosen.diagnostics.push_back(Diagnostic{DiagnosticCode::unsupported_schema_version,
                                                "XRADIO image schema " + named_as +
                                                    " is not one this library reads; it reads " + readable +
                                                    ", and stores that name no version",
                                                "/attributes/schema_version"});
        return chosen;
    }
    if (auto why_not = WhyNot002(store); why_not) {
        chosen.row = RowFor("");
        chosen.diagnostics.push_back(Diagnostic{DiagnosticCode::unversioned_schema,
                                                "The store names schema " + chosen.version + ", but " + *why_not +
                                                    ", so it is read as XRADIO 1.2 wrote image stores",
                                                "/attributes/schema_version"});
        return chosen;
    }
    chosen.row = row;
    return chosen;
}

}  // namespace carta::zarr::internal::xradio
