/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_

// Reading a value out of an XRADIO attribute object, without deciding what it means.
//
// Shared rather than file-local because the modules that interpret those attributes are compiled on
// their own, and each needs the same accessors.
//
// Two families, which differ in what absence means. `Attribute...` returns a value, and an empty one
// when the attribute is absent, which suits a descriptor field that defaults to empty. `Member...`
// returns a pointer or an optional, for a caller that needs to know whether the attribute is there
// before going further in: an XRADIO measure nests its number under `data` and its units under
// `attrs`, so most callers walk downwards.
//
// `Member` returns whatever is under the name; `MemberObject`, `MemberArray` and `MemberNumber`
// return nothing unless it is of that kind. Checking the kind here saves a separate
// `&& member->is_object()` at every call site, where forgetting it would throw out of nlohmann.

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <optional>
#include <string>
#include <string_view>

namespace carta::zarr::internal::xradio {

inline std::string Upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

inline std::string AttributeString(const nlohmann::json& attributes, std::string_view name) {
    if (attributes.is_object() && attributes.contains(name) && attributes.at(name).is_string()) {
        return attributes.at(name).get<std::string>();
    }
    return {};
}

inline const nlohmann::json* Member(const nlohmann::json& object, std::string_view name) {
    if (!object.is_object() || !object.contains(name)) {
        return nullptr;
    }
    return &object.at(name);
}

inline const nlohmann::json* MemberObject(const nlohmann::json& object, std::string_view name) {
    const auto* const member = Member(object, name);
    return member != nullptr && member->is_object() ? member : nullptr;
}

inline const nlohmann::json* MemberArray(const nlohmann::json& object, std::string_view name) {
    const auto* const member = Member(object, name);
    return member != nullptr && member->is_array() ? member : nullptr;
}

inline std::optional<double> MemberNumber(const nlohmann::json& object, std::string_view name) {
    const auto* const member = Member(object, name);
    return member != nullptr && member->is_number() ? std::optional<double>(member->get<double>()) : std::nullopt;
}

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_
