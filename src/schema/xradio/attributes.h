/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_

// Reading a value out of an XRADIO attribute object, without deciding what it means.
//
// Shared rather than file-local because the things that interpret those attributes are worth
// compiling on their own, and each of them needs the same four accessors.

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

inline bool HasAttribute(const nlohmann::json& attributes, std::string_view name) {
    return attributes.is_object() && attributes.contains(name);
}

inline std::string AttributeString(const nlohmann::json& attributes, std::string_view name) {
    if (attributes.is_object() && attributes.contains(name) && attributes.at(name).is_string()) {
        return attributes.at(name).get<std::string>();
    }
    return {};
}

// Takes the value rather than an object and a name, because what it reads is usually already in
// hand: an XRADIO measure keeps its number under `data`, which the caller has reached for anyway.
inline std::optional<double> AsNumber(const nlohmann::json& value) {
    return value.is_number() ? std::optional<double>(value.get<double>()) : std::nullopt;
}

inline const nlohmann::json* ObjectMember(const nlohmann::json& object, std::string_view name) {
    if (!object.is_object() || !object.contains(name)) {
        return nullptr;
    }
    return &object.at(name);
}

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_ATTRIBUTES_H_
