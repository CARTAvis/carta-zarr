/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_VERSION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_VERSION_H_

// The one version of XRADIO's image schema this library reads, and whether a store is one (ADR 0023).

#include "carta-zarr/descriptor.h"

#include "../../store.h"

#include <array>
#include <optional>
#include <string>
#include <string_view>

namespace carta::zarr::internal::xradio {

inline constexpr std::string_view kSchemaVersion = "0.0.2";

// XRADIO's names for a dataset's images. Images that are no data group's own are listed in this
// order, and any other name after them, by name.
inline constexpr std::array<std::string_view, 7> kKnownImages{
    "SKY", "SKY_MODEL", "SKY_RESIDUAL", "SKY_DIRTY", "POINT_SPREAD_FUNCTION", "PRIMARY_BEAM", "MASK"};

// What a store's schema_version says, and whether the store is read.
struct SchemaVersionCheck {
    // As the store names it, which is what a probe reports; empty when it names none.
    std::string version;
    // Why the store is not read, or nothing when it is.
    std::optional<Diagnostic> refusal;
};

// A store is read when it names 0.0.2 and conforms to it: it has data groups, the sky image they name
// (if any) is an array, and its frequency coordinate (if any) gives a frame and units. Whether that
// sky is typed "sky" is checked for the image by qualification. A store naming no version or another
// version is refused with unsupported_schema_version, and one naming 0.0.2 that does not conform with
// invalid_metadata.
SchemaVersionCheck CheckSchemaVersion(const Store& store);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_VERSION_H_
