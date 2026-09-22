/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "qualification.h"

#include "flag.h"

#include <algorithm>
#include <string>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// The Fourier-conjugate plane. Variables on it are images in every other respect, which is why they
// are listed with a reason rather than passed over in silence.
constexpr std::array<std::string_view, 5> kApertureAxes{"time", "frequency", "polarization", "u", "v"};

bool HasAllAxes(const zarr_metadata::ArrayMetadata& metadata, const std::array<std::string_view, 5>& axes) {
    return std::all_of(axes.begin(), axes.end(), [&metadata](const auto axis) {
        return zarr_metadata::FindDimensionIndex(metadata, axis).has_value();
    });
}

}  // namespace

NodeQualification QualifyNode(const Result<zarr::ArrayMetadata>& metadata, std::string_view node) {
    if (!metadata) {
        // Said rather than passed over, because this is where a node leaves the listing: nothing
        // here knows whether it was an image, a coordinate or something else in the directory, so a
        // store whose sky variable lands here reports no images at all. Without this, that reads as
        // a store this profile does not recognise, with nothing to say which node went missing.
        //
        // Not a refusal. A malformed array elsewhere in a store does not stop the images that
        // parsed from opening.
        return NodeQualification{false, false, Diagnostic{"unreadable_array", metadata.error().message, std::string(node)}};
    }

    const auto& array = metadata.value();
    if (IsFlag(array)) {
        return NodeQualification{};
    }

    if (HasAllAxes(array, kSkyAxes)) {
        if (zarr_metadata::IsRealDataType(array.data_type)) {
            return NodeQualification{true, true, std::nullopt};
        }
        return NodeQualification{
            true, false,
            Diagnostic{"unsupported_data_type", "Complex sky-plane variables are not openable", std::string(node)}};
    }

    if (HasAllAxes(array, kApertureAxes)) {
        return NodeQualification{true, false,
                                 Diagnostic{"unsupported_coordinate_plane",
                                            "Aperture-plane variables are not openable", std::string(node)}};
    }

    // A coordinate array, a normalization variable, a beam table: not an image, and nothing to say.
    return NodeQualification{};
}

Result<void> RequireQualified(const Store& store, std::string_view image_id) {
    // Against the store's listing rather than against whatever is on disk under that name. A Store
    // is a read-only view, and its listing is the snapshot a dataset's images were enumerated from;
    // a node written beside them afterwards is not one of them, and reading it by name would reach
    // it. What a listing offers and what will open are the same set in both directions.
    const auto& nodes = store.ListNodes();
    if (!nodes) {
        return nodes.error();
    }
    if (std::find(nodes.value().begin(), nodes.value().end(), image_id) == nodes.value().end()) {
        return Error{ErrorCode::not_found, "Image variable was not found", std::string(image_id)};
    }

    const auto& metadata = store.ReadArrayMetadata(image_id);
    const auto qualified = QualifyNode(metadata, image_id);
    if (qualified.openable) {
        return {};
    }
    if (qualified.listed) {
        return Error{ErrorCode::unsupported_data_type, qualified.diagnostic->message, std::string(image_id)};
    }

    // A node that is there and will not parse says so, rather than being reported as a name the
    // dataset does not have. Everything else here is a node that is not an image at all -- the flag
    // beside one, the coordinate under it, a beam table -- and a caller that named one asked for an
    // image this dataset does not have.
    if (!metadata) {
        return metadata.error();
    }
    return Error{ErrorCode::not_found, "Image variable was not found", std::string(image_id)};
}

}  // namespace carta::zarr::internal::xradio
