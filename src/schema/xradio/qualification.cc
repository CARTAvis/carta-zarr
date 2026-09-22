/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "qualification.h"

#include "../../zarr/array_metadata.h"
#include "flag.h"

#include <algorithm>
#include <cstddef>
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

// A node that says it is a group, asked only once parsing it as an array has already failed. Read
// from the node's own document rather than from the message that parse produced, which says "Zarr
// node is not an array" for a group and for several other things that are not one.
bool IsGroup(const Store& store, std::string_view node) {
    const auto& metadata = store.ReadNodeMetadata(node);
    return metadata && metadata.value().is_object() && metadata.value().value("node_type", "") == "group";
}

// A dataset stores each coordinate once and every image references it by dimension name, so an
// image whose own extent disagrees with the coordinate it names cannot be described with it -- it
// would be reported with another image's coordinate vector, three channels' worth over seven
// channels of pixels.
//
// A coordinate the dataset does not carry at all is the probe's business rather than this one's:
// every image references it, so its absence closes the dataset instead of one variable.
std::optional<Diagnostic> DisagreementWithCoordinates(const Store& store,
                                                      const zarr_metadata::ArrayMetadata& image,
                                                      std::string_view node) {
    const auto rank = std::min(image.dimension_names.size(), image.shape.size());
    for (std::size_t axis = 0; axis < rank; ++axis) {
        const auto& name = image.dimension_names.at(axis);
        const auto& coordinate = store.ReadArrayMetadata(name);
        if (!coordinate) {
            continue;
        }
        if (coordinate.value().shape.size() != 1 || coordinate.value().shape.front() != image.shape.at(axis)) {
            return Diagnostic{"invalid_metadata",
                              "Image dimension '" + name + "' is not the length of the coordinate of that name",
                              std::string(node)};
        }
    }
    return std::nullopt;
}

}  // namespace

NodeQualification QualifyNode(const Store& store, std::string_view node) {
    const auto& metadata = store.ReadArrayMetadata(node);
    if (!metadata) {
        if (IsGroup(store, node)) {
            // A group is a node in the hierarchy like any other, and nothing is wrong with it: it
            // failed to parse as an array because it never claimed to be one. Passed over, the way
            // everything else that is not an image is passed over.
            return NodeQualification{};
        }
        // Said rather than passed over, because this is where a node leaves the listing: nothing
        // here knows whether it was an image, a coordinate or something else in the directory, so a
        // store whose sky variable lands here reports no images at all. Without this, that reads as
        // a store this profile does not recognise, with nothing to say which node went missing.
        //
        // Not a refusal. A malformed array elsewhere in a store does not stop the images that
        // parsed from opening.
        return NodeQualification{false, false,
                                 Diagnostic{"unreadable_array", metadata.error().message, std::string(node)}, true};
    }

    const auto& array = metadata.value();
    if (IsFlag(array)) {
        return NodeQualification{};
    }

    if (HasAllAxes(array, kSkyAxes)) {
        if (!zarr_metadata::IsRealDataType(array.data_type)) {
            return NodeQualification{
                true, false,
                Diagnostic{"unsupported_data_type", "Complex sky-plane variables are not openable", std::string(node)},
                false};
        }
        if (auto disagreement = DisagreementWithCoordinates(store, array, node); disagreement) {
            return NodeQualification{true, false, std::move(disagreement), true};
        }
        return NodeQualification{true, true, std::nullopt, false};
    }

    if (HasAllAxes(array, kApertureAxes)) {
        return NodeQualification{true, false,
                                 Diagnostic{"unsupported_coordinate_plane",
                                            "Aperture-plane variables are not openable", std::string(node)},
                                 false};
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

    const auto qualified = QualifyNode(store, image_id);
    if (qualified.openable) {
        return {};
    }
    if (qualified.listed) {
        // The code comes from the refusal rather than from the asking. A store that is wrong about
        // an image and a library that does not open one are different answers, and a consumer that
        // shows one to whoever picked the file acts on the difference.
        return Error{qualified.malformed ? ErrorCode::invalid_metadata : ErrorCode::unsupported_data_type,
                     qualified.diagnostic->message, std::string(image_id)};
    }
    if (qualified.diagnostic) {
        // The node is there and will not parse. It says what is wrong with it rather than being
        // reported as a name the dataset does not have.
        return store.ReadArrayMetadata(image_id).error();
    }

    // Not an image, and nothing wrong with it: the flag beside one, the coordinate under it, a beam
    // table, a group. A caller that named one asked for an image this dataset does not have.
    return Error{ErrorCode::not_found, "Image variable was not found", std::string(image_id)};
}

}  // namespace carta::zarr::internal::xradio
