/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "qualification.h"

#include "../../zarr/array_metadata.h"
#include "attributes.h"
#include "coordinates.h"
#include "flag.h"

#include <cstddef>
#include <string>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// Why the flag an image declares cannot mask it, or nothing when it declares none or one that can.
// The listing asks DeclaredFlag, as DetermineFlag does, so that the two agree. With nothing declared
// there is nothing to refuse, because an image with no usable candidate opens unmasked.
std::optional<Diagnostic> UnusableDeclaredFlag(const Store& store, const zarr_metadata::ArrayMetadata& image,
                                               std::string_view node) {
    const auto declared = DeclaredFlag(store, image, node);
    if (declared) {
        return std::nullopt;
    }
    return Diagnostic{
        DiagnosticCode::invalid_metadata,
        "Declared flag '" + declared.error().node_path + "' cannot mask this image: " + declared.error().message,
        std::string(node)};
}

// Why a data group's role for this node and the node's own type disagree, or nothing when they do
// not. A group naming a variable as its `sky`, or in another image role, says what it is, and the
// schema requires the variable to have that role's type: a SKY typed "Intensity" is the store
// contradicting itself, not an image of a type this library does not know.
std::optional<Diagnostic> RoleContradiction(const Store& store, std::string_view node, const std::string& type) {
    const auto* const groups = MemberObject(store.RootAttributes(), "data_groups");
    if (groups == nullptr) {
        return std::nullopt;
    }
    for (const auto& [group_name, group] : groups->items()) {
        if (!group.is_object()) {
            continue;
        }
        for (const auto& [role, value] : group.items()) {
            if (!ImagePlane(role) || !value.is_string() || role == type) {
                continue;
            }
            const auto named = NormalizeNodeName(value.get<std::string>());
            if (named && named.value() == node) {
                return Diagnostic{
                    DiagnosticCode::invalid_metadata,
                    "Data group '" + group_name + "' names it as its " + role + ", but it is typed '" + type + "'",
                    std::string(node)};
            }
        }
    }
    return std::nullopt;
}

}  // namespace

NodeQualification QualifyNode(const Store& store, const NodeEntry& entry) {
    const std::string& node = entry.name;
    // A group is a node like any other, and nothing is wrong with it. It is passed over, as everything
    // else that is not an image is.
    if (entry.kind == NodeKind::group) {
        return NodeQualification{};
    }

    // Diagnosed rather than passed over, because this is where a node leaves the listing: nothing here
    // knows whether it was an image, a coordinate or something else, so a store whose sky variable
    // lands here reports no images at all. Without a diagnostic, that would look like a store this
    // profile does not recognize, with nothing to say which node went missing.
    //
    // Not a refusal. A malformed node elsewhere in a store does not stop the images that parsed from
    // opening.
    //
    // Two codes, for two different problems. unreadable_array is a node that said it was an array and
    // whose metadata will not parse as one; unrecognized_node is a node that did not say what it was
    // (a document that would not parse, or a node_type Zarr does not define), and calling that an
    // array would be a guess.
    if (entry.kind == NodeKind::unrecognized) {
        return NodeQualification{false, false,
                                 Diagnostic{DiagnosticCode::unrecognized_node, entry.reason->message, node}, true};
    }
    const auto& metadata = *entry.array;
    if (!metadata) {
        return NodeQualification{false, false,
                                 Diagnostic{DiagnosticCode::unreadable_array, metadata.error().message, node}, true};
    }

    const auto& array = metadata.value();
    // A variable is an image by its type: schema 0.0.2 types every image. An untyped variable, or one
    // of a type that is not an image type (a flag, a coordinate, a normalization, a beam table), is
    // not an image, whatever its dimensions (ADR 0001). The data groups are checked first, because a
    // flag named as a group's sky contradicts the group as much as any other type does.
    const std::string type = AttributeString(array.attributes, "type");
    if (auto contradiction = RoleContradiction(store, node, type); contradiction) {
        return NodeQualification{true, false, std::move(contradiction), true};
    }
    if (IsFlag(array)) {
        return NodeQualification{};
    }
    const auto plane = ImagePlane(type);
    if (plane) {
        const bool sky = *plane == Plane::sky;
        const std::string on_plane = sky ? "sky-plane" : "aperture-plane";
        // Typed as an image of a plane, it has to be over that plane's coordinates, or it would be
        // described by axes it does not have.
        if (!CarriesPlane(array, *plane)) {
            return NodeQualification{true, false,
                                     Diagnostic{DiagnosticCode::invalid_metadata,
                                                "Typed '" + type + "', an image of the " + (sky ? "sky" : "aperture") +
                                                    " plane, but its dimensions are not time, frequency, "
                                                    "polarization, " +
                                                    (sky ? "l and m" : "u and v"),
                                                std::string(node)},
                                     true};
        }
        // This library does not open an image of booleans: booleans over an image's axes are a flag's
        // values, and the message says a flag is typed 'flag'.
        if (HoldsBooleans(array)) {
            return NodeQualification{
                true, false,
                Diagnostic{DiagnosticCode::unsupported_data_type,
                           "Boolean " + on_plane + " variables are not openable; a flag is typed 'flag'",
                           std::string(node)},
                false};
        }
        // On either plane, this library does not open a complex image, such as a visibility grid
        // (ADR 0024).
        if (!zarr_metadata::IsRealDataType(array.data_type)) {
            return NodeQualification{
                true, false,
                Diagnostic{DiagnosticCode::unsupported_data_type, "Complex " + on_plane + " variables are not openable",
                           std::string(node)},
                false};
        }
        // The five axes are all an image is described by, so a sixth would be read with a selection one
        // rank short, and every read of an image listed as openable would be refused.
        if (array.dimension_names.size() != AxisCount(*plane)) {
            return NodeQualification{true, false,
                                     Diagnostic{DiagnosticCode::invalid_metadata,
                                                std::string(sky ? "Sky" : "Aperture") +
                                                    "-plane variable has dimensions beyond time, frequency, "
                                                    "polarization, " +
                                                    (sky ? "l and m" : "u and v"),
                                                std::string(node)},
                                     true};
        }
        // The sky plane's coordinates belong to the dataset: the probe refuses a dataset missing one.
        // The aperture plane's belong to its images alone, and an aperture-plane variable written
        // without them (by hand, not by XRADIO) cannot be described, so it does not open, rather than
        // closing the sky images beside it.
        if (!sky) {
            for (const auto& coordinate : kCoordinates) {
                if (coordinate.plane != Plane::aperture) {
                    continue;
                }
                const auto& metadata = store.ReadArrayMetadata(coordinate.name);
                const auto checked =
                    metadata ? CheckCoordinate(metadata.value(), coordinate) : Result<void>{metadata.error()};
                if (!checked) {
                    return NodeQualification{
                        true, false,
                        Diagnostic{DiagnosticCode::invalid_metadata,
                                   "Aperture-plane variable's coordinate '" + std::string(coordinate.name) +
                                       "' is missing or malformed: " + checked.error().message,
                                   std::string(node)},
                        true};
                }
            }
        }
        if (auto disagreement = ExtentDisagreement(store, array, node); disagreement) {
            return NodeQualification{true, false, std::move(disagreement), true};
        }
        // A declared flag binds the image to it (DetermineFlag closes an image whose flag cannot mask
        // it), so the listing has to say so too, or it would offer an image that will not open.
        if (auto unusable = UnusableDeclaredFlag(store, array, node); unusable) {
            return NodeQualification{true, false, std::move(unusable), true};
        }
        return NodeQualification{true, true, std::nullopt, false};
    }

    // A coordinate array, a normalization variable, a beam table: not an image, and nothing to say.
    return NodeQualification{};
}

Result<void> RequireQualified(const Store& store, std::string_view image_id) {
    // Checked against the store's listing rather than whatever is on disk under that name. A Store is
    // a read-only view, and its listing is the snapshot a dataset's images were enumerated from; a
    // node written beside them afterwards is not one of them, but reading it by name would reach it.
    // What a listing offers and what will open are exactly the same set.
    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }
    const auto* entry = store.FindNode(image_id);
    if (entry == nullptr) {
        return Error{ErrorCode::not_found, "Image variable was not found", std::string(image_id)};
    }

    const auto qualified = QualifyNode(store, *entry);
    if (qualified.openable) {
        return {};
    }
    if (qualified.listed) {
        // The code comes from the reason for the refusal. A store that is wrong about an image and a
        // library that does not open one are different answers, and a consumer showing one to
        // whoever picked the file acts on the difference.
        return Error{qualified.malformed ? ErrorCode::invalid_metadata : ErrorCode::unsupported_data_type,
                     qualified.diagnostic->message, std::string(image_id)};
    }
    if (qualified.diagnostic) {
        // The node is there and will not parse. Report what is wrong with it rather than a name the
        // dataset does not have.
        return entry->kind == NodeKind::array ? entry->array->error() : *entry->reason;
    }

    // Not an image, and nothing wrong with it: a flag, a coordinate, a beam table, a group. A caller
    // that named one asked for an image this dataset does not have.
    return Error{ErrorCode::not_found, "Image variable was not found", std::string(image_id)};
}

}  // namespace carta::zarr::internal::xradio
