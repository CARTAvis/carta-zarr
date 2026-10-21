/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "flag.h"

#include "attributes.h"

#include <algorithm>
#include <set>
#include <utility>

namespace carta::zarr::internal::xradio {

namespace {

// A node a document names, as the store files it, so that it can be compared with the names the
// store lists: "./FLAG" is FLAG. A name that is no node is kept as written, and reading it reports so.
std::string Canonical(std::string_view node) {
    auto normalized = NormalizeNodeName(node);
    return normalized ? std::move(normalized).value() : std::string(node);
}

// The flags the root's data groups declare: those of groups whose `sky` is this image, and those
// of every other group, which belong to some other image. XRADIO's schema defines a group's `flag`
// as its sky image's, so nothing else in a group is taken to own it.
struct DeclaredFlags {
    std::set<std::string> own;
    std::set<std::string> others;
};

DeclaredFlags FlagsOfDataGroups(const nlohmann::json& root_attributes, std::string_view image_id) {
    DeclaredFlags declared;
    const auto* const groups = MemberObject(root_attributes, "data_groups");
    if (groups == nullptr) {
        return declared;
    }
    for (const auto& [name, group] : groups->items()) {
        if (!group.is_object()) {
            continue;
        }
        const auto flag = AttributeString(group, "flag");
        if (flag.empty()) {
            continue;
        }
        // Canonical on both sides: compared as written, another group's "./FLAG" is not the FLAG the
        // inventory lists, and an image would take a flag another image owns.
        (Canonical(AttributeString(group, "sky")) == image_id ? declared.own : declared.others).insert(Canonical(flag));
    }
    return declared;
}

Result<std::string> RequireDeclaredFlag(const Store& store, const zarr::ArrayMetadata& image,
                                        const std::string& spelt) {
    const std::string declared = Canonical(spelt);
    const auto& flag_array = store.ReadArrayMetadata(declared);
    if (!flag_array) {
        return flag_array.error();
    }
    if (auto usable = RequireUsableFlag(flag_array.value(), image, declared); !usable) {
        return usable.error();
    }
    return declared;
}

}  // namespace

Result<std::string> DeclaredFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id) {
    if (auto declared = AttributeString(image.attributes, "flag"); !declared.empty()) {
        return RequireDeclaredFlag(store, image, declared);
    }
    // XRADIO's writer records a flag here and not on the image. Two groups sharing this sky image
    // may each name its flag; naming two different ones is a contradiction, not a choice to make.
    const auto groups = FlagsOfDataGroups(store.RootAttributes(), image_id);
    if (groups.own.size() > 1) {
        return Error{ErrorCode::invalid_metadata, "data_groups name more than one flag for the image",
                     std::string(image_id)};
    }
    if (groups.own.size() == 1) {
        return RequireDeclaredFlag(store, image, *groups.own.begin());
    }
    return std::string{};
}

Result<std::string> DetermineFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id) {
    // A declared flag is the image's own statement that its pixels need a mask, so a flag that
    // cannot serve as one closes the image rather than opening it unmasked.
    auto declared = DeclaredFlag(store, image, image_id);
    if (!declared || !declared.value().empty()) {
        return declared;
    }
    // A sky image's flag is declared wherever XRADIO writes one, so a sky image that declares none
    // has none. Another image's cannot be declared by a data group -- a group's `flag` is its sky's --
    // and XRADIO names it FLAG_<image>, which is the one array taken, and only if it can mask the
    // image. A flag a data group declares for a sky image is that image's, whatever its name.
    if (AttributeString(image.attributes, "type") == "sky") {
        return std::string{};
    }
    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }
    const std::string named = "FLAG_" + std::string(image_id);
    const auto* const entry = store.FindNode(named);
    if (entry == nullptr || entry->kind != NodeKind::array || !*entry->array ||
        FlagsOfDataGroups(store.RootAttributes(), image_id).others.count(named) != 0 ||
        !RequireUsableFlag(entry->array->value(), image, named)) {
        return std::string{};
    }
    return named;
}

}  // namespace carta::zarr::internal::xradio
