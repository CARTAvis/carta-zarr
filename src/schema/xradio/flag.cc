/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "flag.h"

#include "attributes.h"

#include <algorithm>
#include <set>
#include <utility>

namespace carta::zarr::internal::xradio {

namespace {

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
        (AttributeString(group, "sky") == image_id ? declared.own : declared.others).insert(flag);
    }
    return declared;
}

Result<std::string> RequireDeclaredFlag(const Store& store, const zarr::ArrayMetadata& image,
                                        const std::string& declared) {
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

Result<std::string> DetermineFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id,
                                  std::vector<Diagnostic>& diagnostics) {
    // A declared flag is the image's own statement that its pixels need a mask, so a flag that
    // cannot serve as one closes the image rather than opening it unmasked.
    auto declared = DeclaredFlag(store, image, image_id);
    if (!declared || !declared.value().empty()) {
        return declared;
    }
    const auto groups = FlagsOfDataGroups(store.RootAttributes(), image_id);

    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }
    std::vector<std::string> matching_flags;
    for (const auto& entry : inventory.value()) {
        // Nothing declared one, so this is a guess from the metadata alone. A node that is not an
        // array, or whose metadata will not parse, or that does not match exactly, is simply not
        // this image's mask -- which is not an error in the store, and not this module's to report.
        // A flag a data group declares for another image is that image's, however well it fits.
        if (entry.kind != NodeKind::array || !*entry.array || groups.others.count(entry.name) != 0 ||
            !RequireUsableFlag(entry.array->value(), image, entry.name)) {
            continue;
        }
        matching_flags.push_back(entry.name);
    }
    if (matching_flags.size() == 1) {
        return matching_flags.front();
    }
    if (matching_flags.size() > 1) {
        diagnostics.push_back(Diagnostic{DiagnosticCode::ambiguous_pixel_mask,
                                         "More than one flag variable matches the image; no pixel mask was selected",
                                         std::string(image_id)});
    }
    return std::string{};
}

}  // namespace carta::zarr::internal::xradio
