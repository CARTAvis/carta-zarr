/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "flag.h"

#include "attributes.h"

#include <utility>

namespace carta::zarr::internal::xradio {

Result<std::string> DetermineFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id,
                                  std::vector<Diagnostic>& diagnostics) {
    // A declared flag is the image's own statement that its pixels need a mask, so a flag that
    // cannot serve as one closes the image rather than opening it unmasked.
    if (auto declared = AttributeString(image.attributes, "flag"); !declared.empty()) {
        const auto& flag_array = store.ReadArrayMetadata(declared);
        if (!flag_array) {
            return flag_array.error();
        }
        if (auto usable = RequireUsableFlag(flag_array.value(), image, declared); !usable) {
            return usable.error();
        }
        return declared;
    }

    const auto& nodes = store.ListNodes();
    if (!nodes) {
        return nodes.error();
    }
    std::vector<std::string> matching_flags;
    for (const auto& node : nodes.value()) {
        const auto& flag_array = store.ReadArrayMetadata(node);
        // Nothing declared one, so this is a guess from the metadata alone. A variable that does not
        // match exactly is simply not this image's mask, which is not an error in the store.
        if (!flag_array || !RequireUsableFlag(flag_array.value(), image, node)) {
            continue;
        }
        matching_flags.push_back(node);
    }
    if (matching_flags.size() == 1) {
        return matching_flags.front();
    }
    if (matching_flags.size() > 1) {
        diagnostics.push_back(Diagnostic{"ambiguous_pixel_mask",
                                         "More than one flag variable matches the image; no pixel mask was selected",
                                         std::string(image_id)});
    }
    return std::string{};
}

}  // namespace carta::zarr::internal::xradio
