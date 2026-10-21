/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_

// Which variable supplies an image's pixel mask, and what a variable has to be before it can.
//
// The flag is the variable on disk; the pixel mask is what it becomes once a read folds it into the
// pixels. That rule lives in src/pixel_mask.h. This module only chooses the variable.

#include "carta-zarr/descriptor.h"
#include "carta-zarr/error.h"
#include "carta-zarr/result.h"

#include "../../store.h"
#include "../../zarr/array_metadata.h"

#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

// A flag is distinguished by its `type` attribute rather than by its name: the XRADIO CASA reader
// names variables after CASA's internal masks, so MASK_0 is a flag and MASK, a deconvolution's mask,
// is an ordinary image. Classification must never be driven by variable names. See ADR 0001.
inline bool IsFlag(const zarr::ArrayMetadata& metadata) {
    return metadata.attributes.contains("type") && metadata.attributes.at("type").is_string() &&
           metadata.attributes.at("type").get<std::string>() == "flag";
}

// Whether a variable's values are booleans: a Zarr bool, which is how XRADIO writes every flag. An
// int8, whatever its attributes say, is a number.
inline bool HoldsBooleans(const zarr::ArrayMetadata& metadata) {
    return metadata.data_type == "bool";
}

/**
 * What a variable has to be before it can mask this image's pixels.
 *
 * The read applies the mask element by element against the same selection it read pixels with, so
 * anything less than the image's own dimensions in the image's own order is not a mask of it.
 *
 * Inline, and part of this module's interface, because it is everything a declared flag is checked
 * against: the three reasons to refuse a mask need only two pieces of metadata, with no store.
 */
inline Result<void> RequireUsableFlag(const zarr::ArrayMetadata& flag, const zarr::ArrayMetadata& image,
                                      std::string_view node) {
    if (!IsFlag(flag)) {
        return Error{ErrorCode::invalid_metadata, "Pixel mask variable is not marked as a flag", std::string(node)};
    }
    if (!HoldsBooleans(flag)) {
        return Error{ErrorCode::unsupported_data_type, "Pixel mask variable is not boolean", std::string(node)};
    }
    if (flag.dimension_names != image.dimension_names || flag.shape != image.shape) {
        return Error{ErrorCode::invalid_metadata, "Pixel mask variable does not share the image's dimensions",
                     std::string(node)};
    }
    return {};
}

/**
 * The flag this image declares, checked as able to mask it, or an empty name when it declares none.
 *
 * Declared either by the image's own `flag` attribute or by a root `data_groups` entry whose `sky` is
 * this image; the attribute outranks a group. A declaration that cannot be used (a flag that is
 * missing or unusable, or two groups naming different ones) is an error. DetermineFlag and the
 * listing's qualification both call this, so an image is never listed as openable when its
 * declaration would refuse it on opening.
 */
Result<std::string> DeclaredFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id);

/**
 * The flag variable supplying this image's pixel mask, or an empty name when it has none.
 *
 * A flag is declared either by the image's own `flag` attribute or by a root `data_groups` entry
 * whose `sky` is this image, which is the only place XRADIO's writer records it. A declared flag is
 * binding: an image naming one that cannot be used is closed rather than opened unmasked. Reads apply
 * the mask by default, and treating an unusable mask as no mask would show flagged pixels as valid,
 * which a consumer has no way to notice. With nothing declared, a sky image has no flag, and any
 * other image's flag is the array XRADIO names FLAG_<image>, if it can mask the image. No flag is
 * found by its shape alone. See ADR 0019.
 *
 * Takes the whole Store rather than a narrower interface: what varies underneath is the transport,
 * which already has two adapters, and a seam here would have only one. See ADR 0006.
 */
Result<std::string> DetermineFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_
