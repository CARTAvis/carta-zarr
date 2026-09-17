/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_

// Which variable supplies an image's pixel mask, and what a variable has to be before it can.
//
// The flag is the variable on disk; the pixel mask is what it becomes once a read folds it into the
// pixels. That rule lives in src/pixel_mask.h and is a different question from this one -- this
// module only chooses the variable.

#include "../../store.h"
#include "../../zarr/array_metadata.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/error.h"
#include "carta-zarr/result.h"

#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

// A flag is distinguished by its `type` attribute rather than by its name: the XRADIO CASA reader
// names variables after CASA's internal masks, so MASK_0 is a flag and MASK_DECONVOLVE is an
// ordinary image. Classification must never be driven by variable names. See ADR 0001.
inline bool IsFlag(const zarr::ArrayMetadata& metadata) {
    return metadata.attributes.contains("type") && metadata.attributes.at("type").is_string() &&
           metadata.attributes.at("type").get<std::string>() == "flag";
}

/**
 * What a variable has to be before it can mask this image's pixels.
 *
 * The read applies the mask element by element against the same selection it read pixels with, so
 * anything less than the image's own dimensions in the image's own order is not a mask of it.
 *
 * Inline, and part of this module's interface rather than private to it, because it is the whole of
 * what a declared flag is checked against: the three ways a mask can be refused are answerable from
 * two pieces of metadata, with no store to stand up.
 */
inline Result<void> RequireUsableFlag(const zarr::ArrayMetadata& flag, const zarr::ArrayMetadata& image,
                                      std::string_view node) {
    if (!IsFlag(flag)) {
        return Error{ErrorCode::invalid_metadata, "Pixel mask variable is not marked as a flag", std::string(node)};
    }
    if (flag.data_type != "bool") {
        return Error{ErrorCode::unsupported_data_type, "Pixel mask variable is not boolean", std::string(node)};
    }
    if (flag.dimension_names != image.dimension_names || flag.shape != image.shape) {
        return Error{ErrorCode::invalid_metadata, "Pixel mask variable does not share the image's dimensions",
                     std::string(node)};
    }
    return {};
}

/**
 * The flag variable supplying this image's pixel mask, or an empty name when it has none.
 *
 * A declared flag is binding: an image naming one that cannot serve is closed rather than opened
 * unmasked, because reads apply the mask by default and an unusable mask reported as no mask would
 * show flagged pixels as valid -- the one failure a consumer has no way to notice. With nothing
 * declared the store is inspected instead, and more than one match is refused rather than guessed
 * at, as a diagnostic on the image.
 *
 * Takes the whole Store rather than a narrower interface: what varies underneath is the transport,
 * which already has two adapters, and a seam here would have exactly one. See ADR 0006.
 */
Result<std::string> DetermineFlag(const Store& store, const zarr::ArrayMetadata& image, std::string_view image_id,
                                  std::vector<Diagnostic>& diagnostics);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_FLAG_H_
