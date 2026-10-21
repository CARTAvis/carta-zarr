/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_PIXEL_MASK_H_
#define CARTA_ZARR_SRC_PIXEL_MASK_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace carta::zarr::internal {

// Whether a read of this image folds its pixel mask into the pixels: when the image has one, unless
// the caller declined it.
//
// Image::Read and the pass both ask this, so it is answered once: a read and a reduction that
// disagreed about which pixels exist would be a failure neither reports.
inline bool AppliesPixelMask(const ReadOptions& options, const ImageDescriptor& descriptor) noexcept {
    return options.apply_pixel_mask && !descriptor.pixel_mask_id.empty();
}

/**
 * Drop every pixel its mask excludes, by replacing it with NaN.
 *
 * What a byte of a pixel mask means is the schema profile's decision, as is which flag variable the
 * bytes come from (see DetermineFlag in schema/xradio/flag.h). The rule is written once here rather
 * than at each place that reads pixels. In XRADIO a flag's true means the pixel is flagged (its CASA
 * reader inverts casacore's mask, whose true is a good pixel), so a nonzero byte is a pixel to drop.
 * See ADR 0019.
 *
 * Folding the flag into the pixels lets everything downstream ignore it: a flagged pixel and a NaN
 * pixel mean the same thing to a statistic, to a histogram, and to the consumer's own arithmetic.
 */
inline void ApplyPixelMask(float* pixels, const std::uint8_t* mask, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        if (mask[i] != 0) {
            pixels[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }
}

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_PIXEL_MASK_H_
