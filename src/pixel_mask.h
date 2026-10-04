/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
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
// Stated here, beside what folding it in means, for the reason that is stated once: Image::Read and
// the pass both ask it, and a read that disagreed with a reduction about which pixels exist is the
// one failure neither reports. They asked it in two places, each writing out the same conjunction.
inline bool AppliesPixelMask(const ReadOptions& options, const ImageDescriptor& descriptor) noexcept {
    return options.apply_pixel_mask && descriptor.has_pixel_mask;
}

/**
 * Drop every pixel its mask excludes, by replacing it with NaN.
 *
 * What a byte of a pixel mask means belongs to the schema profile, which is also what chose the
 * flag variable the bytes came from -- see DetermineFlag in schema/xradio/flag.h. The rule
 * is stated once here rather than at each place that reads pixels, so that a profile disagreeing
 * about it has one line to change instead of a search to run. XRADIO writes a flag whose true means
 * the pixel is flagged -- its CASA reader inverts casacore's mask, whose true is a good pixel -- so a
 * nonzero byte is the pixel to drop. See ADR 0019.
 *
 * Folding the flag into the pixels is what lets everything downstream ignore it: a flagged pixel
 * and a NaN pixel mean the same thing to a statistic, to a histogram, and to the consumer's own
 * arithmetic. Image::Read and the pass both do this, and a read that disagreed with a reduction
 * about which pixels exist is the one failure neither reports.
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
