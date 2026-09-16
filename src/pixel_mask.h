/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_PIXEL_MASK_H_
#define CARTA_ZARR_SRC_PIXEL_MASK_H_

#include <cstddef>
#include <cstdint>
#include <limits>

namespace carta::zarr::internal {

/**
 * Drop every pixel its mask excludes, by replacing it with NaN.
 *
 * What a byte of a pixel mask means belongs to the schema profile, which is also what chose the
 * flag variable the bytes came from -- see DeterminePixelMask in schema/xradio/image.cc. The rule
 * is stated once here rather than at each place that reads pixels, so that a profile disagreeing
 * about it has one line to change instead of a search to run. XRADIO writes a flag whose true means
 * the pixel is good, so a zero byte is the pixel to drop.
 *
 * Folding the flag into the pixels is what lets everything downstream ignore it: a flagged pixel
 * and a NaN pixel mean the same thing to a statistic, to a histogram, and to the consumer's own
 * arithmetic. Image::Read and the pass both do this, and a read that disagreed with a reduction
 * about which pixels exist is the one failure neither reports.
 */
inline void ApplyPixelMask(float* pixels, const std::uint8_t* mask, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        if (mask[i] == 0) {
            pixels[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }
}

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_PIXEL_MASK_H_
