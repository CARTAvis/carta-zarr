/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_PIXEL_SOURCE_H_
#define CARTA_ZARR_SRC_PIXEL_SOURCE_H_

#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include "zarr/pixel_selection.h"

#include <cstdint>

namespace carta::zarr::internal {

/**
 * Where an image's pixels and its flag come from.
 *
 * The seam neither the pass nor an ordinary read owns, which is why it is named after neither of
 * their units: a pass reads slabs and a read reads pieces, and both ask the same two questions.
 *
 * Two adapters: StorePixelSource, which reads them from a `Store`, and the synthetic one in
 * tests/support, which computes them from their own coordinates. The second is what lets a walk be
 * asked about a cube of any size, about what it read rather than what it handed over, and about
 * what it does at four thousand pixels square -- none of which a directory tree can answer.
 *
 * Consulted once per slab or once per piece, so the indirect call is paid per megabyte of pixels
 * rather than per pixel. That is the whole reason it can be an interface at all while a visitor
 * cannot: see ADR 0005.
 *
 * `destination` is what the caller has to write into, from where this read's first element lands
 * to the end of the caller's buffer -- not the selection's size restated. A source writes the
 * selection's elements into the front of it and nothing past them, and refuses a destination too
 * short to hold them before writing any. The length used to cross as a separate count that every
 * caller worked out from the selection it was passing, so the check against it compared the
 * selection with itself and could never catch the caller's mistake it was there for.
 */
class PixelSource {
public:
    PixelSource() = default;
    PixelSource(const PixelSource&) = delete;
    PixelSource& operator=(const PixelSource&) = delete;
    PixelSource(PixelSource&&) = delete;
    PixelSource& operator=(PixelSource&&) = delete;
    virtual ~PixelSource() = default;

    virtual Result<void> ReadPixels(const zarr::PixelSelection& selection, BufferView<float> destination,
                                    const ReadControl& control) const = 0;
    // Only called when the caller has established that the image has a flag and that this read
    // applies it. A source that has none may report not_found rather than serve zeroes.
    virtual Result<void> ReadMask(const zarr::PixelSelection& selection, BufferView<std::uint8_t> destination,
                                  const ReadControl& control) const = 0;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_PIXEL_SOURCE_H_
