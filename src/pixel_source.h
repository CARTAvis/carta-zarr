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
 * Shared by the pass and an ordinary read, so it is named after neither: a pass reads slabs and a
 * read reads pieces, and both ask the same two questions.
 *
 * Two adapters: StorePixelSource, which reads them from a `Store`, and the synthetic one in
 * tests/support, which computes them from their own coordinates. The synthetic one lets tests run a
 * walk over a cube of any size, check what it read rather than what it returned, and see what it
 * does at four thousand pixels square, none of which a directory tree on disk can do.
 *
 * It is called once per slab or once per piece, so the indirect call costs once per megabyte of
 * pixels rather than per pixel. That is why it can be an interface while a visitor cannot: see
 * ADR 0005.
 *
 * `destination` runs from where this read's first element lands to the end of the caller's buffer;
 * it is not the selection's size restated. A source writes the selection's elements into the front
 * of it and nothing past them, and refuses a destination too short to hold them before writing any.
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
