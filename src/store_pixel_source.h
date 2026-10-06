/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_STORE_PIXEL_SOURCE_H_
#define CARTA_ZARR_SRC_STORE_PIXEL_SOURCE_H_

// The pass's production adapter, kept apart from the pass itself so that pass.h names no Store.
// Everything under src/reduce/ then compiles without one, which is what lets a reduction be run and
// timed against pixels that were never written down.

#include "pixel_source.h"
#include "store.h"

namespace carta::zarr::internal {

// The adapter that serves production: the image's own data variable, and the flag that masks it.
class StorePixelSource final : public PixelSource {
public:
    StorePixelSource(const Store& store, const ImageDescriptor& descriptor)
        : _store(&store), _descriptor(&descriptor) {}

    Result<void> ReadPixels(const zarr::PixelSelection& selection, BufferView<float> destination,
                            const ReadControl& control) const override {
        return _store->ReadPixelsInto(_descriptor->id, selection, destination, control);
    }

    Result<void> ReadMask(const zarr::PixelSelection& selection, BufferView<std::uint8_t> destination,
                          const ReadControl& control) const override {
        return _store->ReadPixelsInto(_descriptor->pixel_mask_id, selection, destination, control);
    }

private:
    const Store* _store;
    const ImageDescriptor* _descriptor;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_STORE_PIXEL_SOURCE_H_
