/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_STORE_SLAB_SOURCE_H_
#define CARTA_ZARR_SRC_REDUCE_STORE_SLAB_SOURCE_H_

// The pass's production adapter, kept apart from the pass itself so that pass.h names no Store.
// Everything under src/reduce/ then compiles without one, which is what lets a reduction be run and
// timed against pixels that were never written down.

#include "reduce/pass.h"
#include "store.h"

namespace carta::zarr::internal {

// The adapter that serves production: the image's own data variable, and the flag that masks it.
class StorePixelSource final : public PixelSource {
public:
    StorePixelSource(const Store& store, const ImageDescriptor& descriptor)
        : _store(&store), _descriptor(&descriptor) {}

    Result<void> ReadPixels(const zarr::PixelSelection& selection, float* destination, std::size_t elements,
                            const ReadOptions& options) const override {
        return _store->ReadPixelsInto(_descriptor->id, selection, destination, elements, options);
    }

    Result<void> ReadMask(const zarr::PixelSelection& selection, std::uint8_t* destination,
                          std::size_t elements, const ReadOptions& options) const override {
        return _store->ReadPixelsInto(_descriptor->pixel_mask_id, selection, destination, elements, options);
    }

private:
    const Store* _store;
    const ImageDescriptor* _descriptor;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_STORE_SLAB_SOURCE_H_
