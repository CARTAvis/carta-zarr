/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_
#define CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_

#include "axis_map.h"
#include "pixel_source.h"
#include "work_pool.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

namespace carta::zarr::internal {

/**
 * One image, opened, with everything a reduction needs before a request arrives.
 *
 * Four things travelled together to every reduction -- where the pixels come from, what the image
 * is, how it is laid out, and what may run the arithmetic -- and the descriptor travelled twice,
 * once on its own and once inside the source. Three entry points took six or seven arguments of
 * which six were the same six in the same order.
 *
 * It is not only the bundle, because a bundle would just move those arguments rather than absorb
 * them: it also answers what each reduction used to work out again for itself. Where the roles sit
 * among the axes is derived once here rather than three times.
 *
 * An ordinary read is not a reduction and does not come through here. It did, while this was called
 * ReadableImage, and so it was handed an axis map and a pool it never asked anything of, and a read
 * of an image whose axes could not be mapped would have failed for a reason that belonged to the
 * walks. No image that opens today is such an image, which is why that was never seen. ReadInPieces
 * takes the three things it uses; see ADR 0005 for why a read stays outside the pass.
 *
 * Checking a plane selection used to live here too, as `ValidateSpectral`. That was one third of
 * the question in the one place the other two thirds were not; `CheckedPlanes::Of` asks all of it,
 * from the descriptor and the map this already holds.
 *
 * Holds references and an AxisMap by value, so it is cheap to build per call and owns nothing. It
 * must not outlive the source, the descriptor, the geometry or the pool it was built from.
 */
class ReducibleImage {
public:
    // Reports not_implemented when the image's axes cannot be mapped -- an axis with no known role
    // that is not degenerate, or a missing spatial or spectral axis. Built per call, so that
    // failure reaches the caller of the reduction that needed it and no other.
    static Result<ReducibleImage> Of(const PixelSource& source, const ImageDescriptor& descriptor,
                                    const ChunkGeometry& geometry, WorkPool& workers) {
        auto map = MapAxes(descriptor);
        if (!map) {
            return map.error();
        }
        return ReducibleImage(source, descriptor, geometry, workers, map.value());
    }

    const PixelSource& source() const noexcept {
        return *_source;
    }
    const ImageDescriptor& descriptor() const noexcept {
        return *_descriptor;
    }
    const ChunkGeometry& geometry() const noexcept {
        return *_geometry;
    }
    WorkPool& workers() const noexcept {
        return *_workers;
    }
    const AxisMap& map() const noexcept {
        return _map;
    }

private:
    ReducibleImage(const PixelSource& source, const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                  WorkPool& workers, const AxisMap& map)
        : _source(&source), _descriptor(&descriptor), _geometry(&geometry), _workers(&workers), _map(map) {}

    const PixelSource* _source;
    const ImageDescriptor* _descriptor;
    const ChunkGeometry* _geometry;
    WorkPool* _workers;
    AxisMap _map;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_
