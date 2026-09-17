/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_READABLE_IMAGE_H_
#define CARTA_ZARR_SRC_READABLE_IMAGE_H_

#include "axis_map.h"
#include "pixel_source.h"
#include "work_pool.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

namespace carta::zarr::internal {

/**
 * One image, opened, with everything a read or a reduction needs before a request arrives.
 *
 * Four things travelled together to every entry point here -- where the pixels come from, what the
 * image is, how it is laid out, and what may run the arithmetic -- and the descriptor travelled
 * twice, once on its own and once inside the source. Four entry points took six or seven arguments
 * of which six were the same six in the same order.
 *
 * It is not only the bundle, because a bundle would just move those arguments rather than absorb
 * them: it also answers what each entry point used to work out again for itself. Where the roles
 * sit among the axes is derived once here rather than four times, and a spectral range is checked
 * against the image through this rather than against a descriptor and a map the caller had to have
 * in hand.
 *
 * Holds references and an AxisMap by value, so it is cheap to build per call and owns nothing. It
 * must not outlive the source, the descriptor, the geometry or the pool it was built from.
 */
class ReadableImage {
public:
    // Reports not_implemented when the image's axes cannot be mapped -- an axis with no known role
    // that is not degenerate, or a missing spatial or spectral axis. Built per call, so that
    // failure reaches the caller of the operation that needed it, which is where it reached before.
    static Result<ReadableImage> Of(const PixelSource& source, const ImageDescriptor& descriptor,
                                    const ChunkGeometry& geometry, WorkPool& workers) {
        auto map = MapAxes(descriptor);
        if (!map) {
            return map.error();
        }
        return ReadableImage(source, descriptor, geometry, workers, map.value());
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

    // The spectral range an operation was asked for has to fall inside this image. Compared by
    // dividing the room that is left rather than by multiplying out the span: (count - 1) * stride
    // wraps, and a wrapped span passes a check it should fail.
    Result<void> ValidateSpectral(const Range& spectral) const {
        const auto channels = _descriptor->axes.at(_map.spectral).length;
        if (spectral.stride == 0 || spectral.count == 0 || spectral.start >= channels ||
            spectral.count - 1 > (channels - 1 - spectral.start) / spectral.stride) {
            return Error{ErrorCode::invalid_argument, "The spectral range falls outside the image",
                         _descriptor->id};
        }
        return {};
    }

private:
    ReadableImage(const PixelSource& source, const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                  WorkPool& workers, const AxisMap& map)
        : _source(&source), _descriptor(&descriptor), _geometry(&geometry), _workers(&workers), _map(map) {}

    const PixelSource* _source;
    const ImageDescriptor* _descriptor;
    const ChunkGeometry* _geometry;
    WorkPool* _workers;
    AxisMap _map;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_READABLE_IMAGE_H_
