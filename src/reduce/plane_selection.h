/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_
#define CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_

#include "axis_map.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include <cstdint>

namespace carta::zarr::internal {

/**
 * A plane selection that an image has agreed to serve.
 *
 * The three fields of a `PlaneSelection` are one question -- can this image give me these planes --
 * and it was being asked in pieces. The spectral third was checked by the readable image, the other
 * two thirds by a stanza written out once per reduction in two different wordings, and a cube
 * histogram, which has no bins or bounds of its own to check, built a stand-in `HistogramRequest`
 * with `bins = 1` over `[0, 1]` for no reason but to reach one of those stanzas. `cad6c8a` fixed the
 * strided-range guard in three files.
 *
 * So the question is asked here, once, and the answer is a value rather than a `Result<void>`: a
 * `CheckedPlanes` is the only way to say a selection has been checked, and the only way to get one
 * is to check it. `PlanPass` and both reductions take one instead of three loose numbers, which is
 * how they stop being able to plan a pass over planes nobody looked at.
 *
 * Takes a descriptor and an `AxisMap` rather than a `ReadableImage`, because those are what the
 * question is about and what a test can stand up with nothing linked behind it -- ADR 0006.
 */
class CheckedPlanes {
public:
    // Reports invalid_argument for a spectral range that does not fit, for a polarization or time
    // index outside its axis, and for a non-zero index on an axis the image does not have.
    static Result<CheckedPlanes> Of(const ImageDescriptor& descriptor, const AxisMap& map,
                                    const PlaneSelection& planes) {
        const auto& node = descriptor.id;

        // Compared by dividing the room that is left rather than by multiplying out the span:
        // (count - 1) * stride wraps, and a wrapped span passes a check it should fail.
        const auto channels = descriptor.axes.at(map.spectral).length;
        if (planes.spectral.stride == 0 || planes.spectral.count == 0 || planes.spectral.start >= channels ||
            planes.spectral.count - 1 > (channels - 1 - planes.spectral.start) / planes.spectral.stride) {
            return Error{ErrorCode::invalid_argument, "The spectral range falls outside the image", node};
        }

        if (!map.has_polarization && planes.polarization != 0) {
            return Error{ErrorCode::invalid_argument,
                         "The image has no polarization axis, so only polarization 0 selects a plane", node};
        }
        if (map.has_polarization && planes.polarization >= descriptor.axes.at(map.polarization).length) {
            return Error{ErrorCode::invalid_argument, "The polarization index falls outside the image", node};
        }
        if (!map.has_time && planes.time != 0) {
            return Error{ErrorCode::invalid_argument,
                         "The image has no time axis, so only time 0 selects a plane", node};
        }
        if (map.has_time && planes.time >= descriptor.axes.at(map.time).length) {
            return Error{ErrorCode::invalid_argument, "The time index falls outside the image", node};
        }

        return CheckedPlanes(planes);
    }

    const PlaneSelection& selection() const noexcept {
        return _planes;
    }
    const Range& spectral() const noexcept {
        return _planes.spectral;
    }
    // How many planes were selected, which is what every block loop over the selection counts to.
    std::uint64_t count() const noexcept {
        return _planes.spectral.count;
    }

private:
    explicit CheckedPlanes(const PlaneSelection& planes) : _planes(planes) {}

    PlaneSelection _planes;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_
