/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_
#define CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include "axis_map.h"

#include <cstdint>

namespace carta::zarr::internal {

/**
 * A plane selection that an image has agreed to serve.
 *
 * The three fields of a `PlaneSelection` answer one question: can this image give these planes? It
 * is checked here, once, and the answer is a value rather than a `Result<void>`. The only way to get
 * a `CheckedPlanes` is to check a selection, so `PlanPass` and the reductions, which take one, cannot
 * plan a pass over unchecked planes.
 *
 * Takes a descriptor and an `AxisMap` rather than a `ReducibleImage`, because those are what the
 * check is about, and a test can build them with nothing linked behind them (ADR 0006).
 */
class CheckedPlanes {
public:
    // Reports invalid_slice for a spectral range that does not fit, for a polarization or time
    // index outside its axis, and for a non-zero index on an axis the image does not have.
    static Result<CheckedPlanes> Of(const ImageDescriptor& descriptor, const AxisMap& map,
                                    const PlaneSelection& planes) {
        const auto& node = descriptor.id;

        // Divides the room that is left rather than multiplying out the span: (count - 1) * stride
        // can wrap, and a wrapped span would pass a check it should fail.
        const auto channels = descriptor.axes.at(map.spectral).length;
        if (planes.spectral.stride == 0 || planes.spectral.count == 0 || planes.spectral.start >= channels ||
            planes.spectral.count - 1 > (channels - 1 - planes.spectral.start) / planes.spectral.stride) {
            return Error{ErrorCode::invalid_slice, "The spectral range falls outside the image", node};
        }

        if (!map.has_polarization && planes.polarization != 0) {
            return Error{ErrorCode::invalid_slice,
                         "The image has no polarization axis, so only polarization 0 selects a plane", node};
        }
        if (map.has_polarization && planes.polarization >= descriptor.axes.at(map.polarization).length) {
            return Error{ErrorCode::invalid_slice, "The polarization index falls outside the image", node};
        }
        if (!map.has_time && planes.time != 0) {
            return Error{ErrorCode::invalid_slice, "The image has no time axis, so only time 0 selects a plane", node};
        }
        if (map.has_time && planes.time >= descriptor.axes.at(map.time).length) {
            return Error{ErrorCode::invalid_slice, "The time index falls outside the image", node};
        }

        return CheckedPlanes(planes);
    }

    const PlaneSelection& selection() const noexcept { return _planes; }
    const Range& spectral() const noexcept { return _planes.spectral; }
    // How many planes were selected, which is what every block loop over the selection counts to.
    std::uint64_t count() const noexcept { return _planes.spectral.count; }

private:
    explicit CheckedPlanes(const PlaneSelection& planes) : _planes(planes) {}

    PlaneSelection _planes;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PLANE_SELECTION_H_
