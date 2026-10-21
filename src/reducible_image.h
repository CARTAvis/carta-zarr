/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_
#define CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"
#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include "axis_map.h"
#include "pixel_source.h"
#include "reduce/pass_plan.h"
#include "reduce/plane_selection.h"
#include "reduce/task_split.h"
#include "work_pool.h"

namespace carta::zarr::internal {

/**
 * One image, opened, with everything a reduction needs before a request arrives.
 *
 * Every reduction needs where the pixels come from, what the image is, how it is laid out, and what
 * may run the arithmetic. This holds all four, and also works out once what each reduction would
 * otherwise derive for itself, such as where the roles sit among the axes.
 *
 * An ordinary read is not a reduction and does not come through here: ReadInPieces takes only the
 * three things it uses. See ADR 0005 for why a read stays outside the pass.
 *
 * `Plan` checks a plane selection (with `CheckedPlanes::Of`) and plans the pass over it in one step,
 * because every reduction does both.
 *
 * Holds references and an AxisMap by value, so it is cheap to build per call and owns nothing. It
 * must not outlive the source, the descriptor, the geometries or the pool it was built from.
 */
class ReducibleImage {
public:
    // Reports not_implemented when the image's axes cannot be mapped: an axis with no known role that
    // is not degenerate, or a missing spatial or spectral axis. Built per call, so that the failure
    // reaches only the caller of the reduction that needed it.
    static Result<ReducibleImage> Of(const PixelSource& source, const ImageDescriptor& descriptor,
                                     const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                                     WorkPool& workers) {
        auto map = MapAxes(descriptor);
        if (!map) {
            return map.error();
        }
        return ReducibleImage(source, descriptor, geometry, flag_geometry, workers, map.value());
    }

    const PixelSource& source() const noexcept { return *_source; }
    const ImageDescriptor& descriptor() const noexcept { return *_descriptor; }
    const AxisMap& map() const noexcept { return _map; }

    // The pass a reduction makes over `planes`, once they are checked against this image. The plan
    // keeps the checked planes, so it is the one value a reduction needs afterwards. `sample` takes
    // every nth pixel along both spatial axes; one reads them all.
    //
    // Reports invalid_slice for planes this image cannot serve; see CheckedPlanes::Of. The rest of a
    // request (bins, bounds, regions, a sink) stays with the reduction that understands it.
    Result<PassPlan> Plan(const PlaneSelection& planes, std::uint64_t sample, const ReadOptions& options) const {
        const auto checked = CheckedPlanes::Of(*_descriptor, _map, planes);
        if (!checked) {
            return checked.error();
        }
        return PlanPass(*_descriptor, *_geometry, *_flag_geometry, _map, checked.value(), sample, options,
                        _workers->size());
    }

    // How a reduction whose tasks each hold an accumulator of `accumulator_bytes` divides a read
    // among the pool, with `budget_bytes` for all of them together. This is the only way a reduction
    // reaches the pool, so ADR 0005's cap applies wherever tasks run. See TaskSplit.
    TaskSplit Split(std::size_t budget_bytes, std::size_t accumulator_bytes) const {
        return TaskSplit(*_workers, budget_bytes, accumulator_bytes);
    }

private:
    ReducibleImage(const PixelSource& source, const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                   const ChunkGeometry& flag_geometry, WorkPool& workers, const AxisMap& map)
        : _source(&source),
          _descriptor(&descriptor),
          _geometry(&geometry),
          _flag_geometry(&flag_geometry),
          _workers(&workers),
          _map(map) {}

    const PixelSource* _source;
    const ImageDescriptor* _descriptor;
    const ChunkGeometry* _geometry;
    // The flag's geometry, which a masked pass also sizes its reads by. See DescribedImage.
    const ChunkGeometry* _flag_geometry;
    WorkPool* _workers;
    AxisMap _map;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCIBLE_IMAGE_H_
