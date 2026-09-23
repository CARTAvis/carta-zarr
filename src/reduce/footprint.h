/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_FOOTPRINT_H_
#define CARTA_ZARR_SRC_REDUCE_FOOTPRINT_H_

#include <cstdint>

namespace carta::zarr::internal {

/**
 * One spatial footprint a run of slabs is read over, in the pass's own axes.
 *
 * `chunks` is what that footprint occupies in one chunk of the spectral axis, and it is doing two
 * jobs: it sizes the slab, because the budget is over chunk data rather than over the pixels kept,
 * and it is the unit progress is counted in. A caller that gets it wrong reads the right pixels and
 * reports a bar that lies.
 *
 * In a header of its own because two modules speak it: the walk reads a footprint, and the
 * occupancy is what forms a reduction's. Neither should have to include the other to say what one
 * is.
 */
struct SlabFootprint {
    std::uint64_t u_start = 0;
    std::uint64_t u_count = 0;
    std::uint64_t u_stride = 1;
    std::uint64_t v_start = 0;
    std::uint64_t v_count = 0;
    std::uint64_t v_stride = 1;
    std::uint64_t chunks = 1;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_FOOTPRINT_H_
