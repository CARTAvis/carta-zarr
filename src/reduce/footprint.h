/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_FOOTPRINT_H_
#define CARTA_ZARR_SRC_REDUCE_FOOTPRINT_H_

#include <cstdint>

namespace carta::zarr::internal {

/**
 * One spatial footprint a run of slabs is read over, in the pass's own axes.
 *
 * `chunks` is how many chunks the footprint covers in one chunk of the spectral axis. It does two
 * jobs: it sizes the slab, because the budget counts chunk data rather than the pixels kept, and it
 * is the unit progress is counted in. If it is wrong, the right pixels are read but the reported
 * progress is wrong.
 *
 * In a header of its own because two modules use it: the walk reads a footprint, and the occupancy
 * builds a reduction's. Neither has to include the other.
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
