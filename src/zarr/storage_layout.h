/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_ZARR_STORAGE_LAYOUT_H_
#define CARTA_ZARR_SRC_ZARR_STORAGE_LAYOUT_H_

#include <cstdint>
#include <string>
#include <vector>

namespace carta::zarr::internal::zarr {

// How an array is laid out, in the order its own dimensions are stored. What a consumer is told is
// ChunkGeometry, the same facts permuted into logical order; this is what that is built from, and
// it is not reported beside it, so a consumer never has to pick between one fact in two orders.
struct StorageLayout {
    std::vector<std::uint64_t> chunk_shape;
    std::vector<std::uint64_t> shard_shape;
    std::string compressor;
    bool sharded = false;
};

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_STORAGE_LAYOUT_H_
