/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_REGION_RUNS_H_
#define CARTA_ZARR_SRC_REDUCE_REGION_RUNS_H_

// A region's raster as the runs a reduction walks.
//
// Runs are what the walk would rather have, for the two reasons RegionMask gives: they say which
// chunks a region occupies without reading its raster again, and every pixel of a run is selected,
// so a run is accumulated by the loop an unmasked region uses. They are made here, along whichever
// spatial axis the store varies fastest, so that no caller has to be told which that is; for
// XRADIO, whose m axis is stored last, it is y.
//
// Neither encoder here reads the raster any other way than a row at a time. Along x that is the
// obvious scan, skipping eight bytes at once where nothing changes. Along y every place one row
// differs from the row above is where some column's run starts or stops, and a counting sort by
// column puts those edges back into runs without ever stepping down a column. Both take about 2 ms
// on a 7763x4742 raster where the column-at-a-time walk took 26.

#include <cstdint>
#include <vector>

namespace carta::zarr::internal {

// One region's runs. Line l -- a row of the raster when the runs lie along x, a column when they
// lie along y -- owns the runs at [offsets[l], offsets[l+1]), and run k is the half-open range
// [runs[2k], runs[2k+1]) along the line, counted from the region's own start. Runs within a line are
// disjoint and ascending.
struct RegionRuns {
    std::vector<std::uint32_t> runs;
    std::vector<std::uint64_t> offsets;
};

// The runs of a `width` x `height` raster laid out row-major with x fastest, along x when `u_is_x`
// and along y otherwise. Any non-zero byte selects its pixel.
//
// Returns false, and leaves `into` empty, when the raster breaks into more runs than it is worth:
// more than one per line and more than one per sixteen pixels. A run costs eight bytes and a turn of
// the loop that reads it, so a raster that fragmented is cheaper read as a raster, and a
// checkerboard would otherwise cost several times its own size. No region CARTA draws comes near
// it -- a rectangle, an ellipse or a polygon is one run a line, or a few -- so this is a bound
// against an arbitrary raster, not a tuning.
bool RunsAlongU(const std::uint8_t* raster, std::uint64_t width, std::uint64_t height, bool u_is_x, RegionRuns& into);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_REGION_RUNS_H_
