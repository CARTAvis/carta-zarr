/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// Which chunks a set of regions occupies, stated directly.
//
// Every claim here used to be reached by running a whole spectral reduction over a 4 x 5 x 2 x 3
// fixture and counting how many times a result was handed over -- two tests that said in their own
// comments that they stop testing anything if the chunk shape ever changes. Two of the refusals
// below were reachable from no test at all.
//
// This target links nothing. The question is regions and a chunk shape in, an index out.

#include "reduce/occupancy.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ErrorCode;
using carta::zarr::RegionMask;
using carta::zarr::internal::Occupancy;

using carta::zarr::testing::Require;

RegionMask Box(std::uint64_t x, std::uint64_t y, std::uint64_t width, std::uint64_t height) {
    RegionMask region;
    region.x_start = x;
    region.y_start = y;
    region.width = width;
    region.height = height;
    return region;
}

Occupancy Built(const std::vector<RegionMask>& regions, std::uint64_t chunk_u, std::uint64_t chunk_v,
                AxisRole fastest = AxisRole::spatial_x) {
    auto built = Occupancy::Of(regions.data(), regions.size(), chunk_u, chunk_v, fastest, "TEST");
    Require(static_cast<bool>(built),
            "Occupancy::Of failed: " + (built ? std::string{} : built.error().message));
    return std::move(built.value());
}

// The cells of the grid that at least one region occupies, as "cu,cv" pairs, read back through
// RegionsTouching rather than through the index that holds them.
std::string Occupied(const Occupancy& occupancy) {
    std::string text;
    for (std::uint64_t cv = 0; cv < occupancy.rows(); ++cv) {
        for (std::uint64_t cu = 0; cu < occupancy.columns(); ++cu) {
            if (occupancy.RegionsTouching(occupancy.chunk_cu0() + cu, occupancy.chunk_cv0() + cv).size > 0) {
                text += (text.empty() ? "" : " ") + std::to_string(cu) + "," + std::to_string(cv);
            }
        }
    }
    return text;
}

// A raster whose pixel (x, y) of the bounding box is set when the chunk it falls in is on the
// diagonal of the chunk grid.
std::vector<std::uint8_t> DiagonalRaster(std::uint64_t width, std::uint64_t height, std::uint64_t chunk) {
    std::vector<std::uint8_t> raster(static_cast<std::size_t>(width * height), 0);
    for (std::uint64_t y = 0; y < height; ++y) {
        for (std::uint64_t x = 0; x < width; ++x) {
            if (x / chunk == y / chunk) {
                raster.at(static_cast<std::size_t>((y * width) + x)) = 1;
            }
        }
    }
    return raster;
}

void TestABoxOccupiesEveryChunkItsBoundingBoxTouches() {
    // x in [2, 8) and y in [3, 8) over a 4 x 4 chunk grid: two chunk columns, two chunk rows, and
    // no mask to narrow either.
    const auto occupancy = Built({Box(2, 3, 6, 5)}, 4, 4);

    Require(occupancy.columns() == 2 && occupancy.rows() == 2, "the chunk grid was not 2 x 2");
    Require(occupancy.u0() == 2 && occupancy.u1() == 8 && occupancy.v0() == 3 && occupancy.v1() == 8,
            "the bounding box was not the region's own");
    Require(Occupied(occupancy) == "0,0 1,0 0,1 1,1", "a box did not occupy its whole bounding box");
    Require(occupancy.LayerChunks() == 4, "a 2 x 2 box did not report four chunks in a layer");
}

void TestThePlacementFollowsTheFastestSpatialAxis() {
    const auto straight = Built({Box(2, 3, 6, 5)}, 4, 4, AxisRole::spatial_x);
    const auto& as_written = straight.regions().at(0);
    Require(as_written.u_start == 2 && as_written.v_start == 3 && as_written.u_size == 6 &&
                as_written.v_size == 5 && as_written.mask_u_stride == 1 && as_written.mask_v_stride == 6,
            "x was not placed on u when the store varies x fastest");

    // The same region on a store that varies y fastest: x and y swap, and so do the steps through
    // the raster the caller wrote in its own order.
    const auto swapped = Built({Box(2, 3, 6, 5)}, 4, 4, AxisRole::spatial_y);
    const auto& placed = swapped.regions().at(0);
    Require(placed.u_start == 3 && placed.v_start == 2 && placed.u_size == 5 && placed.v_size == 6 &&
                placed.mask_u_stride == 6 && placed.mask_v_stride == 1,
            "y was not placed on u when the store varies y fastest");
}

void TestAMaskNarrowsTheOccupancyBelowTheBoundingBox() {
    // The whole reason this module exists: a cut along the diagonal has a bounding box of sixteen
    // chunks and occupies four of them.
    const auto raster = DiagonalRaster(16, 16, 4);
    auto region = Box(0, 0, 16, 16);
    region.mask = raster.data();

    const auto occupancy = Built({region}, 4, 4);

    Require(occupancy.columns() == 4 && occupancy.rows() == 4, "the bounding box was not four chunks square");
    Require(Occupied(occupancy) == "0,0 1,1 2,2 3,3", "the mask did not narrow the occupancy to the diagonal");
    Require(occupancy.LayerChunks() == 4, "a diagonal did not report one chunk per row");
}

void TestRunsAndARasterSayTheSameThing() {
    const auto raster = DiagonalRaster(16, 16, 4);
    auto rastered = Box(0, 0, 16, 16);
    rastered.mask = raster.data();

    // The same selection as runs: row y holds the one run [4 * (y / 4), 4 * (y / 4) + 4).
    std::vector<std::uint32_t> runs;
    std::vector<std::uint64_t> offsets{0};
    for (std::uint64_t y = 0; y < 16; ++y) {
        runs.push_back(static_cast<std::uint32_t>(4 * (y / 4)));
        runs.push_back(static_cast<std::uint32_t>((4 * (y / 4)) + 4));
        offsets.push_back(runs.size() / 2);
    }
    auto run_length = Box(0, 0, 16, 16);
    run_length.row_runs = runs.data();
    run_length.row_run_offsets = offsets.data();

    const auto from_raster = Built({rastered}, 4, 4);
    const auto from_runs = Built({run_length}, 4, 4);

    Require(Occupied(from_raster) == Occupied(from_runs),
            "runs and a raster describing one selection disagreed about the chunks it occupies");
    Require(from_raster.offsets() == from_runs.offsets() && from_raster.entries() == from_runs.entries(),
            "runs and a raster reached different indexes for one selection");
    Require(from_raster.LayerChunks() == from_runs.LayerChunks(),
            "runs and a raster disagreed about how many chunks a layer occupies");
}

void TestTheIncidencesOfOneChunkAreContiguousAndInRegionOrder() {
    // Three regions over one 8 x 8 chunk grid of 4 x 4 chunks. Region 0 covers everything, region 1
    // the top-left chunk, region 2 the bottom-right one. The counting sort is what puts 0 before 1
    // in the first chunk's entries, and the accumulation reads them in that order.
    const auto occupancy = Built({Box(0, 0, 8, 8), Box(0, 0, 4, 4), Box(4, 4, 4, 4)}, 4, 4);

    const auto top_left = occupancy.RegionsTouching(0, 0);
    Require(top_left.size == 2 && top_left.data[0] == 0 && top_left.data[1] == 1,
            "the top-left chunk's regions were not 0 then 1");
    const auto bottom_right = occupancy.RegionsTouching(1, 1);
    Require(bottom_right.size == 2 && bottom_right.data[0] == 0 && bottom_right.data[1] == 2,
            "the bottom-right chunk's regions were not 0 then 2");
    const auto top_right = occupancy.RegionsTouching(1, 0);
    Require(top_right.size == 1 && top_right.data[0] == 0, "the top-right chunk was touched by more than region 0");

    // The same claim against the index itself: offsets rise, and there is one entry per incidence.
    const auto& offsets = occupancy.offsets();
    Require(offsets.size() == 5, "a 2 x 2 grid did not produce five offsets");
    for (std::size_t cell = 0; cell + 1 < offsets.size(); ++cell) {
        Require(offsets.at(cell) <= offsets.at(cell + 1), "the offsets were not non-decreasing");
    }
    Require(offsets.back() == occupancy.entries().size(), "the last offset did not account for every entry");
    Require(occupancy.entries().size() == 6, "four chunks of region 0 plus one each of 1 and 2 is six incidences");
}

void TestAMaskThatSelectsNothingOccupiesNothing() {
    const std::vector<std::uint8_t> raster(16 * 16, 0);
    auto region = Box(0, 0, 16, 16);
    region.mask = raster.data();

    const auto occupancy = Built({region}, 4, 4);

    Require(occupancy.columns() == 4 && occupancy.rows() == 4,
            "a region that selects nothing still has the bounding box it was given");
    Require(occupancy.entries().empty(), "a mask of zeroes produced incidences");
    Require(Occupied(occupancy).empty(), "a mask of zeroes occupied a chunk");
    Require(occupancy.LayerChunks() == 0, "a mask of zeroes did not report an empty layer");
    for (const auto& runs : occupancy.runs_per_row()) {
        Require(runs.empty(), "a mask of zeroes produced a column run");
    }
}

void TestRunsAlongTheOtherAxisAreRefused() {
    std::vector<std::uint32_t> runs{0, 4};
    std::vector<std::uint64_t> offsets{0, 1, 1, 1, 1};
    auto region = Box(0, 0, 4, 4);
    region.row_runs = runs.data();
    region.row_run_offsets = offsets.data();
    region.run_axis = AxisRole::spatial_y;

    const std::vector<RegionMask> regions{region};
    const auto refused = Occupancy::Of(regions.data(), regions.size(), 4, 4, AxisRole::spatial_x, "TEST");
    Require(!refused && refused.error().code == ErrorCode::invalid_argument,
            "runs along the axis the store does not vary fastest were accepted");
}

// kMaxChunkIncidences, the other refusal, is not here, and the measurement is why: the smallest
// input that trips it is 64 regions over a 1024 x 1024 chunk grid plus one more, which took 2.06 s
// and 815 MB of resident memory to reach the check. The refusal exists to stop an allocation
// nothing can serve, and a test of it has to make that allocation first. It was run once by hand
// against this file's Occupancy::Of and reported "The regions together touch more chunks than one
// reduction can index"; that is recorded here rather than paid for on every ctest run.

void TestAGridTooLargeToIndexIsRefused() {
    // 65,536 chunks on a side is 2^32 cells, one more than a cell number can hold. Refused before
    // anything is allocated for it.
    const std::vector<RegionMask> regions{Box(0, 0, 65536, 65536)};
    const auto refused = Occupancy::Of(regions.data(), regions.size(), 1, 1, AxisRole::spatial_x, "TEST");
    Require(!refused && refused.error().code == ErrorCode::invalid_argument,
            "a grid of 2^32 chunks was accepted");
}

}  // namespace

int main() {
    try {
        TestABoxOccupiesEveryChunkItsBoundingBoxTouches();
        TestThePlacementFollowsTheFastestSpatialAxis();
        TestAMaskNarrowsTheOccupancyBelowTheBoundingBox();
        TestRunsAndARasterSayTheSameThing();
        TestTheIncidencesOfOneChunkAreContiguousAndInRegionOrder();
        TestAMaskThatSelectsNothingOccupiesNothing();
        TestRunsAlongTheOtherAxisAreRefused();
        TestAGridTooLargeToIndexIsRefused();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "occupancy test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
