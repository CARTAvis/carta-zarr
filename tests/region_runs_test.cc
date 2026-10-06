/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// Turning a region's raster into the runs a reduction walks.
//
// Every answer here is checked against the obvious encoder -- walk each line, note where selection
// starts and stops -- which is exactly the one the library does not use: along y it reads the raster
// a column at a time, a stride of a whole row per byte, and that is the cost the real encoders exist
// to avoid. The oracle is allowed to be slow because it only has to be right.

#include "reduce/region_runs.h"

#include "support/check.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace {

using carta::zarr::internal::RegionRuns;
using carta::zarr::internal::RunsAlongU;

using carta::zarr::testing::Require;

// Lines are rows when u is x and columns when it is y; runs are ranges along each line.
RegionRuns Oracle(const std::vector<std::uint8_t>& raster, std::uint64_t width, std::uint64_t height, bool u_is_x) {
    RegionRuns out;
    out.offsets.push_back(0);
    const std::uint64_t lines = u_is_x ? height : width;
    const std::uint64_t length = u_is_x ? width : height;
    for (std::uint64_t line = 0; line < lines; ++line) {
        const auto at = [&](std::uint64_t i) {
            return raster.at(u_is_x ? (line * width) + i : (i * width) + line) != 0;
        };
        for (std::uint64_t i = 0; i < length;) {
            if (!at(i)) {
                ++i;
                continue;
            }
            const auto begin = i;
            while (i < length && at(i)) {
                ++i;
            }
            out.runs.push_back(static_cast<std::uint32_t>(begin));
            out.runs.push_back(static_cast<std::uint32_t>(i));
        }
        out.offsets.push_back(out.runs.size() / 2);
    }
    return out;
}

void RequireSame(const RegionRuns& expected, const RegionRuns& actual, const std::string& what) {
    Require(actual.offsets == expected.offsets, what + ": the lines own different runs");
    Require(actual.runs == expected.runs, what + ": the runs differ");
}

// Checks the runs, not whether they are worth taking: the raster is given empty rows below it until
// it has room for sixteen pixels a run either way, which is what RunsAlongU asks of any raster before
// it takes it. Only the refusal test below is about that bound.
void RequireBothWays(std::vector<std::uint8_t> raster, std::uint64_t width, std::uint64_t height,
                     const std::string& what) {
    const auto runs_of = [&](bool u_is_x) { return Oracle(raster, width, height, u_is_x).runs.size() / 2; };
    const std::uint64_t most = std::max(runs_of(true), runs_of(false));
    while (width * height < 16 * most) {
        raster.resize(raster.size() + width, 0);
        ++height;
    }
    for (const bool u_is_x : {true, false}) {
        RegionRuns runs;
        const std::string way = what + (u_is_x ? " along x" : " along y");
        Require(RunsAlongU(raster.data(), width, height, u_is_x, runs), way + ": should be taken as runs");
        RequireSame(Oracle(raster, width, height, u_is_x), runs, way);
    }
}

// A raster drawn from a picture, one string per row, '#' selected.
std::vector<std::uint8_t> Drawn(const std::vector<std::string>& rows) {
    std::vector<std::uint8_t> raster;
    for (const auto& row : rows) {
        for (const char c : row) {
            raster.push_back(c == '#' ? 1 : 0);
        }
    }
    return raster;
}

void TestAShapeWithAHoleAndAnEmptyLine() {
    const std::vector<std::string> picture{
        "..####....#####..", ".######..#######.", ".................",
        "##..##..##..##..#", "#################", "......#..........",
    };
    RequireBothWays(Drawn(picture), picture.front().size(), picture.size(), "holes");
}

// Widths either side of the eight bytes the x encoder reads at a time, so that every way a row can
// end -- inside a word, on its edge, one past it -- is a way some run ends.
void TestEveryWayALineCanEnd() {
    for (const std::uint64_t width : {1, 7, 8, 9, 15, 16, 17, 31}) {
        for (const std::uint64_t height : {1, 3, 8, 9}) {
            std::vector<std::uint8_t> raster(width * height, 0);
            for (std::uint64_t y = 0; y < height; ++y) {
                for (std::uint64_t x = 0; x < width; ++x) {
                    // A different pattern per row: runs that start at 0, end at the edge, and
                    // straddle each word boundary somewhere.
                    raster.at((y * width) + x) = ((x + y) % 5 < 3 || x + 1 == width) ? 1 : 0;
                }
            }
            RequireBothWays(raster, width, height, std::to_string(width) + "x" + std::to_string(height));
        }
    }
}

// The shape a region most often is, at sizes thin in either direction and one of each that is not
// a multiple of anything.
void TestEllipsesOfEveryShape() {
    for (const auto& [width, height] : std::vector<std::pair<std::uint64_t, std::uint64_t>>{
             {1, 1}, {2, 1}, {5, 40}, {40, 5}, {97, 61}, {128, 128}, {300, 7}}) {
        std::vector<std::uint8_t> raster(width * height, 0);
        for (std::uint64_t y = 0; y < height; ++y) {
            for (std::uint64_t x = 0; x < width; ++x) {
                const double dx = (static_cast<double>(x) + 0.5 - (width / 2.0)) / (width / 2.0);
                const double dy = (static_cast<double>(y) + 0.5 - (height / 2.0)) / (height / 2.0);
                raster.at((y * width) + x) = (dx * dx) + (dy * dy) <= 1.0 ? 1 : 0;
            }
        }
        RequireBothWays(raster, width, height, "ellipse " + std::to_string(width) + "x" + std::to_string(height));
    }
}

void TestAFullAndAnEmptyRaster() {
    RequireBothWays(std::vector<std::uint8_t>(13 * 11, 1), 13, 11, "full");
    RequireBothWays(std::vector<std::uint8_t>(13 * 11, 0), 13, 11, "empty");
}

// The mask's contract is that any non-zero byte selects, not that a selected byte is one.
void TestAnyNonZeroByteSelects() {
    // Blank rows below give the fragmented third row room, as RequireBothWays does.
    std::vector<std::string> picture{"..####....", ".#######..", "#.#.#.#.#."};
    picture.resize(24, "..........");
    auto ones = Drawn(picture);
    auto mixed = ones;
    for (std::size_t i = 0; i < mixed.size(); ++i) {
        if (mixed.at(i) != 0) {
            mixed.at(i) = static_cast<std::uint8_t>(i % 2 == 0 ? 0xFF : 0x80 + (i % 7));
        }
    }
    for (const bool u_is_x : {true, false}) {
        RegionRuns from_ones;
        RegionRuns from_mixed;
        Require(RunsAlongU(ones.data(), 10, 24, u_is_x, from_ones), "ones should be taken as runs");
        Require(RunsAlongU(mixed.data(), 10, 24, u_is_x, from_mixed), "mixed bytes should be taken as runs");
        RequireSame(from_ones, from_mixed, u_is_x ? "mixed bytes along x" : "mixed bytes along y");
    }
}

// A raster that breaks into more runs than it is worth is refused, and left as it was found: the
// reduction then reads the raster as a raster. A checkerboard is the worst case, a run per pixel.
void TestAFragmentedRasterIsLeftARaster() {
    const std::uint64_t side = 64;
    std::vector<std::uint8_t> board(side * side, 0);
    for (std::uint64_t y = 0; y < side; ++y) {
        for (std::uint64_t x = 0; x < side; ++x) {
            board.at((y * side) + x) = (x + y) % 2 == 0 ? 1 : 0;
        }
    }
    for (const bool u_is_x : {true, false}) {
        RegionRuns runs;
        Require(!RunsAlongU(board.data(), side, side, u_is_x, runs), "a checkerboard should not be taken as runs");
        Require(runs.runs.empty() && runs.offsets.empty(), "and nothing is left behind");
    }
}

// However small a region, a convex one is one run a line and is always worth taking -- nine pixels
// are not room for sixteen a run, and that must not matter. Asked directly, with no rows added.
void TestASmallConvexRegionIsAlwaysRuns() {
    const auto cross = Drawn({".#.", "###", ".#."});
    for (const bool u_is_x : {true, false}) {
        RegionRuns runs;
        Require(RunsAlongU(cross.data(), 3, 3, u_is_x, runs), "a three-pixel cross should be taken as runs");
        RequireSame(Oracle(cross, 3, 3, u_is_x), runs, "a three-pixel cross");
    }
}

}  // namespace

int main() {
    try {
        TestAShapeWithAHoleAndAnEmptyLine();
        TestEveryWayALineCanEnd();
        TestEllipsesOfEveryShape();
        TestAFullAndAnEmptyRaster();
        TestAnyNonZeroByteSelects();
        TestAFragmentedRasterIsLeftARaster();
        TestASmallConvexRegionIsAlwaysRuns();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "region runs test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
