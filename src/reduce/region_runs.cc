/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "region_runs.h"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <limits>

namespace carta::zarr::internal {
namespace {

constexpr std::uint64_t kLeastPixelsPerRun = 16;
constexpr std::uint64_t kLowBits = 0x0101010101010101ULL;
constexpr std::uint64_t kHighBits = 0x8080808080808080ULL;

std::uint64_t Word(const std::uint8_t* at) {
    std::uint64_t word = 0;
    std::memcpy(&word, at, sizeof(word));
    return word;
}

// Whether no byte of the word is zero, which is "every pixel selected" whatever the non-zero bytes
// happen to be.
bool NoZeroByte(std::uint64_t word) {
    return ((word - kLowBits) & ~word & kHighBits) == 0;
}

void Discard(RegionRuns& into) {
    into.runs.clear();
    into.offsets.clear();
}

// Along x: each row is a line, read in order.
bool AlongX(const std::uint8_t* raster, std::uint64_t width, std::uint64_t height, std::uint64_t most_runs,
            RegionRuns& into) {
    into.offsets.reserve(static_cast<std::size_t>(height) + 1);
    into.offsets.push_back(0);
    for (std::uint64_t y = 0; y < height; ++y) {
        const std::uint8_t* row = raster + (y * width);
        bool inside = false;
        std::uint64_t begin = 0;
        std::uint64_t x = 0;
        while (x < width) {
            // Eight bytes at once where nothing changes: all unselected outside a run, all selected
            // inside one.
            if (x + sizeof(std::uint64_t) <= width) {
                const std::uint64_t word = Word(row + x);
                if (inside ? NoZeroByte(word) : word == 0) {
                    x += sizeof(std::uint64_t);
                    continue;
                }
            }
            const bool selected = row[x] != 0;
            if (selected != inside) {
                if (selected) {
                    begin = x;
                } else {
                    into.runs.push_back(static_cast<std::uint32_t>(begin));
                    into.runs.push_back(static_cast<std::uint32_t>(x));
                }
                inside = selected;
            }
            ++x;
        }
        if (inside) {
            into.runs.push_back(static_cast<std::uint32_t>(begin));
            into.runs.push_back(static_cast<std::uint32_t>(width));
        }
        if (into.runs.size() / 2 > most_runs) {
            Discard(into);
            return false;
        }
        into.offsets.push_back(into.runs.size() / 2);
    }
    return true;
}

// Along y: each column is a line, and the raster is still read a row at a time. A byte whose
// selection differs from the byte above it is an edge of its column's run -- where one starts, or
// where one stopped at the row before -- and the rows arrive in order, so each column's edges do
// too. A counting sort by column keeps that order, and consecutive edges of a column are its runs.
bool AlongY(const std::uint8_t* raster, std::uint64_t width, std::uint64_t height, std::uint64_t most_runs,
            RegionRuns& into) {
    const std::uint64_t most_edges = 2 * most_runs;
    std::vector<std::uint32_t> edge_columns;
    std::vector<std::uint32_t> edge_rows;
    const std::vector<std::uint8_t> nothing(static_cast<std::size_t>(width), 0);
    const std::uint8_t* above = nothing.data();
    // One row past the last, compared with nothing, closes every run still open at the bottom.
    for (std::uint64_t y = 0; y <= height; ++y) {
        const std::uint8_t* row = y < height ? raster + (y * width) : nothing.data();
        std::uint64_t x = 0;
        while (x < width) {
            if (x + sizeof(std::uint64_t) <= width && Word(row + x) == Word(above + x)) {
                x += sizeof(std::uint64_t);
                continue;
            }
            if ((row[x] != 0) != (above[x] != 0)) {
                edge_columns.push_back(static_cast<std::uint32_t>(x));
                edge_rows.push_back(static_cast<std::uint32_t>(y));
            }
            ++x;
        }
        if (edge_columns.size() > most_edges) {
            return false;
        }
        above = row;
    }

    into.offsets.assign(static_cast<std::size_t>(width) + 1, 0);
    for (const auto column : edge_columns) {
        ++into.offsets.at(static_cast<std::size_t>(column) + 1);
    }
    for (std::size_t column = 0; column < width; ++column) {
        into.offsets.at(column + 1) += into.offsets.at(column);
    }
    into.runs.resize(edge_rows.size());
    std::vector<std::uint64_t> cursor(into.offsets.begin(), into.offsets.end() - 1);
    for (std::size_t k = 0; k < edge_rows.size(); ++k) {
        into.runs.at(static_cast<std::size_t>(cursor.at(edge_columns.at(k))++)) = edge_rows.at(k);
    }
    // Counted in edges so far; two edges make a run.
    for (auto& offset : into.offsets) {
        offset /= 2;
    }
    return true;
}

}  // namespace

bool RunsAlongU(const std::uint8_t* raster, std::uint64_t width, std::uint64_t height, bool u_is_x, RegionRuns& into) {
    Discard(into);
    const std::uint64_t lines = u_is_x ? height : width;
    // A run is stored in 32 bits, and so is a column while the y encoder sorts its edges. A region
    // wider or taller than that is not one anyone draws, and is read as the raster it is.
    if (std::max(width, height) > std::numeric_limits<std::uint32_t>::max()) {
        return false;
    }
    const std::uint64_t most_runs = std::max(lines, (width * height) / kLeastPixelsPerRun);
    const bool taken =
        u_is_x ? AlongX(raster, width, height, most_runs, into) : AlongY(raster, width, height, most_runs, into);
    if (!taken) {
        Discard(into);
    }
    return taken;
}

}  // namespace carta::zarr::internal
