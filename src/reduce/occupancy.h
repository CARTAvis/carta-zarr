/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_OCCUPANCY_H_
#define CARTA_ZARR_SRC_REDUCE_OCCUPANCY_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace carta::zarr::internal {

// One region placed on the axes the walk reads in: u is the spatial axis the store varies fastest
// and v is the other, so that a plane arrives with u contiguous and never has to be transposed on
// the way in. A caller's x and y are mapped onto these once, at the top of the reduction.
//
// Placed rather than "as the walk sees it", which is what it used to be called: a pass never knows
// what a region is, so this is what a region looks like after the placement rather than something
// the walk holds.
struct PlacedRegion {
    std::uint64_t u_start = 0;
    std::uint64_t v_start = 0;
    std::uint64_t u_size = 0;
    std::uint64_t v_size = 0;
    // Steps through the raster for one pixel of u and of v. One of them is 1; which one depends on
    // whether the caller's rows run along u or across it.
    const std::uint8_t* mask = nullptr;
    std::uint64_t mask_u_stride = 1;
    std::uint64_t mask_v_stride = 1;
    // Runs along u, indexed by v. Null unless the caller supplied runs; Occupancy::Of refuses a
    // region whose runs go the other way, so anything that reaches here already runs along u.
    const std::uint32_t* runs = nullptr;
    const std::uint64_t* run_offsets = nullptr;
};

// A maximal run of consecutive chunk columns that at least one region touches, in bounding-box
// grid coordinates.
struct ColumnRun {
    std::uint64_t first = 0;
    std::uint64_t last = 0;  // inclusive

    bool operator==(const ColumnRun& other) const {
        return first == other.first && last == other.last;
    }
};

// The regions touching one chunk, as a view into the index that holds them. Valid for as long as
// the Occupancy is.
//
// A view rather than a pair of iterators because the accumulation reads one per chunk cell and
// wants to write a range-for; the same shape BufferView carries, for the same reason.
struct RegionRefs {
    const std::uint32_t* data = nullptr;
    std::size_t size = 0;

    const std::uint32_t* begin() const noexcept {
        return data;
    }
    const std::uint32_t* end() const noexcept {
        return data + size;
    }
};

/**
 * Which chunks a set of regions occupies, and which of them touch each one.
 *
 * Occupies rather than covers: the bounding box of a thin cut laid along the diagonal is the whole
 * image, while the cut touches one chunk per row. Bucketing by the box would read every chunk to
 * reach the band, so the region's own mask or runs decide instead -- one pass over bytes the caller
 * already holds, against the chunks they would otherwise stand for.
 *
 * Built once per reduction and read for the whole of it. Without it the walk is chunks x regions
 * intersection tests -- 6.5 x 10^8 for a WSU diagonal -- which would put the region count back on
 * the critical path the spectral reduction exists to clear.
 *
 * It takes the caller's regions in the caller's own x and y, because placing them onto the walk's
 * axes is part of the same question: which axis the store varies fastest is what decides both where
 * a region lands and whether its runs are the kind worth taking. Said in one place, a PlacedRegion
 * that exists is one whose runs have been agreed to.
 *
 * Takes chunk_u, chunk_v and the fastest spatial axis rather than a PassPlan, because those are
 * what the question is about and what a test can stand up with nothing linked behind it -- ADR 0006,
 * and the same reasoning CheckedPlanes gives for taking a descriptor rather than a ReadableImage.
 */
class Occupancy {
public:
    // Reports invalid_argument for a region whose runs go along the axis the store does not vary
    // fastest, and for a region set touching more chunks than one reduction can index.
    //
    // Assumes at least one region: an empty set is refused a step earlier, where the rest of the
    // request is checked.
    static Result<Occupancy> Of(const RegionMask* regions, std::size_t region_count, std::uint64_t chunk_u,
                                std::uint64_t chunk_v, AxisRole fastest_spatial_axis, const std::string& node);

    // The caller's regions on the walk's axes, in the order they were given. The accumulation
    // indexes this by what RegionsTouching hands back.
    const std::vector<PlacedRegion>& regions() const noexcept {
        return _regions;
    }

    // The occupied runs of each chunk row, which is what the walk reads instead of the bounding box.
    const std::vector<std::vector<ColumnRun>>& runs_per_row() const noexcept {
        return _runs_per_row;
    }

    // The chunks one spectral layer of the whole region set occupies. Zero when the regions select
    // nothing at all, which a mask of zeroes does.
    std::uint64_t LayerChunks() const noexcept {
        return _layer_chunks;
    }

    std::size_t Cell(std::uint64_t chunk_cu, std::uint64_t chunk_cv) const noexcept {
        return static_cast<std::size_t>(((chunk_cv - _chunk_cv0) * _columns) + (chunk_cu - _chunk_cu0));
    }

    // The regions touching one chunk. Read once per chunk cell by the accumulation, so it is inline
    // and does nothing but two lookups.
    RegionRefs RegionsTouching(std::uint64_t chunk_cu, std::uint64_t chunk_cv) const {
        const auto cell = Cell(chunk_cu, chunk_cv);
        const auto first = _offsets.at(cell);
        const auto last = _offsets.at(cell + 1);
        return RegionRefs{_entries.data() + first, static_cast<std::size_t>(last - first)};
    }

    // The compressed-row index itself. Exposed because it is what Of promises -- the incidences of
    // one chunk are contiguous and in region order, which is what the counting sort is for -- and a
    // test has nothing else to say that against.
    const std::vector<std::uint64_t>& offsets() const noexcept {
        return _offsets;
    }
    const std::vector<std::uint32_t>& entries() const noexcept {
        return _entries;
    }

    // The union bounding box of every region, in pixels of the walk's own axes, and the chunk grid
    // covering it. The band walk clamps its footprints against these.
    std::uint64_t u0() const noexcept { return _u0; }
    std::uint64_t v0() const noexcept { return _v0; }
    std::uint64_t u1() const noexcept { return _u1; }
    std::uint64_t v1() const noexcept { return _v1; }
    std::uint64_t chunk_cu0() const noexcept { return _chunk_cu0; }
    std::uint64_t chunk_cv0() const noexcept { return _chunk_cv0; }
    std::uint64_t columns() const noexcept { return _columns; }
    std::uint64_t rows() const noexcept { return _rows; }

private:
    Occupancy() = default;

    std::vector<PlacedRegion> _regions;
    std::uint64_t _u0 = 0;
    std::uint64_t _v0 = 0;
    std::uint64_t _u1 = 0;
    std::uint64_t _v1 = 0;
    std::uint64_t _chunk_cu0 = 0;
    std::uint64_t _chunk_cv0 = 0;
    std::uint64_t _columns = 0;
    std::uint64_t _rows = 0;
    std::vector<std::uint64_t> _offsets;
    std::vector<std::uint32_t> _entries;
    std::vector<std::vector<ColumnRun>> _runs_per_row;
    std::uint64_t _layer_chunks = 0;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_OCCUPANCY_H_
