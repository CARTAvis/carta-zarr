/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "plan.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <string>
#include <unordered_set>
#include <utility>

namespace carta::zarr::bench {

namespace {

constexpr std::uint64_t kQuietNaN32 = 0x7FC00000u;
constexpr std::uint64_t kQuietNaN64 = 0x7FF8000000000000ull;
constexpr std::uint64_t kFnvPrime = 0x100000001B3ull;

}  // namespace

std::vector<Draw> DistinctSample(Stream& stream, std::uint64_t pool, std::uint64_t count) {
    std::vector<Draw> draws(count);
    if (pool == 0) {
        return draws;
    }
    constexpr std::uint64_t kMaterialize = std::uint64_t{1} << 20;
    if (pool <= kMaterialize || count > pool / 2) {
        // A partial Fisher-Yates shuffle, which only ever touches the first min(count, pool) slots.
        std::vector<std::uint64_t> values(pool);
        std::iota(values.begin(), values.end(), std::uint64_t{0});
        const auto distinct = std::min(count, pool);
        for (std::uint64_t index = 0; index < distinct; ++index) {
            std::swap(values[index], values[index + stream.Below(pool - index)]);
        }
        const auto repeats = count > pool ? count - pool : 0;
        for (std::uint64_t index = 0; index < count; ++index) {
            draws[index] = {values[index % pool], index < repeats || index >= pool};
        }
        return draws;
    }
    std::unordered_set<std::uint64_t> taken;
    for (auto& draw : draws) {
        do {
            draw.value = stream.Below(pool);
        } while (!taken.insert(draw.value).second);
    }
    return draws;
}

void Mix(std::uint64_t& hash, std::uint64_t bits, unsigned bytes) {
    for (unsigned byte = 0; byte < bytes; ++byte) {
        hash ^= (bits >> (8 * byte)) & 0xFFu;
        hash *= kFnvPrime;
    }
}

Result<CubeAxes> CubeAxes::Of(const ImageDescriptor& descriptor) {
    const auto& axes = descriptor.axes;
    const auto x = AxisIndex(axes, AxisRole::spatial_x);
    const auto y = AxisIndex(axes, AxisRole::spatial_y);
    const auto spectral = AxisIndex(axes, AxisRole::spectral);
    if (!x || !y || !spectral) {
        return Error{ErrorCode::invalid_argument,
                     "the bench reads cubes, and image " + descriptor.id + " has no " +
                         (!x   ? "spatial x"
                          : !y ? "spatial y"
                               : "spectral") +
                         " axis",
                     descriptor.id};
    }
    CubeAxes cube;
    cube.rank = axes.size();
    cube.x = *x;
    cube.y = *y;
    cube.spectral = *spectral;
    cube.polarization = AxisIndex(axes, AxisRole::polarization);
    cube.time = AxisIndex(axes, AxisRole::time);
    cube.width = axes[*x].length;
    cube.height = axes[*y].length;
    cube.channels = axes[*spectral].length;
    cube.polarizations = cube.polarization ? axes[*cube.polarization].length : 1;
    return cube;
}

ChunkBox ChunkBox::Spanning(const std::array<std::uint64_t, 4>& from, const std::array<std::uint64_t, 4>& to,
                            const CubeAxes& axes, const std::vector<std::uint64_t>& chunk_shape) {
    const std::array<std::optional<std::size_t>, 4> index{axes.x, axes.y, axes.spectral, axes.polarization};
    ChunkBox box;
    for (std::size_t axis = 0; axis < index.size(); ++axis) {
        std::uint64_t chunk = 1;
        if (index[axis] && *index[axis] < chunk_shape.size()) {
            chunk = std::max<std::uint64_t>(chunk_shape[*index[axis]], 1);
        }
        box.first[axis] = from[axis] / chunk;
        box.last[axis] = (std::max(to[axis], from[axis] + 1) - 1) / chunk;
    }
    return box;
}

std::uint64_t Fingerprint(const float* values, std::size_t count) {
    std::uint64_t hash = kFnvOffset;
    for (std::size_t index = 0; index < count; ++index) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &values[index], sizeof(bits));
        Mix(hash, std::isnan(values[index]) ? kQuietNaN32 : bits, sizeof(bits));
    }
    return hash;
}

std::uint64_t Fingerprint(const double* values, std::size_t count) {
    std::uint64_t hash = kFnvOffset;
    for (std::size_t index = 0; index < count; ++index) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &values[index], sizeof(bits));
        Mix(hash, std::isnan(values[index]) ? kQuietNaN64 : bits, sizeof(bits));
    }
    return hash;
}

}  // namespace carta::zarr::bench
