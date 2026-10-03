/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_BENCH_PLAN_H_
#define CARTA_ZARR_BENCH_PLAN_H_

// What every mode plans with, and nothing any one mode decides: the cube's axes, an operation, the
// generator positions are drawn from, the chunks an operation spans, and the fingerprint of what it
// read.
//
// Positions are a function of the seed, the mode, the trial and the cube's logical shape -- never of
// its layout -- so every layout of one cube is asked for the same pixels, and a difference in time is
// the layout's. They come from a generator written out here rather than from <random>, whose
// distributions are free to differ between standard libraries, so the Mac and a Linux server agree.

#include "options.h"

#include <carta-zarr/carta_zarr.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace carta::zarr::bench {

// The axes the operations address, as logical indices into ImageDescriptor::axes. An image with no
// polarization or time axis reads as one with a single plane of each.
struct CubeAxes {
    std::size_t rank = 0;
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t spectral = 0;
    std::optional<std::size_t> polarization;
    std::optional<std::size_t> time;

    std::uint64_t width = 1;
    std::uint64_t height = 1;
    std::uint64_t channels = 1;
    std::uint64_t polarizations = 1;

    // The image's axes, or an error naming the role it lacks.
    static Result<CubeAxes> Of(const ImageDescriptor& descriptor);
};

// One operation, as planned before the trial starts.
struct Operation {
    Mode mode = Mode::plane;
    std::uint64_t x = 0;
    std::uint64_t y = 0;
    std::uint64_t width = 0;
    std::uint64_t height = 0;
    std::uint64_t channel = 0;
    std::uint64_t channel_count = 0;
    std::uint64_t polarization = 0;
    // Whether another operation of this trial, in any process, reads the same position: there were
    // more operations than distinct positions to give them.
    bool overlap = false;
    // Whether it reads a chunk that an earlier operation of its own process read, or that any
    // operation of another process reads. Its time may then be the page cache's rather than the
    // storage's -- the earlier read brought the chunk in, or the other process is bringing it in
    // now -- so it is not a first touch however fresh its cache pool. See MarkSharedChunks.
    bool shares_chunks = false;

    // Where it reads, for the CSV: semicolon-separated, so that it stays one field.
    std::string Describe() const;
};

// Which process of which trial a plan is for, and of how many.
struct PlanSeed {
    std::uint64_t seed = 1;
    unsigned trial = 0;
    unsigned processes = 1;
    unsigned process_index = 0;
};

// SplitMix64: small, fast, and the same sequence everywhere.
class Stream {
public:
    Stream(std::uint64_t seed, Mode mode, unsigned trial, std::uint64_t salt)
        : _state(seed * 0x9E3779B97F4A7C15ull + (static_cast<std::uint64_t>(mode) << 48) +
                 (static_cast<std::uint64_t>(trial) << 16) + salt) {}

    std::uint64_t Next() {
        std::uint64_t z = (_state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }

    // Uniform in [0, bound), by the multiply-shift of Lemire without its rejection step: the bias is
    // bound / 2^64, far below anything a benchmark could notice.
    std::uint64_t Below(std::uint64_t bound) {
        if (bound <= 1) {
            return 0;
        }
        return static_cast<std::uint64_t>((static_cast<unsigned __int128>(Next()) * bound) >> 64);
    }

private:
    std::uint64_t _state;
};

struct Draw {
    std::uint64_t value = 0;
    bool overlap = false;
};

// `count` values from [0, pool), distinct while the pool lasts and repeating it after. Each prefix
// of the result is the same whatever `count` is, which is what keeps process 0's positions fixed as
// processes are added.
std::vector<Draw> DistinctSample(Stream& stream, std::uint64_t pool, std::uint64_t count);

// The operations of a positioned mode: one sequence of `pool` positions drawn for the whole trial,
// process p taking the p-th run of `ops` of them, so processes never share one while there are
// enough, and process 0 reads the same positions whatever the process count. `place` turns a draw
// into an operation, drawing whatever else it needs from `details`, and is called for every
// operation of the trial in order, other processes' included, so that each process's details are the
// same whatever process it is.
template <typename Place>
std::vector<Operation> PlanDraws(Mode mode, const PlanSeed& at, unsigned ops, std::uint64_t pool, Place place) {
    // Two streams: one for which position, one for everything about it. Kept apart so that adding a
    // draw to one never shifts the other.
    Stream positions(at.seed, mode, at.trial, 1);
    Stream details(at.seed, mode, at.trial, 2);
    const std::uint64_t total = static_cast<std::uint64_t>(at.processes) * ops;
    const std::uint64_t first = static_cast<std::uint64_t>(at.process_index) * ops;
    const auto draws = DistinctSample(positions, pool, total);

    std::vector<Operation> operations(ops);
    for (auto& operation : operations) {
        operation.mode = mode;
    }
    for (std::uint64_t index = 0; index < total; ++index) {
        Operation operation;
        operation.mode = mode;
        operation.overlap = draws[index].overlap;
        place(draws[index].value, details, operation);
        if (index >= first && index < first + ops) {
            operations[index - first] = operation;
        }
    }
    return operations;
}

// The chunks an operation reads, as a range of chunk indices along each axis it moves on: the spatial
// two, the spectral and the polarization. Every other axis is read at 0 by every operation.
struct ChunkBox {
    std::array<std::uint64_t, 4> first{};
    std::array<std::uint64_t, 4> last{};

    // The chunks that hold the pixels from `from` up to `to`, along x, y, spectral and polarization.
    static ChunkBox Spanning(const std::array<std::uint64_t, 4>& from, const std::array<std::uint64_t, 4>& to,
                             const CubeAxes& axes, const std::vector<std::uint64_t>& chunk_shape);

    bool Meets(const ChunkBox& other) const {
        for (std::size_t axis = 0; axis < first.size(); ++axis) {
            if (last[axis] < other.first[axis] || other.last[axis] < first[axis]) {
                return false;
            }
        }
        return true;
    }
};

// FNV-1a over a run of values, with every NaN hashed as the one quiet NaN.
inline constexpr std::uint64_t kFnvOffset = 0xCBF29CE484222325ull;
void Mix(std::uint64_t& hash, std::uint64_t bits, unsigned bytes);
std::uint64_t Fingerprint(const float* values, std::size_t count);
std::uint64_t Fingerprint(const double* values, std::size_t count);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_PLAN_H_
