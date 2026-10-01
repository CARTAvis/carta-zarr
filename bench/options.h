/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_BENCH_OPTIONS_H_
#define CARTA_ZARR_BENCH_OPTIONS_H_

// What carta-zarr-bench was asked to do, parsed from its command line.

#include <carta-zarr/descriptor.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace carta::zarr::bench {

// The ways a CARTA user reads a cube, each measured on its own.
enum class Mode {
    plane,           // one whole l x m plane at a random channel
    spectrum,        // every channel at one random pixel
    region,          // ReduceSpectral over a box covering 5% of the plane
    cube_histogram,  // ComputeCubeHistogram over this process's share of the channels
    open,            // Context::Create, Dataset::Open and OpenImage, with nothing cached
};

const char* ModeName(Mode mode) noexcept;
std::optional<Mode> ParseMode(std::string_view name) noexcept;

// How many operations one process makes in one trial, when --ops does not say: as many as keep a
// trial of each mode about as long as the others on a cube of tens of gigabytes.
unsigned DefaultOps(Mode mode) noexcept;

// How the caches between a trial and the storage are emptied before it.
enum class ColdMethod {
    command,      // the administrator's --drop-cache-cmd, the only one that can reach a server's cache
    drop_caches,  // /proc/sys/vm/drop_caches, which needs root and Linux
    fadvise,      // posix_fadvise(DONTNEED) on every file of the dataset, which needs Linux
    off,          // nothing: a warm measurement
};

const char* ColdMethodName(ColdMethod method) noexcept;
std::optional<ColdMethod> ParseColdMethod(std::string_view name) noexcept;

struct RunOptions {
    std::string dataset;
    // Empty for the dataset's default image.
    std::string image_id;
    std::vector<Mode> modes;
    unsigned trials = 5;
    // Operations per process per trial; DefaultOps when unset.
    std::optional<unsigned> ops;
    std::uint64_t seed = 1;
    ContextOptions context;
    std::size_t read_budget_bytes = 0;
    unsigned processes = 1;
    // Unset for auto: the command when there is one, then drop_caches, then fadvise, then off.
    std::optional<ColdMethod> cold;
    std::string drop_cache_command;
    std::chrono::seconds trial_timeout{600};
    // Empty for stdout.
    std::string csv_path;
    std::string label;
    bool resume = false;

    unsigned OpsFor(Mode mode) const noexcept {
        return ops.value_or(DefaultOps(mode));
    }
};

struct ProbeOptions {
    std::string dataset;
    std::string image_id;
};

// Asked for help, or asked wrongly. `error` says which: help goes to stdout and exits 0.
struct Usage {
    std::string message;
    bool error = false;
};

using Command = std::variant<RunOptions, ProbeOptions, Usage>;

Command ParseCommandLine(int argc, const char* const* argv);

// A size in bytes, with an optional binary suffix: 4096, 512K, 64M, 4G, 1T, and the same with iB or
// B after them. Every multiple is a power of 1024, because what these size is memory.
std::optional<std::size_t> ParseSize(std::string_view text) noexcept;

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_OPTIONS_H_
