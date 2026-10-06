/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_BENCH_MODE_H_
#define CARTA_ZARR_BENCH_MODE_H_

// The ways a CARTA user reads a cube, each everything about itself in one place: the settings that
// shape it and the columns they are written to, where its operations read, which chunks one spans,
// what a process has to have before the clock starts, the reading itself, and what is written down
// after it.
//
// Each mode used to be a case in a switch, in five files: its name and its default count of
// operations, its positions, how a position is written, which chunks it spans, whether it reads
// through a cache pool of its own, how it reads, its settings columns, and what a trial opens for it.
// Adding a column to animations touched seven files. A mode is now one file under modes/ and a line
// in Workload::For, and its name and default count in options.cc, which the command line needs before
// there is a Workload to ask.

#include "options.h"
#include "plan.h"
#include "record.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <carta-zarr/carta_zarr.h>

namespace carta::zarr::bench {

// One process's operations of one trial, run one at a time.
class Runner {
public:
    Runner() = default;
    Runner(const Runner&) = delete;
    Runner& operator=(const Runner&) = delete;
    Runner(Runner&&) = delete;
    Runner& operator=(Runner&&) = delete;
    virtual ~Runner() = default;

    // What an operation needs that is not part of what it measures, done before the clock starts.
    virtual Result<void> Prepare(const Operation& operation) {
        (void)operation;
        return {};
    }

    // The operation, timed: how many elements of the cube it covered.
    virtual Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) = 0;

    // A fingerprint of what the last operation read, taken after the clock stops, and only of what
    // every layout of the same pixels must agree on exactly: the pixels a read returns, with every NaN
    // one NaN; and for a reduction or a histogram the pixel counts and extremes, which are exact, but
    // not sums, whose rounding depends on the order the chunks were visited in.
    virtual std::uint64_t Fingerprint() const = 0;

    // The columns of `result` only this mode fills, whatever came of the operation: the bytes the
    // elements it covered are stored in, and anything else it has to say. After the clock stops.
    virtual void Record(Row& result) const = 0;
};

// One process's plan for one trial, and what runs it.
struct Planned {
    std::vector<Operation> plan;
    std::unique_ptr<Runner> runner;
};

class Workload {
public:
    // The one place every mode is named, each handed its own settings and the ones every mode reads.
    static std::unique_ptr<Workload> For(Mode mode, const RunOptions& options);

    Workload(const Workload&) = delete;
    Workload& operator=(const Workload&) = delete;
    Workload(Workload&&) = delete;
    Workload& operator=(Workload&&) = delete;
    virtual ~Workload() = default;

    Mode mode() const noexcept { return _mode; }

    // The settings columns only this mode has, which the run key is made from as well. Every other
    // mode leaves them empty.
    virtual void WriteSettings(Row& row) const { (void)row; }

    // The operations one process makes in one trial, a function of the seed, the trial and the cube's
    // logical shape and never of its layout. See PlanDraws.
    virtual std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const = 0;

    // Where an operation reads, for the CSV: semicolon-separated, so that it stays one field.
    virtual std::string Describe(const Operation& operation) const = 0;

    // The chunks an operation reads, given the image's chunk shape in logical axis order, or nothing
    // for a mode that is not timing a first touch. See MarkSharedChunks.
    //
    // Chunks and not shards: of two operations reading different chunks of one shard, the later finds
    // the shard's index in the page cache, sixteen bytes a chunk, and its chunks wherever storage left
    // them. Whether reading ahead in the shard's file brought them in is the storage's, and not
    // something the bench can know; counting by shard would leave a layout sharded across every
    // channel with one first touch a trial to rank it by.
    virtual std::optional<ChunkBox> ChunksRead(const Operation& operation, const CubeAxes& axes,
                                               const std::vector<std::uint64_t>& chunk_shape) const {
        (void)operation;
        (void)axes;
        (void)chunk_shape;
        return std::nullopt;
    }

    // Everything a process has to have before it says it is ready: for a mode that reads a cube, the
    // image opened through a context of its own, `row`'s image columns filled in, every process's plan
    // drawn so that this one's can be told which chunks another reads, and a runner over the image.
    // An error's message is the whole of what the row says.
    virtual Result<Planned> Begin(unsigned trial, unsigned process_index, Row& row) const;

    // A runner over an image already open, for a mode that reads a cube. The context is the one the
    // image was opened through, and the one a cache pool of the runner's own is made from.
    virtual Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const = 0;

protected:
    Workload(Mode mode, const RunOptions& options) : _mode(mode), _options(options) {}

    const RunOptions& options() const noexcept { return _options; }

private:
    Mode _mode;
    // Every mode's settings are in it, but a mode is handed its own as well, and reads those.
    RunOptions _options;
};

// Sets shares_chunks on every operation of a trial, given each process's plan, indexed by process,
// and the image's chunk shape in logical axis order. An operation that reads a chunk first is never
// marked: of two that share one, the earlier in its process is the first touch, unless the other is
// another process's, which may be reading it at the same moment. A mode whose ChunksRead says nothing
// is left alone.
void MarkSharedChunks(const Workload& workload, std::vector<std::vector<Operation>>& plans, const CubeAxes& axes,
                      const std::vector<std::uint64_t>& chunk_shape);

// The modes, each in modes/, made by Workload::For.
std::unique_ptr<Workload> PlaneWorkload(const RunOptions& options);
std::unique_ptr<Workload> SpectrumWorkload(const RunOptions& options);
std::unique_ptr<Workload> RegionWorkload(const RunOptions& options, RegionSettings region);
std::unique_ptr<Workload> CubeHistogramWorkload(const RunOptions& options, HistogramMethod histogram);
std::unique_ptr<Workload> OpenWorkload(const RunOptions& options);
std::unique_ptr<Workload> AnimationWorkload(const RunOptions& options, AnimationSettings animation);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_MODE_H_
