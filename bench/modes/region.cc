/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// region: ReduceSpectral over a box covering a share of the plane, every channel of it, timing when
// its results arrive and, when asked to, how long it takes to stop.

#include "cube.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace carta::zarr::bench {

namespace {

// How many boxes of `side` (a share of the axis) fit along an axis of `length` pixels: at 5% of the
// plane a side is 0.2236 of the axis, so four, each in a cell of 0.25 with room to move.
std::uint64_t GridCells(std::uint64_t length, double side) {
    const auto cells = static_cast<std::uint64_t>(std::floor(1.0 / side));
    return std::clamp<std::uint64_t>(cells, 1, std::max<std::uint64_t>(length, 1));
}

// The chunks of `length` pixels from `start`, in chunks of `chunk`.
std::uint64_t ChunksSpanned(std::uint64_t start, std::uint64_t length, std::uint64_t chunk) {
    return length == 0 ? 0 : (start + length - 1) / chunk - start / chunk + 1;
}

// The shortest spelling of `value` that reads back as it. Settings columns are part of the run key,
// and at four decimals 0.000001 and 0.000049 were both 0.0000 -- region boxes of 4 and 36 elements on
// a 1024-square plane -- and a resume skipped the second as done. max_digits10 always reads back; the
// shortest that does keeps 0.05 from being written 0.050000000000000003.
std::string Shortest(double value) {
    std::array<char, 32> text{};
    for (int digits = 1; digits <= std::numeric_limits<double>::max_digits10; ++digits) {
        std::snprintf(text.data(), text.size(), "%.*g", digits, value);
        if (std::strtod(text.data(), nullptr) == value) {
            break;
        }
    }
    return text.data();
}

using Clock = std::chrono::steady_clock;

// When a reduction's sink was called, as the columns of the same names say: every call, and the
// finished blocks among them. A gap is between two arrivals, never from the start, which is what
// block_first_s says, so that a slow first block does not hide how often the rest came.
class HandOvers {
public:
    void Start(Clock::time_point start) {
        *this = HandOvers();
        _start = start;
    }

    void Arrived(bool complete) {
        const auto now = Clock::now();
        if (_handovers > 0) {
            _handover_gap_max_s = std::max(_handover_gap_max_s, Seconds(now - _last_handover));
        }
        ++_handovers;
        _last_handover = now;
        if (complete) {
            if (_blocks == 0) {
                _block_first_s = Seconds(now - _start);
            } else {
                _block_gap_max_s = std::max(_block_gap_max_s, Seconds(now - _last_block));
            }
            ++_blocks;
            _last_block = now;
        }
    }

    void Record(Row& result) const {
        result.handovers = _handovers;
        result.blocks = _blocks;
        if (_blocks > 0) {
            result.block_first_s = _block_first_s;
        }
        if (_blocks > 1) {
            result.block_gap_max_s = _block_gap_max_s;
        }
        if (_handovers > 1) {
            result.handover_gap_max_s = _handover_gap_max_s;
        }
    }

private:
    static double Seconds(Clock::duration span) { return std::chrono::duration<double>(span).count(); }

    Clock::time_point _start;
    Clock::time_point _last_handover;
    Clock::time_point _last_block;
    unsigned _handovers = 0;
    unsigned _blocks = 0;
    double _block_first_s = 0.0;
    double _handover_gap_max_s = 0.0;
    double _block_gap_max_s = 0.0;
};

class RegionRunner final : public CubeRunner {
public:
    RegionRunner(Image image, CubeAxes axes, RegionSettings settings)
        : CubeRunner(std::move(image), axes), _settings(settings) {
        const auto& chunks = this->image().chunk_geometry().chunk_shape;
        const auto along = [&](std::size_t axis) {
            return axis < chunks.size() ? std::max<std::uint64_t>(chunks[axis], 1) : std::uint64_t{1};
        };
        _chunk_x = along(axes.x);
        _chunk_y = along(axes.y);
        _chunk_depth = along(axes.spectral);
        // Sized here so that no operation pays for growing it.
        _exact.reserve(axes.channels * 4);
    }

    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) override {
        const auto start = Clock::now();
        _hand_overs.Start(start);
        _stop_s.reset();
        _exact.clear();
        const RegionMask region{operation.x, operation.y, operation.width, operation.height, {}};
        SpectralReduceRequest request;
        request.planes.spectral = Range{0, axes().channels, 1};
        request.emit_every_channels = EmitHint(operation);
        request.planes.polarization = operation.polarization;
        request.regions = {&region, 1};
        // What carta-backend asks of a region's spectral profile, sum_sq_dev included: sigma is made
        // from it, and it is one more accumulator per channel and one more pass of arithmetic per
        // pixel, so a run without it timed a lighter reduction than the one being stood in for.
        request.statistics = Statistic::num_pixels | Statistic::nan_count | Statistic::sum | Statistic::sum_sq |
                             Statistic::min | Statistic::max | Statistic::sum_sq_dev;

        // Told to stop as a backend tells it: the flag its storage operations check, and the sink.
        const bool stops = _settings.stop_after_s > 0.0;
        // In seconds rather than as a time point, which a stop far enough off would overflow.
        const auto elapsed = [&](Clock::time_point now) { return std::chrono::duration<double>(now - start).count(); };
        const auto stopped = [&] { return stops && elapsed(Clock::now()) >= _settings.stop_after_s; };
        ReadOptions told = options;
        if (stops) {
            told.control.cancellation_requested = [&, before = options.control.cancellation_requested] {
                return stopped() || (before && before());
            };
        }

        _exact.assign(axes().channels * 4, std::numeric_limits<double>::quiet_NaN());
        const auto sink = [&](const SpectralBlock& block) {
            _hand_overs.Arrived(block.complete);
            if (!block.complete) {
                return !stopped();
            }
            for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
                const auto totals = block.Totals(0, channel);
                auto* exact = &_exact[(block.first_channel + channel) * 4];
                exact[0] = totals.num_pixels;
                exact[1] = totals.nan_count;
                exact[2] = totals.min;
                exact[3] = totals.max;
            }
            // A stop that comes with the last block comes too late to stop anything, but it still waited
            // for that block, which stop_s counts below.
            return block.first_channel + block.channel_count == axes().channels || !stopped();
        };
        auto reduced = image().ReduceSpectral(request, sink, told);
        const auto end = Clock::now();
        // Judged at the return, not by the clock now, which may have passed the stop since.
        const bool past_stop = stops && elapsed(end) >= _settings.stop_after_s;
        if (!reduced) {
            _exact.clear();
            // A stop this runner asked for is what it measures, not an error. The trial's deadline
            // still is one.
            const bool deadline = end >= options.control.deadline;
            if (reduced.error().code == ErrorCode::cancelled && past_stop && !deadline) {
                _stop_s = elapsed(end) - _settings.stop_after_s;
                return std::uint64_t{0};
            }
            return std::move(reduced).error();
        }
        // Told to stop while its last slab was read or accumulated, it finished instead, and the stop
        // waited for that as long as a cancelled one would have.
        if (past_stop) {
            _stop_s = elapsed(end) - _settings.stop_after_s;
        }
        return operation.width * operation.height * axes().channels;
    }

    // Each channel's pixel count, NaN count and extremes, which are exact; not its sums. A stopped
    // operation fingerprints as nothing read, since how far it got depends on the layout.
    std::uint64_t Fingerprint() const override { return bench::Fingerprint(_exact.data(), _exact.size()); }

    void Record(Row& result) const override {
        CubeRunner::Record(result);
        _hand_overs.Record(result);
        result.stop_s = _stop_s;
    }

private:
    // The emit_every_channels that asks for blocks of about _settings.emit_chunks pixel chunks, or 0.
    std::uint32_t EmitHint(const Operation& operation) const {
        if (_settings.emit_chunks == 0) {
            return 0;
        }
        const std::uint64_t footprint =
            std::max<std::uint64_t>(1, ChunksSpanned(operation.x, operation.width, _chunk_x) *
                                           ChunksSpanned(operation.y, operation.height, _chunk_y));
        const std::uint64_t layers = std::max<std::uint64_t>(1, _settings.emit_chunks / footprint);
        return static_cast<std::uint32_t>(
            std::min<std::uint64_t>(layers * _chunk_depth, std::numeric_limits<std::uint32_t>::max()));
    }

    RegionSettings _settings;
    std::uint64_t _chunk_x = 1;
    std::uint64_t _chunk_y = 1;
    std::uint64_t _chunk_depth = 1;
    std::vector<double> _exact;
    HandOvers _hand_overs;
    std::optional<double> _stop_s;
};

class Region final : public Workload {
public:
    Region(const RunOptions& options, RegionSettings settings) : Workload(Mode::region, options), _settings(settings) {}

    void WriteSettings(Row& row) const override {
        row.region_fraction = Shortest(_settings.fraction);
        row.region_emit_chunks = std::to_string(_settings.emit_chunks);
        row.region_stop_after_s = Shortest(_settings.stop_after_s);
    }

    // A box covers the fraction of the plane, square in pixels' proportion to the plane, and sits in a
    // cell of a grid of as many such boxes as fit along each side, so that no two of a trial's boxes
    // overlap while the grid has cells enough.
    std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const override {
        const double side = std::sqrt(_settings.fraction);
        const auto cells_x = GridCells(axes.width, side);
        const auto cells_y = GridCells(axes.height, side);
        const auto cell_width = axes.width / cells_x;
        const auto cell_height = axes.height / cells_y;
        const auto box_width =
            std::clamp<std::uint64_t>(std::llround(side * static_cast<double>(axes.width)), 1, cell_width);
        const auto box_height =
            std::clamp<std::uint64_t>(std::llround(side * static_cast<double>(axes.height)), 1, cell_height);
        return PlanDraws(
            mode(), at, ops, cells_x * cells_y, [&](std::uint64_t value, Stream& details, Operation& operation) {
                operation.width = box_width;
                operation.height = box_height;
                operation.x = (value % cells_x) * cell_width + details.Below(cell_width - box_width + 1);
                operation.y = (value / cells_x) * cell_height + details.Below(cell_height - box_height + 1);
                operation.polarization = details.Below(axes.polarizations);
            });
    }

    std::string Describe(const Operation& operation) const override {
        return "pol=" + std::to_string(operation.polarization) + ";l=" + std::to_string(operation.x) + ":" +
               std::to_string(operation.x + operation.width) + ";m=" + std::to_string(operation.y) + ":" +
               std::to_string(operation.y + operation.height);
    }

    std::optional<ChunkBox> ChunksRead(const Operation& operation, const CubeAxes& axes,
                                       const std::vector<std::uint64_t>& chunk_shape) const override {
        return ChunkBox::Spanning(
            {operation.x, operation.y, 0, operation.polarization},
            {operation.x + operation.width, operation.y + operation.height, axes.channels, operation.polarization + 1},
            axes, chunk_shape);
    }

    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        (void)context;
        auto axes = CubeAxes::Of(image.descriptor());
        if (!axes) {
            return std::move(axes).error();
        }
        return std::unique_ptr<Runner>(std::make_unique<RegionRunner>(std::move(image), *axes, _settings));
    }

private:
    RegionSettings _settings;
};

}  // namespace

std::unique_ptr<Workload> RegionWorkload(const RunOptions& options, RegionSettings region) {
    return std::make_unique<Region>(options, region);
}

}  // namespace carta::zarr::bench
