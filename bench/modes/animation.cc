/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// animation: consecutive planes from a random channel at a frame rate, through the context's cache,
// which its frames share as the backend's do. Where plane and spectrum time a first touch, this is
// where reuse is the point, and what it measures.

#include "cube.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <optional>
#include <thread>
#include <vector>

#include <carta-zarr/read_ahead.h>

namespace carta::zarr::bench {

namespace {

class AnimationRunner final : public CubeRunner {
public:
    // How an animation's frames went, the first apart: it is a cold read whatever the layout, and what a
    // layout decides is how the frames after it go.
    struct FrameStats {
        // The frames played, which is the channels the animation read: as many as were asked for, or every
        // channel of a cube with fewer.
        unsigned frames = 0;
        double first_s = 0.0;
        // The read times of every frame after the first.
        double median_s = 0.0;
        double max_s = 0.0;
        // Frames after the first not ready by the end of their turn, and the most any was late by. Zero
        // when the frames are read back to back, which gives them no turn to miss.
        unsigned late = 0;
        double late_max_s = 0.0;
        // Prefetches of the next run of chunks started, and how many of them the animation caught up with
        // before they had finished: ReadAheadStats's prefetches and caught_up. Zero without prefetch, or
        // when the cache cannot hold two runs. Prefetches stop once a frame is late while one is under way,
        // so the first is the more telling.
        unsigned prefetches = 0;
        unsigned late_prefetches = 0;
    };

    AnimationRunner(Image image, CubeAxes axes, AnimationSettings settings)
        : CubeRunner(std::move(image), axes), _fps(settings.fps), _prefetch(settings.prefetch) {
        // Sized here so that no operation pays for growing it.
        _pixels.resize(axes.width * axes.height);
    }

    // Plane after plane, as carta-backend serves a playing animation: each frame its own read, through
    // the context's cache, so that a frame finds what the one before it decoded when they share chunks.
    //
    // Played at a frame rate, a frame is read no earlier than its turn and is late if it is not ready by
    // the end of it; a frame that runs over pushes the ones after it back, as a viewer waiting on it would
    // see.
    //
    // With prefetch, reading ahead is the library's ReadAhead, as carta-backend's is: told after every
    // frame when it began, whether it was late and which planes come next, it decodes the next run of
    // chunks on a thread of its own while this run's frames play from the cache -- one at a time, none once
    // a frame is late while one is under way, and none at all unless the cache the frames read through
    // holds two runs. A context left to TensorStore's cache holds nothing, so without --cache-bytes nothing
    // is read ahead, and the reason is said on stderr. A prefetch the animation catches up with is counted
    // as late, but is no reason to stop: the frame that caught it waits for the decode under way rather
    // than starting another. See ReadAhead and ADR 0016.
    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) override {
        _pixel_count = 0;
        _frames.reset();
        _frame_stats.reset();
        using Clock = ReadAhead::Clock;
        std::uint64_t hash = kFnvOffset;
        std::uint64_t elements = 0;
        Operation frame = operation;
        const auto& axes = this->axes();
        const auto period = _fps > 0.0
                              ? std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(1.0 / _fps))
                              : Clock::duration::zero();
        const auto plane_at = [&](std::uint64_t channel) {
            ReadRequest request;
            request.axes.assign(axes.rank, Range{0, 1, 1});
            request.axes[axes.x] = Range{0, axes.width, 1};
            request.axes[axes.y] = Range{0, axes.height, 1};
            request.axes[axes.spectral] = Range{channel, 1, 1};
            if (axes.polarization) {
                request.axes[*axes.polarization] = Range{operation.polarization, 1, 1};
            }
            return AnimatedPlane{0, std::move(request)};
        };

        std::optional<ReadAhead> reading;
        if (_prefetch) {
            auto made = ReadAhead::For({{image(), options}});
            if (made) {
                reading = std::move(made).value();
            } else {
                std::fprintf(stderr, "carta-zarr-bench: an animation reads nothing ahead: %s\n",
                             made.error().message.c_str());
            }
        }

        FrameStats stats;
        stats.frames = static_cast<unsigned>(operation.channel_count);
        std::vector<double> reads;
        const auto start = Clock::now();
        auto turn = start;
        for (std::uint64_t index = 0; index < operation.channel_count; ++index) {
            frame.channel = operation.channel + index;
            if (period > Clock::duration::zero()) {
                std::this_thread::sleep_until(turn);
            }
            const auto began = Clock::now();
            auto read = ReadPixels(image(), axes, frame, true, _pixels.data(), options);
            const auto done = Clock::now();
            if (!read) {
                return std::move(read).error();
            }
            _pixel_count = *read;
            elements += *read;
            const bool late = period > Clock::duration::zero() && done > turn + period;
            if (reading) {
                std::vector<std::vector<AnimatedPlane>> upcoming;
                const auto left = operation.channel_count - index - 1;
                const auto upcoming_frames = std::min<std::uint64_t>(left, ReadAhead::kUpcomingFrames);
                for (std::uint64_t ahead = 1; ahead <= upcoming_frames; ++ahead) {
                    upcoming.push_back({plane_at(frame.channel + ahead)});
                }
                reading->Served(began, late, {plane_at(frame.channel)}, upcoming);
            }
            // Played at a frame rate, only the last frame is fingerprinted: a plane of a large cube takes
            // longer to hash than a frame's turn, and hashing every one would make every frame late.
            if (period == Clock::duration::zero() || index + 1 == operation.channel_count) {
                Mix(hash, bench::Fingerprint(_pixels.data(), _pixel_count), sizeof(hash));
            }

            const double seconds = std::chrono::duration<double>(done - began).count();
            if (index == 0) {
                stats.first_s = seconds;
            } else {
                reads.push_back(seconds);
                if (late) {
                    ++stats.late;
                    stats.late_max_s =
                        std::max(stats.late_max_s, std::chrono::duration<double>(done - turn - period).count());
                }
            }
            // The next frame's turn follows this one's, or this frame's end when it ran over.
            turn = std::max(turn + period, done);
        }
        if (reading) {
            const auto ahead = reading->stats();
            stats.prefetches = ahead.prefetches;
            stats.late_prefetches = ahead.caught_up;
        }
        if (!reads.empty()) {
            std::sort(reads.begin(), reads.end());
            stats.median_s = reads[reads.size() / 2];
            stats.max_s = reads.back();
        }
        _frame_stats = stats;
        _frames = hash;
        return elements;
    }

    // Its frames', each fingerprinted as it is read and the fingerprints hashed together.
    std::uint64_t Fingerprint() const override {
        return _frames.value_or(bench::Fingerprint(_pixels.data(), _pixel_count));
    }

    // How its frames went, when it played them.
    void Record(Row& result) const override {
        CubeRunner::Record(result);
        if (const auto& frames = _frame_stats) {
            result.frames_played = frames->frames;
            result.frame_first_s = frames->first_s;
            result.frame_median_s = frames->median_s;
            result.frame_max_s = frames->max_s;
            result.late_frames = frames->late;
            result.late_max_s = frames->late_max_s;
            result.prefetches = frames->prefetches;
            result.late_prefetches = frames->late_prefetches;
        }
    }

private:
    double _fps;
    bool _prefetch;
    std::vector<float> _pixels;
    std::size_t _pixel_count = 0;
    std::optional<std::uint64_t> _frames;
    std::optional<FrameStats> _frame_stats;
};

class Animation final : public Workload {
public:
    Animation(const RunOptions& options, AnimationSettings settings)
        : Workload(Mode::animation, options), _settings(settings) {}

    void WriteSettings(Row& row) const override {
        row.animation_frames = std::to_string(_settings.frames);
        std::array<char, 32> fps{};
        std::snprintf(fps.data(), fps.size(), "%g", _settings.fps);
        row.animation_fps = fps.data();
        row.animation_prefetch = _settings.prefetch ? "true" : "false";
    }

    // An animation reads `frames` consecutive channels from a start in a cell of channels of its own,
    // so that no two of a trial's animations play the same planes while there are enough channels for
    // them.
    std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const override {
        const auto frames = std::clamp<std::uint64_t>(_settings.frames, 1, std::max<std::uint64_t>(axes.channels, 1));
        const auto runs = std::max<std::uint64_t>(axes.channels / frames, 1);
        const auto run_cell = axes.channels / runs;
        return PlanDraws(mode(), at, ops, runs * axes.polarizations,
                         [&](std::uint64_t value, Stream& details, Operation& operation) {
                             operation.channel_count = frames;
                             operation.channel = (value % runs) * run_cell + details.Below(run_cell - frames + 1);
                             operation.polarization = value / runs;
                         });
    }

    std::string Describe(const Operation& operation) const override { return DescribeChannels(operation); }

    // Not a first touch, so ChunksRead says nothing: it reuses chunks on purpose.

    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        (void)context;
        auto axes = CubeAxes::Of(image.descriptor());
        if (!axes) {
            return std::move(axes).error();
        }
        return std::unique_ptr<Runner>(std::make_unique<AnimationRunner>(std::move(image), *axes, _settings));
    }

private:
    AnimationSettings _settings;
};

}  // namespace

std::unique_ptr<Workload> AnimationWorkload(const RunOptions& options, AnimationSettings animation) {
    return std::make_unique<Animation>(options, animation);
}

}  // namespace carta::zarr::bench
