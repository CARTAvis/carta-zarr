/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_READ_AHEAD_H_
#define CARTA_ZARR_READ_AHEAD_H_

#include "carta-zarr/carta_zarr.h"
#include "carta-zarr/export.h"
#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include <chrono>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace carta::zarr {

// One plane a frame of an animation shows: which of the images its read-ahead was made for, by
// position, and the read of the plane.
struct AnimatedPlane {
    std::size_t image = 0;
    ReadRequest request;
};

// What reading ahead has done so far.
struct ReadAheadStats {
    // Runs decoded ahead, or started to be: one for each image a prefetch was for.
    unsigned prefetches = 0;
    // Of those, the ones a frame reached before they had finished. That frame waited for the decode
    // under way rather than starting its own, so such a prefetch is not wasted: it only started too
    // late to hide all of its run.
    unsigned caught_up = 0;
    // Whether reading ahead has stopped, because a frame was late while a prefetch was under way or
    // because it was cancelled. It does not start again.
    bool stopped = false;
    // Whether a prefetch is being decoded now.
    bool under_way = false;
};

// Reading ahead of a playing animation: decoding the next run of chunks of each animated image while
// the frames of this run play from the cache, so that the frame entering it does not stall. See
// GLOSSARY.md, Run and Read-ahead, and ADR 0016.
//
// A frame that enters a run decodes all of it, which costs what jumping to the channel does, and every
// frame after it in the run is served from the cache in a few milliseconds. On Lustre, a 7763 x 4742
// cube in 512 x 512 x 16 chunks stalled for up to 400 ms every sixteenth frame; decoding the next run
// in the time the frames of this one left over hid every stall for one viewer.
//
// It reads ahead only while there is spare time. One prefetch runs at a time, and once a frame is
// late while one runs, reading ahead stops: with eight viewers animating that cube at once,
// prefetches that kept ahead of every run still doubled the time of the frames played from the cache.
// A frame that reaches a run still being decoded waits for that decode rather than starting its own;
// see Image::Prefetch.
//
// A handle like Image: copies share it, and the last of them to go cancels the prefetch under way and
// waits for it. Nothing it starts outlives it.
//
// Threads: call Served from one thread at a time, the one that serves the frames, even through two
// copies. Cancel and stats may be called from any thread at the same time. A prefetch runs on a
// thread of its own and calls each image's cancellation_requested from there, so that callable must be
// safe to call from another thread.
class CARTA_ZARR_EXPORT ReadAhead final {
public:
    using Clock = std::chrono::steady_clock;

    // How many frames ahead the next run is looked for, and so how many upcoming frames a caller
    // needs to pass. A run deeper than this is decoded ahead once it is this near, which at CARTA's 5
    // frames a second is 13 s before it is needed. Measured, not configured: ADR 0014.
    static constexpr std::size_t kUpcomingFrames = 64;

    ReadAhead(const ReadAhead&) = default;
    ReadAhead& operator=(const ReadAhead&) = default;
    // A copy, so that the handle moved from still refers to what it did. See ADR 0013.
    ReadAhead(ReadAhead&& other) noexcept : ReadAhead(other) {}
    ReadAhead& operator=(ReadAhead&& other) noexcept { return *this = other; }
    ~ReadAhead();

    // Reading ahead of the animated `images`, each read with its own options. Runs are decoded into
    // the cache that the image's frames are read through. AnimatedPlane::image is an index into
    // `images`.
    //
    // Fails with buffer_too_small, naming the images and the size needed, unless each cache holds two
    // runs of every image that reads through it: the one playing and the one decoded ahead. Otherwise
    // the second would evict the first. A run of 512 x 512 x 4 chunks of the cube above is 589 MB; of
    // 512 x 512 x 16, 2.4 GB. Fails with invalid_argument when there are no images, or one has no plane
    // to animate.
    static Result<ReadAhead> For(const std::vector<std::pair<Image, ReadOptions>>& images);

    // Call after every frame with when it began, whether it was late by the caller's own measure, the
    // planes it showed, and the planes of the frames to come, nearest first. Only the first
    // kUpcomingFrames of those are looked at. Returns at once; a prefetch it starts runs on a thread
    // of its own.
    //
    // A late frame that began while a prefetch was under way stops all reading ahead. Otherwise, if
    // none is under way, it starts one that decodes, for each image shown, the next run after the one
    // shown now, unless that run was already decoded ahead. Never fails: reading ahead that cannot be
    // done is skipped.
    void Served(Clock::time_point began, bool late, const std::vector<AnimatedPlane>& shown,
                const std::vector<std::vector<AnimatedPlane>>& upcoming);

    // Stops the prefetch under way without waiting for it, and starts nothing more. For an animation
    // that has stopped, or whose image is closing.
    void Cancel();

    ReadAheadStats stats() const;

private:
    class Impl;
    explicit ReadAhead(std::shared_ptr<Impl> impl);

    std::shared_ptr<Impl> _impl;
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_READ_AHEAD_H_
