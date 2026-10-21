/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_READ_AHEAD_H_
#define CARTA_ZARR_SRC_READ_AHEAD_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"
#include "carta-zarr/read_ahead.h"
#include "carta-zarr/result.h"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace carta::zarr::internal {

// The chunks a read decodes, as the first and the last chunk it touches along each axis. Two reads
// with the same box decode the same chunks, which is what makes them one run; see GLOSSARY.md.
//
// Every axis counts alike. A plane is one element along the spectrum and along the polarization, and
// the chunk it lies in along each decides its run: a chunk deep in Stokes puts several polarizations
// in one run, as a chunk deep in frequency puts several channels in one.
struct Run {
    std::vector<std::uint64_t> first;
    std::vector<std::uint64_t> last;

    bool operator==(const Run& other) const { return first == other.first && last == other.last; }
    bool operator!=(const Run& other) const { return !(*this == other); }
};

// The run a read is in. An axis with no chunk shape to speak of is one chunk.
Run RunOf(const ChunkGeometry& geometry, const ReadRequest& request);

// The read that covers the run `plane` is in: along every axis, the whole of the chunks it touches,
// clipped to the axis. A prefetch asks for this rather than the plane, because decoding the plane
// alone decodes every pixel chunk of its run but the flag of only that plane, and a flag chunked
// shallower than the pixels keeps the run's other planes in chunks of their own.
ReadRequest RunRequest(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const ReadRequest& plane);

// What one run of a whole plane of this image holds once decoded: the plane rounded out to whole
// chunks along its two spatial axes, one chunk deep along every other axis, at the decoded size. When
// the read applies the pixel mask, add the flag chunks the run crosses, in the flag's own chunks at a
// byte per element. A flag may be chunked differently from its image; one stored in a single chunk
// counts whole for every run. `flag_geometry` without a chunk shape is taken to be the pixels'. Zero
// for an image without two spatial axes, which has no plane to animate.
std::uint64_t PlaneRunBytes(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                            const ChunkGeometry& flag_geometry, bool apply_mask);

// The cache an image's reads keep their chunks in, and how much it holds. Images reading through the
// same cache share it, so it has to be large enough for all of them.
struct CacheShare {
    // Equal for two images whose reads go through the same cache.
    const void* identity = nullptr;
    std::uint64_t bytes = 0;
};

// What reading ahead asks of one image: which run a plane is in, what a run holds and what holds it,
// and the decode itself.
//
// The seam between an image's real answers and a test's. A test's prefetch can block until the test
// releases it, which is the only way to test "a frame began while one was under way" without racing
// a real decode.
class RunSource {
public:
    RunSource() = default;
    RunSource(const RunSource&) = delete;
    RunSource& operator=(const RunSource&) = delete;
    RunSource(RunSource&&) = delete;
    RunSource& operator=(RunSource&&) = delete;
    virtual ~RunSource() = default;

    virtual Run RunOf(const ReadRequest& plane) const = 0;
    virtual std::uint64_t PlaneRunBytes() const = 0;
    virtual CacheShare Cache() const = 0;
    // Named in what is said when reading ahead is declined.
    virtual std::string Name() const = 0;
    // Decodes the run `plane` is in, so that a read of any plane of it finds its chunks decoded.
    // Blocks until it has, or `cancelled` says to stop, or it fails; says whether it finished.
    virtual bool Prefetch(const ReadRequest& plane, const std::function<bool()>& cancelled) const = 0;
};

// Reading ahead of an animation over the images `sources` answer for. What ReadAhead is, behind its
// handle; see ReadAhead and ADR 0016 for the policy.
class ReadingAhead {
public:
    using Clock = std::chrono::steady_clock;

    static constexpr std::size_t kUpcomingFrames = ReadAhead::kUpcomingFrames;

    // Reading ahead over `sources`. Fails with buffer_too_small, naming the sources, when a cache
    // cannot hold two runs of every source that reads through it, and with invalid_argument when
    // there are none or one cannot say what a run of it holds.
    static Result<std::unique_ptr<ReadingAhead>> For(std::vector<std::shared_ptr<const RunSource>> sources);

    // Cancels what is under way and waits for it.
    ~ReadingAhead();
    ReadingAhead(const ReadingAhead&) = delete;
    ReadingAhead& operator=(const ReadingAhead&) = delete;
    ReadingAhead(ReadingAhead&&) = delete;
    ReadingAhead& operator=(ReadingAhead&&) = delete;

    // Called after a frame that began at `began` and showed `shown`, with whether it was late by the
    // caller's measure. `upcoming` is the planes of the frames to come, nearest first; only the
    // first kUpcomingFrames are looked at.
    //
    // A late frame that began while a prefetch was under way stops all reading ahead. Otherwise, if
    // none is under way, it starts one that decodes, for each image shown, the first run in
    // `upcoming` after the one shown now, unless that run was already decoded ahead.
    void Served(Clock::time_point began, bool late, const std::vector<AnimatedPlane>& shown,
                const std::vector<std::vector<AnimatedPlane>>& upcoming);

    // Stops what is under way without waiting for it, and starts nothing more.
    void Cancel();

    ReadAheadStats Stats() const;
    bool UnderWay() const;

private:
    explicit ReadingAhead(std::vector<std::shared_ptr<const RunSource>> sources);
    void Join();

    std::vector<std::shared_ptr<const RunSource>> _sources;
    // Only Served, Cancel and the destructor touch the worker thread, and a caller calls those from
    // the one thread that plays the animation, or after it has stopped. The rest is shared with the
    // worker.
    std::thread _worker;
    std::atomic<bool> _cancelled{false};

    mutable std::mutex _mutex;
    bool _under_way = false;
    Clock::time_point _started;
    Clock::time_point _finished;
    // The runs the prefetch under way, or the last one, is decoding: by image.
    std::map<std::size_t, Run> _prefetching;
    bool _counted_catch = false;
    // The run last decoded ahead for each image, so that it is not decoded again.
    std::map<std::size_t, Run> _requested;
    ReadAheadStats _stats;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_READ_AHEAD_H_
