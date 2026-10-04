/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "read_ahead.h"

#include "chunk_blocks.h"

#include <algorithm>
#include <map>
#include <utility>

namespace carta::zarr::internal {
namespace {

// A size as a message says it: in MiB once there are any, so that a cache of gigabytes reads as one.
std::string Size(std::uint64_t bytes) {
    return bytes < (1ULL << 20) ? std::to_string(bytes) + " bytes" : std::to_string(bytes >> 20) + " MiB";
}

// Chunks along one axis of `length` elements: one for an axis with no chunk shape to speak of.
std::uint64_t ChunksAlong(std::uint64_t length, std::uint64_t chunk) {
    return chunk == 0 ? 1 : std::max<std::uint64_t>(1, (length + chunk - 1) / chunk);
}

// The most flag chunks one pixel chunk along an axis crosses: one when the two are chunked alike, more
// when the flag's are shorter or their boundaries fall inside a pixel chunk. Asked of every pixel
// chunk rather than bounded, since a bound that assumes the worst alignment doubles the count of two
// chunk shapes that line up; an axis is at most a few thousand chunks long.
std::uint64_t FlagChunksAcrossARun(std::uint64_t length, std::uint64_t pixel_chunk, std::uint64_t flag_chunk) {
    if (flag_chunk == 0 || length == 0) {
        return 1;
    }
    const std::uint64_t step = pixel_chunk == 0 ? length : pixel_chunk;
    std::uint64_t most = 1;
    for (std::uint64_t first = 0; first < length; first += step) {
        const std::uint64_t last = std::min(first + step, length) - 1;
        most = std::max(most, (last / flag_chunk) - (first / flag_chunk) + 1);
    }
    return most;
}

}  // namespace

Run RunOf(const ChunkGeometry& geometry, const ReadRequest& request) {
    Run run;
    run.first.reserve(request.axes.size());
    run.last.reserve(request.axes.size());
    for (std::size_t axis = 0; axis < request.axes.size(); ++axis) {
        const auto chunk = axis < geometry.chunk_shape.size() ? geometry.chunk_shape[axis] : 0;
        const auto& range = request.axes[axis];
        const auto step = std::max<std::uint64_t>(1, range.stride);
        const auto last = range.start + (range.count == 0 ? 0 : (range.count - 1) * step);
        run.first.push_back(chunk == 0 ? 0 : range.start / chunk);
        run.last.push_back(chunk == 0 ? 0 : last / chunk);
    }
    return run;
}

std::uint64_t PlaneRunBytes(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                            const ChunkGeometry& flag_geometry, bool apply_mask) {
    const auto x = AxisIndex(descriptor.axes, AxisRole::spatial_x);
    const auto y = AxisIndex(descriptor.axes, AxisRole::spatial_y);
    if (!x || !y) {
        return 0;
    }
    const auto chunk = [](const ChunkGeometry& of, std::size_t axis) {
        return axis < of.chunk_shape.size() ? of.chunk_shape[axis] : 0;
    };
    const std::uint64_t chunks = ChunksAlong(descriptor.axes[*x].length, chunk(geometry, *x)) *
                                 ChunksAlong(descriptor.axes[*y].length, chunk(geometry, *y));
    const std::uint64_t pixels = chunks * DecodedChunkBytes(descriptor, geometry);
    if (!apply_mask) {
        return pixels;
    }
    const auto& flag = flag_geometry.chunk_shape.empty() ? geometry : flag_geometry;
    std::uint64_t flag_chunks = 1;
    for (std::size_t axis = 0; axis < descriptor.axes.size(); ++axis) {
        const auto length = descriptor.axes[axis].length;
        const auto flag_chunk = chunk(flag, axis);
        if (axis == *x || axis == *y) {
            flag_chunks *= ChunksAlong(length, flag_chunk);
        } else {
            flag_chunks *= FlagChunksAcrossARun(length, chunk(geometry, axis), flag_chunk);
        }
    }
    return pixels + (flag_chunks * ChunkElements(flag));
}

Result<std::unique_ptr<ReadingAhead>> ReadingAhead::For(std::vector<std::shared_ptr<const RunSource>> sources) {
    if (sources.empty()) {
        return Error{ErrorCode::invalid_argument, "Reading ahead needs an image to read ahead of", {}};
    }
    // What each cache has to hold -- two runs of every image reading through it, the one playing and
    // the one decoded ahead -- and of which images, for saying so.
    struct Demand {
        std::uint64_t needed = 0;
        std::uint64_t holds = 0;
        std::string images;
    };
    std::map<const void*, Demand> caches;
    for (const auto& source : sources) {
        if (!source || source->PlaneRunBytes() == 0) {
            return Error{ErrorCode::invalid_argument,
                         "Reading ahead needs a plane of every image, and " +
                             (source ? source->Name() : std::string("one image")) + " has none",
                         source ? source->Name() : std::string{}};
        }
        const auto cache = source->Cache();
        auto& demand = caches[cache.identity];
        demand.needed += 2 * source->PlaneRunBytes();
        demand.holds = cache.bytes;
        demand.images += (demand.images.empty() ? "" : ", ") + source->Name();
    }
    for (const auto& [identity, demand] : caches) {
        if (demand.needed > demand.holds) {
            return Error{ErrorCode::buffer_too_small,
                         "Reading ahead needs room for two runs of chunks of each image in the cache it reads "
                         "through: " +
                             demand.images + " need " + Size(demand.needed) + " and their cache holds " +
                             Size(demand.holds),
                         {}};
        }
    }
    return std::unique_ptr<ReadingAhead>(new ReadingAhead(std::move(sources)));
}

ReadingAhead::ReadingAhead(std::vector<std::shared_ptr<const RunSource>> sources) : _sources(std::move(sources)) {}

ReadingAhead::~ReadingAhead() {
    Cancel();
    Join();
}

void ReadingAhead::Served(Clock::time_point began, bool late, const std::vector<AnimatedPlane>& shown,
                          const std::vector<std::vector<AnimatedPlane>>& upcoming) {
    std::vector<std::pair<std::shared_ptr<const RunSource>, ReadRequest>> work;
    {
        const std::scoped_lock lock(_mutex);
        // Whether the frame shared the machine with a prefetch: one had started by the time it began
        // and had not finished by then.
        const bool overlapped = !_prefetching.empty() && _started <= began && (_under_way || _finished > began);
        if (late && overlapped) {
            _stats.stopped = true;
        }
        if (overlapped && !_counted_catch) {
            for (const auto& plane : shown) {
                const auto prefetching = _prefetching.find(plane.image);
                if (prefetching != _prefetching.end() && plane.image < _sources.size() &&
                    _sources[plane.image]->RunOf(plane.request) == prefetching->second) {
                    ++_stats.caught_up;
                    _counted_catch = true;
                    break;
                }
            }
        }
        if (_stats.stopped || _cancelled || _under_way) {
            return;
        }

        std::map<std::size_t, Run> runs;
        const auto looked_at = std::min(upcoming.size(), kUpcomingFrames);
        for (const auto& now : shown) {
            if (now.image >= _sources.size() || runs.count(now.image) != 0) {
                continue;
            }
            const auto& source = _sources[now.image];
            const auto run_now = source->RunOf(now.request);
            for (std::size_t frame = 0; frame < looked_at; ++frame) {
                const auto& planes = upcoming[frame];
                const auto plane = std::find_if(planes.begin(), planes.end(), [&](const AnimatedPlane& candidate) {
                    return candidate.image == now.image;
                });
                if (plane == planes.end()) {
                    continue;
                }
                const auto run = source->RunOf(plane->request);
                if (run == run_now) {
                    continue;
                }
                const auto requested = _requested.find(now.image);
                if (requested == _requested.end() || requested->second != run) {
                    _requested[now.image] = run;
                    runs[now.image] = run;
                    work.emplace_back(source, plane->request);
                }
                break;
            }
        }
        if (work.empty()) {
            return;
        }
        _under_way = true;
        _started = Clock::now();
        _prefetching = std::move(runs);
        _counted_catch = false;
        _stats.prefetches += static_cast<unsigned>(work.size());
    }

    // The last prefetch has finished -- none is under way -- but its thread may not have returned.
    Join();
    try {
        _worker = std::thread([this, work = std::move(work)] {
            const auto cancelled = [this] { return _cancelled.load(); };
            bool failed = false;
            // Nothing above this thread can catch what escapes it -- an exception leaving a thread's
            // function is std::terminate, which takes the consumer with it -- and a prefetch can throw
            // before the image's own guard is reached: copying the options it reads with allocates.
            // Reading ahead is for memory to spare, so a prefetch that ran out stops it.
            try {
                for (const auto& [source, plane] : work) {
                    if (_cancelled) {
                        break;
                    }
                    (void)source->Prefetch(plane, cancelled);
                }
            } catch (...) {
                failed = true;
            }
            const std::scoped_lock lock(_mutex);
            _under_way = false;
            _finished = Clock::now();
            if (failed) {
                _stats.stopped = true;
            }
        });
    } catch (...) {
        // A machine with no thread to spare has no time to spare either, and one without the memory
        // to start a thread none to decode ahead into.
        const std::scoped_lock lock(_mutex);
        _under_way = false;
        _finished = Clock::now();
        _stats.stopped = true;
    }
}

void ReadingAhead::Cancel() {
    _cancelled = true;
    const std::scoped_lock lock(_mutex);
    _stats.stopped = true;
}

ReadAheadStats ReadingAhead::Stats() const {
    const std::scoped_lock lock(_mutex);
    ReadAheadStats stats = _stats;
    stats.under_way = _under_way;
    return stats;
}

bool ReadingAhead::UnderWay() const {
    const std::scoped_lock lock(_mutex);
    return _under_way;
}

void ReadingAhead::Join() {
    if (_worker.joinable()) {
        _worker.join();
    }
}

}  // namespace carta::zarr::internal
