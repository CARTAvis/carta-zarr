/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_PASS_H_
#define CARTA_ZARR_SRC_REDUCE_PASS_H_

#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include "chunk_blocks.h"
#include "pixel_source.h"
#include "reduce/axis_map.h"
#include "zarr/pixel_selection.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace carta::zarr::internal {

/**
 * One read of the pass, handed to the visitor.
 *
 * A pointer and three strides rather than a packed buffer, because the destination comes back in
 * the store's own order and packing it would be the transpose the pass exists to avoid.
 *
 * A read, not a plane: a visitor that splits the work across threads needs a piece big enough to
 * pay for the dispatch, and a plane of a few hundred thousand pixels is not one. A visitor that
 * wants planes loops over `channel_count` itself, which costs it nothing.
 */
struct Slab {
    // Index into the channel range this RunPass was given -- so a visitor accumulating into a block
    // of its own indexes by it directly. Not an index into the pass's spectral selection; see
    // SlabRequest::channel_index.
    std::uint64_t first_channel = 0;
    std::uint64_t channel_count = 0;
    const float* pixels = nullptr;
    std::uint64_t stride_u = 1;
    std::uint64_t stride_v = 1;
    std::uint64_t stride_z = 1;
    std::uint64_t u_count = 0;
    std::uint64_t v_count = 0;
};

/**
 * Everything a pass decides before a byte is read: which axis is contiguous, how wide a band is,
 * how deep a slab goes, and how much a read may decode.
 *
 * Pure, and separated from the pass for the same reason `PlanRowTasks` is separated from the pool:
 * it is the half with an answer worth checking, and checking it needs no store, no transport and no
 * fixture. The read strategy has moved more than once, and every time it moved it moved in three
 * places at once.
 */
class PassPlan;

PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const AxisMap& map,
                  const Range& spectral, std::uint64_t polarization, std::uint64_t time,
                  std::uint64_t sample, const ReadOptions& options);

class PassPlan {
public:
    const ImageDescriptor* descriptor = nullptr;
    AxisMap map;
    // u is the spatial axis the store varies fastest; v is the other one.
    std::size_t axis_u = 0;
    std::size_t axis_v = 0;
    std::uint64_t u_length = 0;
    std::uint64_t v_length = 0;
    std::uint64_t chunk_u = 1;
    std::uint64_t chunk_v = 1;
    // What one chunk costs to decode, counting the flag beside it when this read applies the mask,
    // and how much of that a single read may hold. Both are answers rather than steps towards one:
    // ADR 0005 turns on the first, and the second is the caller's own ceiling when it stated one.
    std::uint64_t chunk_bytes = 1;
    std::size_t slab_budget_bytes = 0;
    std::uint64_t band_rows = 1;
    std::uint64_t layer_chunks = 1;
    Range spectral;
    std::uint64_t polarization = 0;
    std::uint64_t time = 0;
    // Take every nth pixel along both spatial axes. This does not reduce the chunks a read decodes
    // -- a chunk comes back whole however few of its pixels are wanted -- so it pays only when it
    // steps over chunks entirely.
    std::uint64_t sample = 1;
    bool apply_mask = false;

    // How many chunks along the spectrum a run of `channels` selected channels covers.
    //
    // A walk counts its progress in chunks, and the spectral axis is the one where a selection's
    // stride makes that not simply a division. Six call sites wrote this out; they disagreed about
    // nothing, which is the argument for saying it once rather than the argument for leaving it.
    std::uint64_t ChunksFor(std::uint64_t channels) const {
        return (channels + _least_channels - 1) / _least_channels;
    }


    // How many units of `chunks_per_unit` chunks one read's budget affords. Never zero: a budget
    // smaller than a single unit still reads one, because a chunk is the smallest thing that can be
    // decoded and refusing to read is the worse answer.
    std::uint64_t UnitsAffordable(std::uint64_t chunks_per_unit) const {
        return std::max<std::uint64_t>(
            1, slab_budget_bytes / std::max<std::uint64_t>(1, chunks_per_unit * chunk_bytes));
    }

    // How many channels one slab may hold when its spatial footprint occupies `footprint_chunks`
    // chunks of each spectral chunk it touches.
    std::uint64_t SlabChannels(std::uint64_t footprint_chunks) const {
        return std::max<std::uint64_t>(1, UnitsAffordable(footprint_chunks) * _least_channels);
    }

    // The end of a slab that begins at `begin` and would like to be `desired` channels long, moved
    // onto a chunk boundary so that no decode serves two slabs.
    std::uint64_t AlignedSlabEnd(std::uint64_t begin, std::uint64_t desired, std::uint64_t end) const {
        return AlignedBlockEnd(begin, desired, end, spectral.start, spectral.stride, _chunk_depth);
    }

    // How many channels one emitted block may hold.
    //
    // The hint is the caller's, the budgets are the library's, and the chunk alignment is the
    // plan's; the smallest wins and the block reports what it used. Without a hint a block costs one
    // budget of decoded bytes -- the same invariant a piece of Read carries -- so it is free: the
    // block spends whatever the spatial walk left over. A small region leaves almost all of it and
    // the block spans many chunks along the spectrum; a region covering the image spends the budget
    // spatially and the block becomes the single chunk layer the walk is already reading.
    //
    // `layer_chunks` is the chunks one spectral layer of whatever the caller is walking occupies,
    // which is the plan's own for a whole plane and the region set's for a reduction.
    std::uint64_t EmitChannels(std::uint64_t layer_chunks, std::size_t bytes_per_channel,
                               std::uint32_t hint) const;

private:
    friend PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const AxisMap& map,
                             const Range& spectral, std::uint64_t polarization, std::uint64_t time,
                             std::uint64_t sample, const ReadOptions& options);

    // Steps towards the answers above rather than answers themselves, and the two a caller used to
    // divide by itself: the chunk-count rule was written out in six places and the slab-sizing rule
    // in two. Nothing asserts either directly -- what a test has to say about _least_channels it
    // says through ChunksFor, which is the question a caller actually asks.
    std::uint64_t _chunk_depth = 1;
    std::uint64_t _least_channels = 1;
};

// The spectral range a pass was asked for has to fall inside the image. Compared by dividing the
// room that is left rather than by multiplying out the span: (count - 1) * stride wraps, and a
// wrapped span passes a check it should fail.
Result<void> ValidateSpectralRange(const ImageDescriptor& descriptor, const AxisMap& map,
                                   const Range& spectral);

// One slab to read, in the pass's own axes.
struct SlabRequest {
    std::uint64_t u_start = 0;
    std::uint64_t u_count = 0;
    std::uint64_t u_stride = 1;
    std::uint64_t v_start = 0;
    std::uint64_t v_count = 0;
    std::uint64_t v_stride = 1;
    // Which channels to read, as an index into the pass's own spectral selection -- so channel
    // `channel_index` of the image is `spectral.start + channel_index * spectral.stride`. Absolute,
    // unlike Slab::first_channel, which is relative to the range one RunPass was given. They are
    // the same number only when a pass starts at the beginning of the selection, which is why
    // confusing them is invisible until something asks for a later block.
    std::uint64_t channel_index = 0;
    std::uint64_t channel_count = 0;
};

// Reused across slabs, so that a pass allocates once rather than once per read.
struct SlabBuffers {
    std::vector<float> pixels;
    std::vector<std::uint8_t> mask;
};

/**
 * Read one slab and hand back how to walk it.
 *
 * This is what every pass over a cube has in common, whatever order it visits chunks in: ask for
 * the stored dimensions reversed so the plane arrives untransposed, derive the strides of what
 * comes back, read the pixels, and -- when the image has a flag the caller did not decline -- read
 * that too and fold it into the pixels, because a flagged pixel and a NaN pixel mean the same thing
 * to everything downstream.
 *
 * The returned Slab points into `buffers`, so it is valid until the next call with them.
 */
Result<Slab> ReadSlab(const PixelSource& source, const PassPlan& plan, const ReadOptions& options,
                      const SlabRequest& request, SlabBuffers& buffers);

// Samples of `stride` that fall in [begin, end), as a start and a count.
inline void SampledRange(std::uint64_t begin, std::uint64_t end, std::uint64_t stride, std::uint64_t& start,
                         std::uint64_t& count) {
    const std::uint64_t first = (begin + stride - 1) / stride;
    const std::uint64_t last = (end + stride - 1) / stride;
    start = first * stride;
    count = last > first ? last - first : 0;
}

/**
 * One spatial footprint a run of slabs is read over, in the pass's own axes.
 *
 * `chunks` is what that footprint occupies in one chunk of the spectral axis, and it is doing two
 * jobs: it sizes the slab, because the budget is over chunk data rather than over the pixels kept,
 * and it is the unit progress is counted in. A caller that gets it wrong reads the right pixels and
 * reports a bar that lies.
 */
struct SlabFootprint {
    std::uint64_t u_start = 0;
    std::uint64_t u_count = 0;
    std::uint64_t u_stride = 1;
    std::uint64_t v_start = 0;
    std::uint64_t v_count = 0;
    std::uint64_t v_stride = 1;
    std::uint64_t chunks = 1;
};

/**
 * Walks one spatial footprint along the spectrum, a slab at a time.
 *
 * Which footprints to visit is the caller's: a whole-plane pass bands the plane, and a reduction
 * cuts the chunk runs its regions occupy. Those two are genuinely different walks and stay that
 * way. What they had in common was everything inside one footprint -- how deep a slab goes, where
 * it is cut so that no decode serves two slabs, when the caller is told what it has, where
 * cancellation is checked, and how progress is counted -- and that was written out twice, in full,
 * with the subtlety intact in both copies: `reads_done` guards the report so that a footprint taken
 * in a single read still reports once at the end rather than before it starts.
 *
 * `reads_done` and `chunks_done` are the caller's because their spans are: a reduction resets them
 * per emitted block, a whole-plane pass keeps them for the call.
 *
 * `visit` is a template parameter and must not become a `std::function`: the per-pixel loop inlines
 * through it, and 54731c1 measured a quarter of a reduction riding on that. ADR 0005. `report` is
 * one call per slab, which is the same footing `PixelSource` stands on.
 *
 * Holds the buffers, so a walk allocates once rather than once per footprint.
 */
class SlabWalk {
public:
    SlabWalk(const PixelSource& source, const PassPlan& plan, const ReadOptions& options)
        : _source(source), _plan(plan), _options(options) {}

    SlabWalk(const SlabWalk&) = delete;
    SlabWalk& operator=(const SlabWalk&) = delete;

    template <typename Report, typename Visit>
    Result<void> Over(const SlabFootprint& footprint, std::uint64_t begin, std::uint64_t end,
                      std::uint64_t& reads_done, std::uint64_t& chunks_done, Report&& report, Visit&& visit) {
        const std::uint64_t slab_channels = _plan.SlabChannels(footprint.chunks);

        for (std::uint64_t slab_begin = begin; slab_begin < end;) {
            const std::uint64_t slab_end =
                _plan.AlignedSlabEnd(slab_begin, std::min(slab_channels, end - slab_begin), end);
            const std::uint64_t slab_length = slab_end - slab_begin;

            // What is in hand before spending another budget. A footprint that takes one read never
            // gets here, so a small one still reports once -- at the end, through its caller.
            if (reads_done > 0) {
                if (auto ready = report(chunks_done); !ready) {
                    return ready.error();
                }
            }
            ++reads_done;

            if (auto control = zarr::CheckReadControl(_options, _plan.descriptor->id); !control) {
                return control.error();
            }

            SlabRequest slab_request;
            slab_request.u_start = footprint.u_start;
            slab_request.u_count = footprint.u_count;
            slab_request.u_stride = footprint.u_stride;
            slab_request.v_start = footprint.v_start;
            slab_request.v_count = footprint.v_count;
            slab_request.v_stride = footprint.v_stride;
            slab_request.channel_index = slab_begin;
            slab_request.channel_count = slab_length;

            auto slab = ReadSlab(_source, _plan, _options, slab_request, _buffers);
            if (!slab) {
                return slab.error();
            }
            // What was read is absolute; what the visitor indexes by is relative to this walk.
            slab.value().first_channel = slab_begin - begin;

            visit(slab.value());

            chunks_done += footprint.chunks * _plan.ChunksFor(slab_length);
            slab_begin = slab_end;
        }
        return {};
    }

private:
    const PixelSource& _source;
    const PassPlan& _plan;
    const ReadOptions& _options;
    SlabBuffers _buffers;
};

/**
 * Visit every plane of one channel range, a band of chunk rows at a time.
 *
 * `before_read` runs before each read after the first of the range, which is where a caller reports
 * what it has or decides to stop. `visit` receives one `Slab`.
 *
 * Both are template parameters and neither may become a `std::function`: the per-pixel loop inlines
 * through `visit`, and 54731c1 measured a quarter of a reduction riding on that. ADR 0005.
 */
template <typename BeforeRead, typename Visit>
Result<void> RunPass(const PixelSource& source, const PassPlan& plan, const ReadOptions& options,
                     std::uint64_t begin, std::uint64_t end, std::uint64_t& chunks_done,
                     BeforeRead&& before_read, Visit&& visit) {
    SlabWalk walk(source, plan, options);
    std::uint64_t reads_done = 0;

    for (std::uint64_t v_begin = 0; v_begin < plan.v_length;) {
        const std::uint64_t v_end = std::min(plan.v_length, v_begin + (plan.band_rows * plan.chunk_v));
        SlabFootprint band;
        SampledRange(v_begin, v_end, plan.sample, band.v_start, band.v_count);
        band.chunks = std::max<std::uint64_t>(
            1, (((plan.u_length - 1) / plan.chunk_u) + 1) * ((((v_end - v_begin) - 1) / plan.chunk_v) + 1));
        if (band.v_count == 0) {
            chunks_done += band.chunks * plan.ChunksFor(end - begin);
            v_begin = v_end;
            continue;
        }
        SampledRange(0, plan.u_length, plan.sample, band.u_start, band.u_count);
        band.u_stride = plan.sample;
        band.v_stride = plan.sample;

        if (auto walked = walk.Over(band, begin, end, reads_done, chunks_done, before_read, visit); !walked) {
            return walked.error();
        }
        v_begin = v_end;
    }
    return {};
}
}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PASS_H_
