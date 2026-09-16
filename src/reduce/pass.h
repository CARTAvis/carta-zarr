/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_PASS_H_
#define CARTA_ZARR_SRC_REDUCE_PASS_H_

#include "carta-zarr/carta_zarr.h"

#include "chunk_blocks.h"
#include "reduce/axis_map.h"
#include "store.h"
#include "zarr/pixel_reader.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace carta::zarr::internal {

/**
 * Where a pass gets its pixels.
 *
 * One adapter today, which reads them from a `Store`. It is a parameter rather than a `Store&`
 * reached for directly so that a second adapter -- one that answers from memory, and lets a test
 * state a cube the size of a real one without writing it to disk -- arrives without this signature
 * changing. Until that one exists this is a hypothetical seam, not a real one.
 *
 * Consulted once per slab, so the indirect call is paid per megabyte of pixels rather than per
 * pixel. That is the whole reason it can be an interface at all while the visitor cannot: see
 * ADR 0005.
 */
class SlabSource {
public:
    SlabSource() = default;
    SlabSource(const SlabSource&) = delete;
    SlabSource& operator=(const SlabSource&) = delete;
    SlabSource(SlabSource&&) = delete;
    SlabSource& operator=(SlabSource&&) = delete;
    virtual ~SlabSource() = default;

    virtual Result<void> ReadPixels(const zarr::PixelSelection& selection, float* destination,
                                    std::size_t elements, const ReadOptions& options) const = 0;
    // Only called when the plan says the image has a pixel mask to apply.
    virtual Result<void> ReadMask(const zarr::PixelSelection& selection, std::uint8_t* destination,
                                  std::size_t elements, const ReadOptions& options) const = 0;
};

// The adapter that serves production: the image's own data variable, and the flag that masks it.
class StoreSlabSource final : public SlabSource {
public:
    StoreSlabSource(const Store& store, const ImageDescriptor& descriptor)
        : _store(&store), _descriptor(&descriptor) {}

    Result<void> ReadPixels(const zarr::PixelSelection& selection, float* destination, std::size_t elements,
                            const ReadOptions& options) const override {
        return _store->ReadPixelsFloat32(_descriptor->id, selection, destination, elements, options);
    }

    Result<void> ReadMask(const zarr::PixelSelection& selection, std::uint8_t* destination,
                          std::size_t elements, const ReadOptions& options) const override {
        return _store->ReadPixelMaskBytes(_descriptor->pixel_mask_id, selection, destination, elements,
                                          options);
    }

private:
    const Store* _store;
    const ImageDescriptor* _descriptor;
};

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
    // Index into the channel range the pass was given, not an image channel.
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
struct PassPlan {
    const ImageDescriptor* descriptor = nullptr;
    AxisMap map;
    // u is the spatial axis the store varies fastest; v is the other one.
    std::size_t axis_u = 0;
    std::size_t axis_v = 0;
    std::uint64_t u_length = 0;
    std::uint64_t v_length = 0;
    std::uint64_t chunk_u = 1;
    std::uint64_t chunk_v = 1;
    std::uint64_t chunk_depth = 1;
    std::uint64_t least_channels = 1;
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
};

PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const AxisMap& map,
                  const Range& spectral, std::uint64_t polarization, std::uint64_t time,
                  std::uint64_t sample, const ReadOptions& options);

// The spectral range a pass was asked for has to fall inside the image. Compared by dividing the
// room that is left rather than by multiplying out the span: (count - 1) * stride wraps, and a
// wrapped span passes a check it should fail.
Result<void> ValidateSpectralRange(const ImageDescriptor& descriptor, const AxisMap& map,
                                   const Range& spectral);

/**
 * How many channels one emitted block may hold.
 *
 * The hint is the caller's, the budgets are the library's, and the chunk alignment is the pass's;
 * the smallest wins and the block reports what it used. Without a hint a block costs one budget of
 * decoded bytes -- the same invariant a piece of Read carries -- so it is free: the block spends
 * whatever the spatial walk left over. A small region leaves almost all of it and the block spans
 * many chunks along the spectrum; a region covering the image spends the budget spatially and the
 * block becomes the single chunk layer the pass is already reading.
 *
 * `layer_chunks` is the chunks one spectral layer of whatever the caller is walking occupies, which
 * is the plan's for a whole plane and the region set's own for a reduction.
 */
std::uint64_t PlanEmitChannels(const PassPlan& plan, std::uint64_t layer_chunks,
                               std::size_t bytes_per_channel, std::uint32_t hint);

// One slab to read, in the pass's own axes.
struct SlabRequest {
    std::uint64_t u_start = 0;
    std::uint64_t u_count = 0;
    std::uint64_t u_stride = 1;
    std::uint64_t v_start = 0;
    std::uint64_t v_count = 0;
    std::uint64_t v_stride = 1;
    // Index into the channel range the pass was given, and how many of them this slab holds.
    std::uint64_t channel_begin = 0;
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
Result<Slab> ReadSlab(const SlabSource& source, const PassPlan& plan, const ReadOptions& options,
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
 * Visit every plane of one channel range, a band of chunk rows at a time.
 *
 * `before_read` runs before each read after the first of the range, which is where a caller reports
 * what it has or decides to stop. `visit` receives one `Slab`.
 *
 * Both are template parameters and neither may become a `std::function`: the per-pixel loop inlines
 * through `visit`, and 54731c1 measured a quarter of a reduction riding on that. ADR 0005.
 */
template <typename BeforeRead, typename Visit>
Result<void> RunPass(const SlabSource& source, const PassPlan& plan, const ReadOptions& options,
                     std::uint64_t begin, std::uint64_t end, std::uint64_t& chunks_done,
                     BeforeRead&& before_read, Visit&& visit) {
    const auto& node = plan.descriptor->id;
    const Range spectral = plan.spectral;
    SlabBuffers buffers;
    std::uint64_t reads_done = 0;

    for (std::uint64_t v_begin = 0; v_begin < plan.v_length;) {
        const std::uint64_t v_end = std::min(plan.v_length, v_begin + (plan.band_rows * plan.chunk_v));
        std::uint64_t v_start = 0;
        std::uint64_t v_count = 0;
        SampledRange(v_begin, v_end, plan.sample, v_start, v_count);
        const std::uint64_t band_chunks = std::max<std::uint64_t>(
            1, (((plan.u_length - 1) / plan.chunk_u) + 1) * ((((v_end - v_begin) - 1) / plan.chunk_v) + 1));
        if (v_count == 0) {
            chunks_done += band_chunks * ((end - begin + plan.least_channels - 1) / plan.least_channels);
            v_begin = v_end;
            continue;
        }
        const std::uint64_t spectral_chunks =
            std::max<std::uint64_t>(1, plan.slab_budget_bytes / (band_chunks * plan.chunk_bytes));
        const std::uint64_t slab_channels = std::max<std::uint64_t>(1, spectral_chunks * plan.least_channels);

        for (std::uint64_t slab_begin = begin; slab_begin < end;) {
            const std::uint64_t slab_end = AlignedBlockEnd(slab_begin, std::min(slab_channels, end - slab_begin),
                                                           end, spectral.start, spectral.stride, plan.chunk_depth);
            const std::uint64_t slab_length = slab_end - slab_begin;

            if (reads_done > 0) {
                if (auto ready = before_read(chunks_done); !ready) {
                    return ready.error();
                }
            }
            ++reads_done;

            if (auto control = zarr::CheckReadControl(options, node); !control) {
                return control.error();
            }

            std::uint64_t u_start = 0;
            std::uint64_t u_count = 0;
            SampledRange(0, plan.u_length, plan.sample, u_start, u_count);

            SlabRequest slab_request;
            slab_request.u_start = u_start;
            slab_request.u_count = u_count;
            slab_request.u_stride = plan.sample;
            slab_request.v_start = v_start;
            slab_request.v_count = v_count;
            slab_request.v_stride = plan.sample;
            slab_request.channel_begin = slab_begin - begin;
            slab_request.channel_count = slab_length;

            auto slab = ReadSlab(source, plan, options, slab_request, buffers);
            if (!slab) {
                return slab.error();
            }

            visit(slab.value());

            chunks_done += band_chunks * ((slab_length + plan.least_channels - 1) / plan.least_channels);
            slab_begin = slab_end;
        }
        v_begin = v_end;
    }
    return {};
}

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PASS_H_
