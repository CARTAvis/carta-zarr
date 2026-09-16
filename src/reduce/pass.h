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
    const auto& descriptor = *plan.descriptor;
    const auto& node = descriptor.id;
    const auto rank = descriptor.axes.size();
    const Range spectral = plan.spectral;
    std::vector<float> pixels;
    std::vector<std::uint8_t> mask;
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

            ReadRequest read_request;
            read_request.axes.assign(rank, Range{0, 1, 1});
            read_request.axes.at(plan.axis_u) = Range{u_start, u_count, plan.sample};
            read_request.axes.at(plan.axis_v) = Range{v_start, v_count, plan.sample};
            read_request.axes.at(plan.map.spectral) =
                Range{spectral.start + (slab_begin * spectral.stride), slab_length, spectral.stride};
            if (plan.map.has_polarization) {
                read_request.axes.at(plan.map.polarization) = Range{plan.polarization, 1, 1};
            }
            if (plan.map.has_time) {
                read_request.axes.at(plan.map.time) = Range{plan.time, 1, 1};
            }

            auto selection = zarr::BuildSelection(descriptor, read_request);
            if (!selection) {
                return selection.error();
            }
            // Ask for the stored dimensions reversed, which against a destination whose dimension 0
            // is fastest is asking for no transpose at all.
            for (std::size_t i = 0; i < rank; ++i) {
                selection.value().logical_to_stored.at(i) = rank - 1 - i;
            }

            std::vector<std::uint64_t> stored_stride(rank, 1);
            std::uint64_t running = 1;
            for (std::size_t stored = rank; stored-- > 0;) {
                stored_stride.at(stored) = running;
                running *= selection.value().count.at(stored);
            }
            Slab slab;
            slab.first_channel = slab_begin - begin;
            slab.channel_count = slab_length;
            slab.stride_u = stored_stride.at(descriptor.axes.at(plan.axis_u).storage_index);
            slab.stride_v = stored_stride.at(descriptor.axes.at(plan.axis_v).storage_index);
            slab.stride_z = stored_stride.at(descriptor.axes.at(plan.map.spectral).storage_index);
            slab.u_count = u_count;
            slab.v_count = v_count;
            const auto elements = static_cast<std::size_t>(running);

            pixels.resize(elements);
            if (auto read = source.ReadPixels(selection.value(), pixels.data(), pixels.size(), options); !read) {
                return read.error();
            }
            if (plan.apply_mask) {
                mask.resize(elements);
                if (auto read = source.ReadMask(selection.value(), mask.data(), mask.size(), options); !read) {
                    return read.error();
                }
                for (std::size_t i = 0; i < elements; ++i) {
                    if (mask[i] == 0) {
                        pixels[i] = std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }
            slab.pixels = pixels.data();

            visit(slab);

            chunks_done += band_chunks * ((slab_length + plan.least_channels - 1) / plan.least_channels);
            slab_begin = slab_end;
        }
        v_begin = v_end;
    }
    return {};
}

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PASS_H_
