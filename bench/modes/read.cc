/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// plane and spectrum: one whole l x m plane at a random channel, and every channel at one random
// pixel. Each times a first touch, reading through a cache pool of its own made before the clock
// starts, so that what an earlier operation decoded never answers a later one.

#include "cube.h"

#include <algorithm>
#include <optional>
#include <vector>

namespace carta::zarr::bench {

Result<std::size_t> ReadPixels(const Image& image, const CubeAxes& axes, const Operation& operation, bool plane,
                               float* pixels, const ReadOptions& options) {
    ReadRequest request;
    request.axes.assign(axes.rank, Range{0, 1, 1});
    request.axes[axes.x] = plane ? Range{0, axes.width, 1} : Range{operation.x, 1, 1};
    request.axes[axes.y] = plane ? Range{0, axes.height, 1} : Range{operation.y, 1, 1};
    request.axes[axes.spectral] = plane ? Range{operation.channel, 1, 1} : Range{0, axes.channels, 1};
    if (axes.polarization) {
        request.axes[*axes.polarization] = Range{operation.polarization, 1, 1};
    }
    const std::size_t elements = plane ? axes.width * axes.height : axes.channels;
    return image.Read(request, {pixels, elements}, options);
}

std::string DescribeChannels(const Operation& operation) {
    return "pol=" + std::to_string(operation.polarization) + ";chan=" + std::to_string(operation.channel) + ":" +
           std::to_string(operation.channel + operation.channel_count);
}

namespace {

// A plane or a spectrum through a cache pool of `first_touch_bytes` made for each operation, with the
// previous one let go of first, so that two pools' worth of chunks are never held at once. Not a pool
// that keeps nothing: a sharded layout keeps its shard index in the pool, and one that kept nothing
// would read the index again for every chunk of a single read.
class FirstTouchRunner final : public CubeRunner {
public:
    FirstTouchRunner(const Context& context, Image image, CubeAxes axes, bool plane, std::size_t first_touch_bytes)
        : CubeRunner(std::move(image), axes), _context(context), _plane(plane), _first_touch_bytes(first_touch_bytes) {
        // Sized here so that no operation pays for growing it.
        _pixels.resize(std::max(axes.width * axes.height, axes.channels));
    }

    Result<void> Prepare(const Operation& operation) override {
        (void)operation;
        _first_touch.reset();
        auto pool = _context.NewCachePool(_first_touch_bytes);
        if (!pool) {
            return std::move(pool).error();
        }
        _first_touch = *pool;
        return {};
    }

    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) override {
        _pixel_count = 0;
        auto own = options;
        if (_first_touch) {
            own.control.cache_pool = _first_touch;
        }
        auto written = ReadPixels(image(), axes(), operation, _plane, _pixels.data(), own);
        if (!written) {
            return std::move(written).error();
        }
        _pixel_count = *written;
        return static_cast<std::uint64_t>(*written);
    }

    std::uint64_t Fingerprint() const override {
        return bench::Fingerprint(_pixels.data(), _pixel_count);
    }

private:
    Context _context;
    bool _plane;
    std::size_t _first_touch_bytes;
    std::optional<CachePool> _first_touch;
    std::vector<float> _pixels;
    std::size_t _pixel_count = 0;
};

class FirstTouchWorkload final : public Workload {
public:
    FirstTouchWorkload(Mode mode, const RunOptions& options) : Workload(mode, options) {}

    std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const override {
        if (plane()) {
            return PlanDraws(mode(), at, ops, axes.channels * axes.polarizations,
                             [&](std::uint64_t value, Stream&, Operation& operation) {
                                 operation.channel = value % axes.channels;
                                 operation.polarization = value / axes.channels;
                             });
        }
        return PlanDraws(mode(), at, ops, axes.width * axes.height,
                         [&](std::uint64_t value, Stream& details, Operation& operation) {
                             operation.x = value % axes.width;
                             operation.y = value / axes.width;
                             operation.polarization = details.Below(axes.polarizations);
                         });
    }

    std::string Describe(const Operation& operation) const override {
        const auto pol = "pol=" + std::to_string(operation.polarization);
        if (plane()) {
            return pol + ";chan=" + std::to_string(operation.channel);
        }
        return pol + ";l=" + std::to_string(operation.x) + ";m=" + std::to_string(operation.y);
    }

    std::optional<ChunkBox> ChunksRead(const Operation& operation, const CubeAxes& axes,
                                       const std::vector<std::uint64_t>& chunk_shape) const override {
        if (plane()) {
            return ChunkBox::Spanning({0, 0, operation.channel, operation.polarization},
                                      {axes.width, axes.height, operation.channel + 1, operation.polarization + 1},
                                      axes, chunk_shape);
        }
        return ChunkBox::Spanning({operation.x, operation.y, 0, operation.polarization},
                                  {operation.x + 1, operation.y + 1, axes.channels, operation.polarization + 1}, axes,
                                  chunk_shape);
    }

    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        auto axes = CubeAxes::Of(image.descriptor());
        if (!axes) {
            return std::move(axes).error();
        }
        return std::unique_ptr<Runner>(std::make_unique<FirstTouchRunner>(context, std::move(image), *axes, plane(),
                                                                          options().FirstTouchCacheBytes()));
    }

private:
    bool plane() const noexcept {
        return mode() == Mode::plane;
    }
};

}  // namespace

std::unique_ptr<Workload> PlaneWorkload(const RunOptions& options) {
    return std::make_unique<FirstTouchWorkload>(Mode::plane, options);
}

std::unique_ptr<Workload> SpectrumWorkload(const RunOptions& options) {
    return std::make_unique<FirstTouchWorkload>(Mode::spectrum, options);
}

}  // namespace carta::zarr::bench
