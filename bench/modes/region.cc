/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// region: ReduceSpectral over a box covering a share of the plane, every channel of it.

#include "cube.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace carta::zarr::bench {

namespace {

// How many boxes of `side` (a share of the axis) fit along an axis of `length` pixels: at 5% of the
// plane a side is 0.2236 of the axis, so four, each in a cell of 0.25 with room to move.
std::uint64_t GridCells(std::uint64_t length, double side) {
    const auto cells = static_cast<std::uint64_t>(std::floor(1.0 / side));
    return std::clamp<std::uint64_t>(cells, 1, std::max<std::uint64_t>(length, 1));
}

class RegionRunner final : public CubeRunner {
public:
    RegionRunner(Image image, CubeAxes axes) : CubeRunner(std::move(image), axes) {
        // Sized here so that no operation pays for growing it.
        _exact.reserve(axes.channels * 4);
    }

    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) override {
        _exact.clear();
        const RegionMask region{operation.x, operation.y, operation.width, operation.height, {}};
        SpectralReduceRequest request;
        request.planes.spectral = Range{0, axes().channels, 1};
        request.planes.polarization = operation.polarization;
        request.regions = {&region, 1};
        request.statistics = Statistic::num_pixels | Statistic::nan_count | Statistic::sum | Statistic::sum_sq |
                             Statistic::min | Statistic::max;

        _exact.assign(axes().channels * 4, std::numeric_limits<double>::quiet_NaN());
        const auto sink = [this](const SpectralBlock& block) {
            if (!block.complete) {
                return true;
            }
            for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
                const auto totals = block.Totals(0, channel);
                auto* exact = &_exact[(block.first_channel + channel) * 4];
                exact[0] = totals.num_pixels;
                exact[1] = totals.nan_count;
                exact[2] = totals.min;
                exact[3] = totals.max;
            }
            return true;
        };
        if (auto reduced = image().ReduceSpectral(request, sink, options); !reduced) {
            _exact.clear();
            return std::move(reduced).error();
        }
        return operation.width * operation.height * axes().channels;
    }

    // Each channel's pixel count, NaN count and extremes, which are exact; not its sums.
    std::uint64_t Fingerprint() const override {
        return bench::Fingerprint(_exact.data(), _exact.size());
    }

private:
    std::vector<double> _exact;
};

class Region final : public Workload {
public:
    Region(const RunOptions& options, RegionSettings settings) : Workload(Mode::region, options), _settings(settings) {}

    void WriteSettings(Row& row) const override {
        std::array<char, 32> fraction{};
        std::snprintf(fraction.data(), fraction.size(), "%.4f", _settings.fraction);
        row.region_fraction = fraction.data();
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
        return PlanDraws(mode(), at, ops, cells_x * cells_y,
                         [&](std::uint64_t value, Stream& details, Operation& operation) {
                             operation.width = box_width;
                             operation.height = box_height;
                             operation.x = (value % cells_x) * cell_width + details.Below(cell_width - box_width + 1);
                             operation.y =
                                 (value / cells_x) * cell_height + details.Below(cell_height - box_height + 1);
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
        return ChunkBox::Spanning({operation.x, operation.y, 0, operation.polarization},
                                  {operation.x + operation.width, operation.y + operation.height, axes.channels,
                                   operation.polarization + 1},
                                  axes, chunk_shape);
    }

    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        (void)context;
        auto axes = CubeAxes::Of(image.descriptor());
        if (!axes) {
            return std::move(axes).error();
        }
        return std::unique_ptr<Runner>(std::make_unique<RegionRunner>(std::move(image), *axes));
    }

private:
    RegionSettings _settings;
};

}  // namespace

std::unique_ptr<Workload> RegionWorkload(const RunOptions& options, RegionSettings region) {
    return std::make_unique<Region>(options, region);
}

}  // namespace carta::zarr::bench
