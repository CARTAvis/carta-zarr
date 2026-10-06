/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// cube-histogram: a histogram of this process's share of the channels, as carta-backend computes one
// with the --zarr_histogram_method it is given.

#include "cube.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace carta::zarr::bench {

namespace {

// CARTA's automatic bin count for a plane: the square root of its pixel count, and at least two.
std::uint32_t AutomaticBins(const CubeAxes& axes) {
    const auto bins = std::sqrt(static_cast<double>(axes.width) * static_cast<double>(axes.height));
    return static_cast<std::uint32_t>(std::max(bins, 2.0));
}

// Read through a cache pool that keeps nothing, as carta-backend's cube walks do, so that a scan of
// the whole cube neither fills nor finds the context's cache.
class HistogramRunner final : public CubeRunner {
public:
    HistogramRunner(const Context& context, Image image, CubeAxes axes, HistogramMethod method)
        : CubeRunner(std::move(image), axes), _method(method) {
        if (auto pool = context.NewCachePool(0)) {
            _keeping_nothing = *pool;
        }
        // Sized here so that no operation pays for growing it.
        _exact.reserve(axes.channels * 4);
    }

    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& passed) override {
        _exact.clear();
        _counts.clear();
        auto options = passed;
        options.control.cache_pool = _keeping_nothing;
        if (_method.kind == HistogramMethod::Kind::exact) {
            return Exact(operation, options);
        }
        CubeHistogramRequest request;
        request.planes.spectral = Range{operation.channel, operation.channel_count, 1};
        request.planes.polarization = operation.polarization;
        request.bins = AutomaticBins(axes());
        request.spatial_sample = _method.kind == HistogramMethod::Kind::sampled ? _method.stride : 1;
        auto histogram = image().ComputeCubeHistogram(request, options);
        if (!histogram) {
            return std::move(histogram).error();
        }
        const auto& totals = histogram->totals;
        _exact = {totals.num_pixels, totals.nan_count, totals.min, totals.max};
        return axes().width * axes().height * operation.channel_count;
    }

    // The pixel count and extremes, which are exact, and an exact histogram's counts, since binning over
    // fixed bounds is exact; not a one-pass histogram's, which ComputeCubeHistogram says depend on the
    // thread count.
    std::uint64_t Fingerprint() const override {
        auto hash = bench::Fingerprint(_exact.data(), _exact.size());
        for (const auto count : _counts) {
            Mix(hash, count, sizeof(count));
        }
        return hash;
    }

private:
    // As carta-backend computes a cube histogram by default: a reduction over whole planes for the
    // cube's range, then every plane binned over it and the planes added together.
    Result<std::uint64_t> Exact(const Operation& operation, const ReadOptions& options) {
        const RegionMask plane{0, 0, axes().width, axes().height, {}};
        SpectralReduceRequest range;
        range.planes.spectral = Range{operation.channel, operation.channel_count, 1};
        range.planes.polarization = operation.polarization;
        range.regions = {&plane, 1};
        // The statistics carta-backend asks of that reduction, sum_sq_dev included; see RegionRunner.
        range.statistics = Statistic::num_pixels | Statistic::sum | Statistic::sum_sq | Statistic::min |
                           Statistic::max | Statistic::sum_sq_dev;

        double pixels = 0.0;
        double lowest = std::numeric_limits<double>::quiet_NaN();
        double highest = std::numeric_limits<double>::quiet_NaN();
        const auto extremes = [&](const SpectralBlock& block) {
            if (!block.complete) {
                return true;
            }
            for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
                const auto totals = block.Totals(0, channel);
                pixels += totals.num_pixels;
                // fmin and fmax pass over NaN, which is what a plane with nothing finite reports.
                lowest = std::fmin(lowest, totals.min);
                highest = std::fmax(highest, totals.max);
            }
            return true;
        };
        if (auto reduced = image().ReduceSpectral(range, extremes, options); !reduced) {
            return std::move(reduced).error();
        }
        _exact = {pixels, lowest, highest};

        // The backend declines an empty or inverted range, and so does this: there is nothing to bin.
        if (lowest < highest) {
            HistogramRequest bins;
            bins.planes = range.planes;
            bins.bins = AutomaticBins(axes());
            bins.lower = lowest;
            bins.upper = highest;
            _counts.assign(bins.bins, 0);
            const auto add = [this](const HistogramBlock& block) {
                if (!block.complete) {
                    return true;
                }
                for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
                    const auto* counts = block.Counts(channel);
                    for (std::size_t bin = 0; bin < block.bin_count; ++bin) {
                        _counts[bin] += counts[bin];
                    }
                }
                return true;
            };
            if (auto binned = image().ComputeHistogram(bins, add, options); !binned) {
                _counts.clear();
                return std::move(binned).error();
            }
        }
        return axes().width * axes().height * operation.channel_count;
    }

    HistogramMethod _method;
    std::optional<CachePool> _keeping_nothing;
    std::vector<double> _exact;
    std::vector<std::uint64_t> _counts;
};

class CubeHistogram final : public Workload {
public:
    CubeHistogram(const RunOptions& options, HistogramMethod method)
        : Workload(Mode::cube_histogram, options), _method(method) {}

    void WriteSettings(Row& row) const override { row.histogram_method = _method.Spell(); }

    // The operation covers the plane, so process p takes the p-th of the trial's contiguous runs of
    // channels, and each of its operations a polarization of its own.
    std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const override {
        std::vector<Operation> operations(ops);
        // With more processes than channels, each takes one channel and some of them share it.
        const bool shared = at.processes > axes.channels;
        std::uint64_t start = at.process_index % axes.channels;
        std::uint64_t end = start + 1;
        if (!shared) {
            start = axes.channels * at.process_index / at.processes;
            end = axes.channels * (at.process_index + 1) / at.processes;
        }
        Stream details(at.seed, mode(), at.trial, 2);
        const auto base = details.Below(axes.polarizations);
        for (unsigned index = 0; index < ops; ++index) {
            auto& operation = operations[index];
            operation.channel = start;
            operation.channel_count = end - start;
            operation.polarization = (base + index) % axes.polarizations;
            operation.overlap = shared || index >= axes.polarizations;
        }
        return operations;
    }

    std::string Describe(const Operation& operation) const override { return DescribeChannels(operation); }

    std::optional<ChunkBox> ChunksRead(const Operation& operation, const CubeAxes& axes,
                                       const std::vector<std::uint64_t>& chunk_shape) const override {
        return ChunkBox::Spanning(
            {0, 0, operation.channel, operation.polarization},
            {axes.width, axes.height, operation.channel + operation.channel_count, operation.polarization + 1}, axes,
            chunk_shape);
    }

    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        auto axes = CubeAxes::Of(image.descriptor());
        if (!axes) {
            return std::move(axes).error();
        }
        return std::unique_ptr<Runner>(std::make_unique<HistogramRunner>(context, std::move(image), *axes, _method));
    }

private:
    HistogramMethod _method;
};

}  // namespace

std::unique_ptr<Workload> CubeHistogramWorkload(const RunOptions& options, HistogramMethod histogram) {
    return std::make_unique<CubeHistogram>(options, histogram);
}

}  // namespace carta::zarr::bench
