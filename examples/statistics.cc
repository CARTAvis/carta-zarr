/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "support.h"

#include <array>
#include <numeric>

int main(int argc, char** argv) {
    using namespace carta::zarr;
    if (argc != 2) {
        std::cerr << "Usage: carta_zarr_statistics PATH_TO_DATASET.zarr\n";
        return 1;
    }
    const auto context = Context::Create();
    if (!context) {
        return examples::Report(context.error());
    }
    const auto opened = examples::OpenFirstImage(context.value(), argv[1]);
    if (!opened) {
        return examples::Report(opened.error());
    }
    const auto& image = opened.value();
    const auto patch = examples::Patch(image);
    if (!patch) {
        return examples::Report(patch.error());
    }
    const auto& axes = image.descriptor().axes;
    const auto x = *AxisIndex(axes, AxisRole::spatial_x);
    const auto y = *AxisIndex(axes, AxisRole::spatial_y);
    const auto spectral = *AxisIndex(axes, AxisRole::spectral);
    const PlaneSelection planes{{0, std::min<std::uint64_t>(4, axes[spectral].length), 1}, 0, 0};

    // Reduce a patch and its first pixel together. The region views stay alive
    // until the operation returns; no full cube is allocated by the caller.
    const std::array<RegionMask, 2> regions{
        {{0, 0, patch.value().axes[x].count, patch.value().axes[y].count}, {0, 0, 1, 1}}};
    SpectralReduceRequest reduction;
    reduction.planes = planes;
    reduction.regions = {regions.data(), regions.size()};
    reduction.statistics = Statistic::num_pixels | Statistic::sum | Statistic::sum_sq_dev;
    std::uint64_t reduced_channels = 0;
    const auto reduced = image.ReduceSpectral(reduction, [&](const SpectralBlock& block) {
        if (!block.complete) {
            return true;  // Final blocks replace the estimates emitted earlier.
        }
        reduced_channels += block.channel_count;
        for (std::size_t region = 0; region < block.region_count; ++region) {
            for (std::uint64_t i = 0; i < block.channel_count; ++i) {
                const auto totals = block.Totals(region, i);
                const auto channel = planes.spectral.start + (block.first_channel + i) * planes.spectral.stride;
                std::cout << "region=" << region << " channel=" << channel << " finite=" << totals.num_pixels;
                if (totals.num_pixels > 0) {
                    std::cout << " mean=" << totals.sum / totals.num_pixels;
                }
                if (totals.num_pixels > 1) {
                    std::cout << " variance=" << totals.sum_sq_dev / (totals.num_pixels - 1);
                }
                std::cout << '\n';
            }
        }
        return true;
    });
    if (!reduced) {
        return examples::Report(reduced.error());
    }

    // Histograms cover whole planes, independently of the regions above. This
    // range is explicitly chosen by this application; outliers are excluded.
    HistogramRequest histogram;
    histogram.planes = planes;
    histogram.bins = 64;
    histogram.lower = -1e6;
    histogram.upper = 1e6;
    std::uint64_t histogram_channels = 0;
    const auto binned = image.ComputeHistogram(histogram, [&](const HistogramBlock& block) {
        if (!block.complete) {
            return true;
        }
        histogram_channels += block.channel_count;
        for (std::uint64_t i = 0; i < block.channel_count; ++i) {
            const auto* counts = block.Counts(i);
            const auto total = std::accumulate(counts, counts + block.bin_count, std::uint64_t{0});
            const auto channel = planes.spectral.start + (block.first_channel + i) * planes.spectral.stride;
            std::cout << "histogram channel=" << channel << " counted=" << total << '\n';
        }
        return true;
    });
    if (!binned) {
        return examples::Report(binned.error());
    }

    CubeHistogramRequest cube;
    cube.planes = planes;
    cube.bins = 64;
    const auto combined = image.ComputeCubeHistogram(cube);
    if (!combined) {
        return examples::Report(combined.error());
    }
    const auto& result = combined.value();
    std::cout << "cube finite=" << result.totals.num_pixels << " bins=" << result.counts.size();
    if (result.totals.num_pixels > 0) {
        std::cout << " range=[" << result.totals.min << ", " << result.totals.max << ']';
    }
    std::cout << '\n';
    // Check that the streaming consumers received all requested final channels.
    const bool complete = reduced_channels == planes.spectral.count && histogram_channels == planes.spectral.count &&
                          result.counts.size() == cube.bins;
    return complete ? 0 : 1;
}
