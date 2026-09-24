/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What stands in for a regression test on the pass.
//
// Nothing in the suite asserts a duration, and the one change that made the accumulation loop
// vectorise was worth about a quarter of a whole-region reduction while measuring as exactly
// nothing in a Debug build -- see ADR 0005. So the guard against collapsing the three walks into
// one and quietly paying that back is this: run it against build-release/ before and after, and
// compare.
//
// It is not a ctest. It reports numbers; it has no opinion about them, because the fixtures in
// this repository are far too small to set a threshold that would mean anything. Point it at a
// real dataset when the answer matters:
//
//     ./carta_zarr_pass_timing [dataset-path] [image-id] [threads] [repeats]
//
// Reported per entry point: the minimum over the repeats, which is the statistic that survives a
// noisy machine, and the median beside it so that a suspiciously lonely minimum is visible.

#include <carta-zarr/carta_zarr.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Timing {
    double best_ms = 0.0;
    double median_ms = 0.0;
    bool ok = false;
    std::string failure;
};

template <typename Body>
Timing TimeIt(unsigned int repeats, Body&& body) {
    std::vector<double> samples;
    samples.reserve(repeats);
    for (unsigned int i = 0; i < repeats; ++i) {
        const auto started = Clock::now();
        const auto outcome = body();
        const auto elapsed = Clock::now() - started;
        if (!outcome.empty()) {
            return Timing{0.0, 0.0, false, outcome};
        }
        samples.push_back(std::chrono::duration<double, std::milli>(elapsed).count());
    }
    std::sort(samples.begin(), samples.end());
    Timing timing;
    timing.ok = true;
    timing.best_ms = samples.front();
    timing.median_ms = samples.at(samples.size() / 2);
    return timing;
}

void Report(const std::string& label, const Timing& timing) {
    std::cout << std::left << std::setw(22) << label;
    if (!timing.ok) {
        std::cout << "FAILED  " << timing.failure << '\n';
        return;
    }
    std::cout << std::right << std::fixed << std::setprecision(3) << std::setw(10) << timing.best_ms
              << " ms   median " << std::setw(10) << timing.median_ms << " ms\n";
}

std::uint64_t AxisLength(const carta::zarr::ImageDescriptor& descriptor, carta::zarr::AxisRole role) {
    for (const auto& axis : descriptor.axes) {
        if (axis.role == role) {
            return axis.length;
        }
    }
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string fixture = argc > 1 ? argv[1] : CARTA_ZARR_PIXEL_FIXTURE_WIDE;
    const std::string image_id = argc > 2 ? argv[2] : "SKY";
    const unsigned int threads = argc > 3 ? static_cast<unsigned int>(std::strtoul(argv[3], nullptr, 10)) : 4;
    const unsigned int repeats = argc > 4 ? static_cast<unsigned int>(std::strtoul(argv[4], nullptr, 10)) : 9;

    if (!std::filesystem::exists(fixture)) {
        std::cerr << "no such dataset: " << fixture << '\n';
        return 2;
    }

    carta::zarr::OpenOptions open_options;
    open_options.decode_threads = threads;
    const auto context = carta::zarr::Context::Create(open_options);
    if (!context) {
        std::cerr << "Context::Create failed: " << context.error().message << '\n';
        return 2;
    }
    const auto dataset = carta::zarr::Dataset::Open(context.value(), fixture);
    if (!dataset) {
        std::cerr << "Dataset::Open failed: " << dataset.error().message << '\n';
        return 2;
    }
    const auto image = dataset.value().OpenImage(image_id);
    if (!image) {
        std::cerr << "OpenImage failed: " << image.error().message << '\n';
        return 2;
    }
    const auto& sky = image.value();
    const auto& descriptor = sky.descriptor();

    const auto x = AxisLength(descriptor, carta::zarr::AxisRole::spatial_x);
    const auto y = AxisLength(descriptor, carta::zarr::AxisRole::spatial_y);
    const auto channels = AxisLength(descriptor, carta::zarr::AxisRole::spectral);

    std::cout << "dataset  " << fixture << "\nimage    " << image_id << "  " << x << " x " << y << " x "
              << channels << " channels, " << threads << " threads, " << repeats << " repeats\n\n";

    // A whole-cube read of one polarization.
    const auto read = TimeIt(repeats, [&]() -> std::string {
        carta::zarr::ReadRequest request;
        request.axes.resize(descriptor.axes.size());
        std::size_t elements = 1;
        for (std::size_t i = 0; i < descriptor.axes.size(); ++i) {
            const auto& axis = descriptor.axes.at(i);
            const bool whole = axis.role == carta::zarr::AxisRole::spatial_x ||
                               axis.role == carta::zarr::AxisRole::spatial_y ||
                               axis.role == carta::zarr::AxisRole::spectral;
            request.axes.at(i) = {0, whole ? axis.length : 1, 1};
            elements *= request.axes.at(i).count;
        }
        std::vector<float> destination(elements);
        carta::zarr::BufferView<float> view{destination.data(), destination.size()};
        const auto outcome = sky.Read(request, view);
        return outcome ? std::string{} : outcome.error().message;
    });

    // Two overlapping regions over every channel, which is what the multi-region pass exists for.
    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, x, y, nullptr},
        {x / 4, y / 4, x / 2, y / 2, nullptr},
    };
    const auto reduce_over = [&](const std::vector<carta::zarr::RegionMask>& over) {
        return TimeIt(repeats, [&]() -> std::string {
            carta::zarr::SpectralReduceRequest request;
            request.planes.spectral = {0, channels, 1};
            request.regions = over.data();
            request.region_count = over.size();
            request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::nan_count |
                                 carta::zarr::Statistic::sum | carta::zarr::Statistic::sum_sq |
                                 carta::zarr::Statistic::min | carta::zarr::Statistic::max;
            const auto outcome = sky.ReduceSpectral(request, [](const carta::zarr::SpectralBlock&) { return true; });
            return outcome ? std::string{} : outcome.error().message;
        });
    };
    const auto reduce = reduce_over(regions);

    // The whole plane as one region: the case ADR 0005's figures are for.
    const auto reduce_whole = reduce_over({regions.front()});

    // Sixty-four one-pixel strips, half of them along x and half along y. The per-pixel loop is the
    // same as above; what this grows is the work done once per row of a region -- folding a row's
    // totals into its slots, above all -- which two large regions hardly exercise. A region's rows
    // run along whichever spatial axis the store varies fastest, so taking both orientations makes
    // half the strips a row per pixel whichever axis that is.
    std::vector<carta::zarr::RegionMask> boxes;
    for (std::uint64_t i = 0; i < 32; ++i) {
        boxes.push_back({(i * x) / 32, 0, 1, y, nullptr});
        boxes.push_back({0, (i * y) / 32, x, 1, nullptr});
    }
    const auto reduce_boxes = reduce_over(boxes);

    // One ellipse filling the plane, which is what a masked region is: once with runs made here
    // beforehand -- what a caller that keeps them between calls hands over -- and once as the raster
    // alone, which the reduction has to turn into runs itself on every call. The difference between
    // the two is what a caller gives up by not keeping its own.
    std::vector<std::uint8_t> ellipse(x * y, 0);
    for (std::uint64_t row = 0; row < y; ++row) {
        for (std::uint64_t column = 0; column < x; ++column) {
            const double dx = (static_cast<double>(column) + 0.5 - (x / 2.0)) / (x / 2.0);
            const double dy = (static_cast<double>(row) + 0.5 - (y / 2.0)) / (y / 2.0);
            ellipse.at((row * x) + column) = (dx * dx) + (dy * dy) <= 1.0 ? 1 : 0;
        }
    }
    const auto fastest = sky.chunk_geometry().fastest_spatial_axis;
    const bool along_y = fastest == carta::zarr::AxisRole::spatial_y;
    std::vector<std::uint32_t> runs;
    std::vector<std::uint64_t> run_offsets{0};
    for (std::uint64_t line = 0; line < (along_y ? x : y); ++line) {
        const std::uint64_t length = along_y ? y : x;
        const auto at = [&](std::uint64_t i) {
            return ellipse.at(along_y ? (i * x) + line : (line * x) + i) != 0;
        };
        for (std::uint64_t i = 0; i < length;) {
            if (!at(i)) {
                ++i;
                continue;
            }
            const auto begin = i;
            while (i < length && at(i)) {
                ++i;
            }
            runs.push_back(static_cast<std::uint32_t>(begin));
            runs.push_back(static_cast<std::uint32_t>(i));
        }
        run_offsets.push_back(runs.size() / 2);
    }
    carta::zarr::RegionMask as_raster{0, 0, x, y, ellipse.data()};
    carta::zarr::RegionMask as_runs = as_raster;
    as_runs.row_runs = runs.data();
    as_runs.row_run_offsets = run_offsets.data();
    as_runs.run_axis = fastest;
    const auto reduce_runs = reduce_over({as_runs});
    const auto reduce_raster = reduce_over({as_raster});

    const auto histogram = TimeIt(repeats, [&]() -> std::string {
        carta::zarr::HistogramRequest request;
        request.planes.spectral = {0, channels, 1};
        request.bins = 1024;
        request.lower = -1.0e9;
        request.upper = 1.0e9;
        const auto outcome = sky.ComputeHistogram(request, [](const carta::zarr::HistogramBlock&) { return true; });
        return outcome ? std::string{} : outcome.error().message;
    });

    const auto cube = TimeIt(repeats, [&]() -> std::string {
        carta::zarr::CubeHistogramRequest request;
        request.planes.spectral = {0, channels, 1};
        request.bins = 1024;
        const auto outcome = sky.ComputeCubeHistogram(request);
        return outcome ? std::string{} : outcome.error().message;
    });

    Report("Read", read);
    Report("ReduceSpectral", reduce);
    Report("ReduceSpectral x1", reduce_whole);
    Report("ReduceSpectral x64", reduce_boxes);
    Report("ReduceSpectral runs", reduce_runs);
    Report("ReduceSpectral raster", reduce_raster);
    Report("ComputeHistogram", histogram);
    Report("ComputeCubeHistogram", cube);

    const bool all_ok = read.ok && reduce.ok && reduce_whole.ok && reduce_boxes.ok && reduce_runs.ok &&
                        reduce_raster.ok && histogram.ok && cube.ok;
    return all_ok ? 0 : 1;
}
