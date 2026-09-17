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
    const auto reduce = TimeIt(repeats, [&]() -> std::string {
        carta::zarr::SpectralReduceRequest request;
        request.planes.spectral = {0, channels, 1};
        request.regions = regions.data();
        request.region_count = regions.size();
        request.statistics = carta::zarr::Statistic::num_pixels | carta::zarr::Statistic::nan_count |
                             carta::zarr::Statistic::sum | carta::zarr::Statistic::sum_sq |
                             carta::zarr::Statistic::min | carta::zarr::Statistic::max;
        const auto outcome = sky.ReduceSpectral(request, [](const carta::zarr::SpectralBlock&) { return true; });
        return outcome ? std::string{} : outcome.error().message;
    });

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
    Report("ComputeHistogram", histogram);
    Report("ComputeCubeHistogram", cube);

    const bool all_ok = read.ok && reduce.ok && histogram.ok && cube.ok;
    return all_ok ? 0 : 1;
}
