/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "support.h"

#include <vector>

#include <carta-zarr/read_ahead.h>

int main(int argc, char** argv) {
    using namespace carta::zarr;
    if (argc != 2) {
        std::cerr << "Usage: carta_zarr_animation PATH_TO_DATASET.zarr\n";
        return 1;
    }
    ContextOptions context_options;
    context_options.cache_bytes = std::size_t{64} << 20;
    const auto context = Context::Create(context_options);
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
    const auto elements = examples::Elements(patch.value());
    if (!elements) {
        return examples::Report(elements.error());
    }
    const auto& axes = image.descriptor().axes;
    const auto spectral = *AxisIndex(axes, AxisRole::spectral);
    const auto frames = std::min<std::uint64_t>(4, axes[spectral].length);
    const auto frame_request = [&](std::uint64_t channel) {
        auto request = patch.value();
        request.axes[spectral].start = channel;
        return request;
    };

    ReadOptions options;
    auto made = ReadAhead::For({{image, options}});
    if (!made) {
        // A large dataset can need more than this example's 64 MiB cache.
        // Report the required capacity instead of silently disabling read-ahead.
        return examples::Report(made.error());
    }
    auto ahead = made.value();
    std::vector<float> pixels(elements.value());
    for (std::uint64_t channel = 0; channel < frames; ++channel) {
        const auto began = ReadAhead::Clock::now();
        const auto request = frame_request(channel);
        const auto read = image.Read(request, {pixels.data(), pixels.size()}, options);
        if (!read) {
            return examples::Report(read.error());
        }
        if (read.value() != elements.value()) {
            std::cerr << "The frame is incomplete\n";
            return 1;
        }
        std::cout << "frame channel=" << channel << " pixels=" << read.value() << '\n';
        // A real viewer displays pixels here and measures lateness against its
        // frame schedule. This CLI shows the API lifecycle, not a speed benchmark.
        std::vector<std::vector<AnimatedPlane>> upcoming;
        for (std::uint64_t next = channel + 1; next < frames; ++next) {
            upcoming.push_back({AnimatedPlane{0, frame_request(next)}});
        }
        ahead.Served(began, false, {AnimatedPlane{0, request}}, upcoming);
    }
    ahead.Cancel();
    const auto stats = ahead.stats();
    std::cout << "Read-ahead prefetches=" << stats.prefetches << " stopped=" << stats.stopped << '\n';
    // The last controller is destroyed before image/resources leave scope.
    return stats.stopped ? 0 : 1;
}
