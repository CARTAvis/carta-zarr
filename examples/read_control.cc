/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "support.h"

#include <atomic>
#include <chrono>
#include <vector>

int main(int argc, char** argv) {
    using namespace carta::zarr;
    if (argc != 2) {
        std::cerr << "Usage: carta_zarr_read_control PATH_TO_DATASET.zarr\n";
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
    auto patch = examples::Patch(image);
    if (!patch) {
        return examples::Report(patch.error());
    }
    auto& request = patch.value();
    const auto& axes = image.descriptor().axes;
    const auto spectral = *AxisIndex(axes, AxisRole::spectral);
    request.axes[spectral].start = axes[spectral].length > 1 ? 1 : 0;
    const auto elements = examples::Elements(request);
    if (!elements) {
        return examples::Report(elements.error());
    }
    std::vector<float> pixels(elements.value());

    // A one-off scan uses its own zero-byte pool, leaving the session's cache
    // available for interactive reads. Cancellation callbacks may use an atomic
    // flag set by the application's event loop or another thread.
    const auto pool = context.value().CreateCachePool(0);
    if (!pool) {
        return examples::Report(pool.error());
    }
    std::atomic<bool> stop{false};
    ReadOptions options;
    options.control.cache_pool = pool.value();
    options.control.cancellation_requested = [&stop] { return stop.load(); };
    options.control.deadline = std::chrono::steady_clock::now() + std::chrono::seconds(30);
    options.read_budget_bytes = std::size_t{16} << 20;
    std::size_t completed = 0;
    const auto read =
        image.Read(request, {pixels.data(), pixels.size()}, options, [&](std::size_t written, std::size_t total) {
            completed = written;
            std::cout << "progress " << written << '/' << total << '\n';
            return true;
        });
    if (!read) {
        // Only what progress last reported is final; this example keeps none of it.
        return examples::Report(read.error());
    }
    if (read.value() != elements.value() || completed != elements.value()) {
        std::cerr << "The read did not report its complete output\n";
        return 1;
    }

    // Demonstrate cancellation deterministically, without racing a fast read.
    stop = true;
    const auto cancelled = image.Read(request, {pixels.data(), pixels.size()}, options);
    if (cancelled || cancelled.error().code != ErrorCode::cancelled) {
        std::cerr << "The stopped read did not report cancellation\n";
        return 1;
    }
    std::cout << "Cancellation reported; the destination is discarded\n";

    // Prefetch must retain chunks to benefit a following read. It uses the
    // context's shared cache here, rather than the scan's zero-byte pool.
    const auto prefetched = image.Prefetch(request);
    if (!prefetched) {
        return examples::Report(prefetched.error());
    }
    const auto reread = image.Read(request, {pixels.data(), pixels.size()});
    if (!reread) {
        return examples::Report(reread.error());
    }
    std::cout << "Prefetched " << prefetched.value() << " chunks; read channel " << request.axes[spectral].start
              << '\n';
    return reread.value() == elements.value() ? 0 : 1;
}
