/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <carta-zarr/carta_zarr.h>

namespace {

int Report(const carta::zarr::Error& error) {
    std::cerr << carta::zarr::ErrorCodeName(error.code) << ": " << error.message;
    if (!error.node_path.empty()) {
        std::cerr << " [" << error.node_path << ']';
    }
    std::cerr << '\n';
    return 1;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace carta::zarr;
    if (argc != 2) {
        std::cerr << "Usage: carta_zarr_first_read PATH_TO_DATASET.zarr\n";
        return 1;
    }

    const auto context = Context::Create();
    if (!context) {
        return Report(context.error());
    }
    const auto dataset = Dataset::Open(context.value(), argv[1]);
    if (!dataset) {
        return Report(dataset.error());
    }

    const auto& descriptor = dataset.value().descriptor();
    std::cout << "schema=" << descriptor.schema_id << " version=" << descriptor.schema_version << '\n';
    for (const auto& diagnostic : descriptor.diagnostics) {
        std::cout << "  " << DiagnosticCodeName(diagnostic.code) << ": " << diagnostic.message << '\n';
    }
    const auto& entries = descriptor.images;
    for (const auto& entry : entries) {
        std::cout << entry.id << " openable=" << entry.openable << " role=" << entry.image_role << '\n';
        for (const auto& diagnostic : entry.diagnostics) {
            std::cout << "  " << DiagnosticCodeName(diagnostic.code) << ": " << diagnostic.message << '\n';
        }
    }
    const auto chosen =
        std::find_if(entries.begin(), entries.end(), [](const ImageEntry& entry) { return entry.openable; });
    if (chosen == entries.end()) {
        std::cerr << "No openable image was listed\n";
        return 1;
    }

    const auto image = dataset.value().OpenImage(chosen->id);
    if (!image) {
        return Report(image.error());
    }
    const auto& image_descriptor = image.value().descriptor();
    if (!image_descriptor.sub_type.empty()) {
        std::cout << "sub_type=" << image_descriptor.sub_type << '\n';
    }
    if (image_descriptor.spectral) {
        const auto& spectral = *image_descriptor.spectral;
        std::cout << "spectral frame=" << spectral.system;
        if (spectral.increment) {
            std::cout << " increment=" << *spectral.increment << ' ' << spectral.unit;
        }
        std::cout << '\n';
    }
    const auto& axes = image_descriptor.axes;
    const auto x = AxisIndex(axes, AxisRole::spatial_x);
    const auto y = AxisIndex(axes, AxisRole::spatial_y);
    if (!x || !y) {
        std::cerr << "The image has no pair of spatial axes\n";
        return 1;
    }
    for (const auto& axis : axes) {
        if (axis.length == 0) {
            std::cerr << "An axis is empty: " << axis.name << '\n';
            return 1;
        }
    }

    ReadRequest request;
    request.axes.resize(axes.size(), {0, 1, 1});
    request.axes[*x].count = std::min<std::uint64_t>(64, axes[*x].length);
    request.axes[*y].count = std::min<std::uint64_t>(64, axes[*y].length);

    std::size_t elements = 1;
    for (const auto& range : request.axes) {
        if (range.count > std::numeric_limits<std::size_t>::max() / elements) {
            std::cerr << "The selected shape does not fit in memory\n";
            return 1;
        }
        elements *= static_cast<std::size_t>(range.count);
    }
    std::vector<float> pixels(elements);
    const auto read = image.value().Read(request, {pixels.data(), pixels.size()});
    if (!read) {
        // A failed read may have written part of the buffer. Do not consume it.
        return Report(read.error());
    }
    const auto finite = std::count_if(pixels.begin(), pixels.end(), [](float pixel) { return std::isfinite(pixel); });
    std::cout << "Read " << read.value() << " pixels from " << chosen->id << "; " << finite << " are finite\n";
    return read.value() == elements ? 0 : 1;
}
