/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_EXAMPLES_SUPPORT_H_
#define CARTA_ZARR_EXAMPLES_SUPPORT_H_

#include <algorithm>
#include <iostream>
#include <limits>

#include <carta-zarr/carta_zarr.h>

namespace examples {

inline int Report(const carta::zarr::Error& error) {
    std::cerr << carta::zarr::ErrorCodeName(error.code) << ": " << error.message;
    if (!error.node_path.empty()) {
        std::cerr << " [" << error.node_path << ']';
    }
    std::cerr << '\n';
    return 1;
}

inline carta::zarr::Result<carta::zarr::Image> OpenFirstImage(const carta::zarr::Context& context, const char* path) {
    using namespace carta::zarr;
    const auto dataset = Dataset::Open(context, path);
    if (!dataset) {
        return dataset.error();
    }
    for (const auto& entry : dataset.value().descriptor().images) {
        if (entry.openable) {
            // Image retains the Store and Context after this local Dataset ends.
            return dataset.value().OpenImage(entry.id);
        }
    }
    return Error{ErrorCode::invalid_metadata, "No openable image was listed", path};
}

// A bounded spatial patch at the first channel, polarization and time.
inline carta::zarr::Result<carta::zarr::ReadRequest> Patch(const carta::zarr::Image& image) {
    using namespace carta::zarr;
    const auto& axes = image.descriptor().axes;
    const auto x = AxisIndex(axes, AxisRole::spatial_x);
    const auto y = AxisIndex(axes, AxisRole::spatial_y);
    const auto spectral = AxisIndex(axes, AxisRole::spectral);
    if (!x || !y || !spectral) {
        return Error{ErrorCode::invalid_argument, "The example needs spatial and spectral axes", image.descriptor().id};
    }
    for (const auto& axis : axes) {
        if (axis.length == 0) {
            return Error{ErrorCode::invalid_slice, "An axis is empty", axis.name};
        }
    }
    ReadRequest request;
    request.axes.assign(axes.size(), {0, 1, 1});
    request.axes[*x].count = std::min<std::uint64_t>(64, axes[*x].length);
    request.axes[*y].count = std::min<std::uint64_t>(64, axes[*y].length);
    return request;
}

inline carta::zarr::Result<std::size_t> Elements(const carta::zarr::ReadRequest& request) {
    using namespace carta::zarr;
    std::size_t elements = 1;
    for (const auto& range : request.axes) {
        if (range.count == 0 || range.count > std::numeric_limits<std::size_t>::max() / elements) {
            return Error{ErrorCode::invalid_argument, "The selected shape is empty or too large", {}};
        }
        elements *= static_cast<std::size_t>(range.count);
    }
    return elements;
}

}  // namespace examples

#endif  // CARTA_ZARR_EXAMPLES_SUPPORT_H_
