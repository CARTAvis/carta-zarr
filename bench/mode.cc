/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "mode.h"

#include <utility>

namespace carta::zarr::bench {

std::unique_ptr<Workload> Workload::For(Mode mode, const RunOptions& options) {
    switch (mode) {
        case Mode::plane:
            return PlaneWorkload(options);
        case Mode::spectrum:
            return SpectrumWorkload(options);
        case Mode::region:
            return RegionWorkload(options, options.region);
        case Mode::cube_histogram:
            return CubeHistogramWorkload(options, options.histogram);
        case Mode::open:
            return OpenWorkload(options);
        case Mode::animation:
            return AnimationWorkload(options, options.animation);
    }
    return nullptr;
}

Result<Planned> Workload::Begin(unsigned trial, unsigned process_index, Row& row) const {
    const auto failed = [](const std::string& what, const Error& error) {
        return Error{error.code, what + bench::Describe(error), error.node_path};
    };
    auto context = Context::Create(_options.context);
    if (!context) {
        return failed("Context::Create: ", context.error());
    }
    auto dataset = Dataset::Open(*context, _options.dataset);
    if (!dataset) {
        return failed("Dataset::Open: ", dataset.error());
    }
    auto id = _options.image_id;
    if (id.empty()) {
        if (!dataset->descriptor().default_image_id) {
            return Error{ErrorCode::not_found, "the dataset lists no image that opens", _options.dataset};
        }
        id = *dataset->descriptor().default_image_id;
    }
    auto image = dataset->OpenImage(id);
    if (!image) {
        return failed("OpenImage: ", image.error());
    }
    const auto axes = CubeAxes::Of(image->descriptor());
    if (!axes) {
        return failed("", axes.error());
    }
    DescribeImage(row, *image);

    // Every process's plan, to tell which of this one's operations share a chunk with another's.
    const unsigned ops = _options.OpsFor(_mode);
    std::vector<std::vector<Operation>> plans;
    for (unsigned process = 0; process < _options.processes; ++process) {
        plans.push_back(Plan(*axes, PlanSeed{_options.seed, trial, _options.processes, process}, ops));
    }
    MarkSharedChunks(*this, plans, *axes, image->chunk_geometry().chunk_shape);

    auto runner = MakeRunner(*context, std::move(image).value());
    if (!runner) {
        return failed("", runner.error());
    }
    return Planned{std::move(plans[process_index]), std::move(runner).value()};
}

void MarkSharedChunks(const Workload& workload, std::vector<std::vector<Operation>>& plans, const CubeAxes& axes,
                      const std::vector<std::uint64_t>& chunk_shape) {
    struct Placed {
        std::size_t process;
        std::size_t index;
        ChunkBox box;
    };
    std::vector<Placed> placed;
    for (std::size_t process = 0; process < plans.size(); ++process) {
        for (std::size_t index = 0; index < plans[process].size(); ++index) {
            if (const auto box = workload.ChunksRead(plans[process][index], axes, chunk_shape)) {
                placed.push_back({process, index, *box});
            }
        }
    }
    for (const auto& one : placed) {
        bool shares = false;
        for (const auto& other : placed) {
            const bool earlier_here = other.process == one.process && other.index < one.index;
            if ((earlier_here || other.process != one.process) && one.box.Meets(other.box)) {
                shares = true;
                break;
            }
        }
        plans[one.process][one.index].shares_chunks = shares;
    }
}

}  // namespace carta::zarr::bench
