/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// open: Context::Create, Dataset::Open and OpenImage, with nothing cached. The one mode that reads no
// cube: each operation brings its own context and image, so a process opens nothing before the clock
// starts, and its rows say what was opened only once something has been.

#include "../mode.h"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace carta::zarr::bench {

namespace {

class OpenRunner final : public Runner {
public:
    OpenRunner(ContextOptions context, std::string dataset, std::string image_id)
        : _context(std::move(context)), _dataset(std::move(dataset)), _image_id(std::move(image_id)) {}

    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options) override {
        (void)operation;
        (void)options;
        _lengths.clear();
        auto context = Context::Create(_context);
        if (!context) {
            return std::move(context).error();
        }
        auto dataset = Dataset::Open(*context, _dataset);
        if (!dataset) {
            return std::move(dataset).error();
        }
        std::string id = _image_id;
        if (id.empty()) {
            if (!dataset->descriptor().default_image_id) {
                return Error{ErrorCode::not_found, "the dataset lists no image that opens", _dataset};
            }
            id = *dataset->descriptor().default_image_id;
        }
        auto image = dataset->OpenImage(id);
        if (!image) {
            return std::move(image).error();
        }
        _image = std::move(image).value();
        for (const auto& axis : _image->descriptor().axes) {
            _lengths.push_back(static_cast<double>(axis.length));
        }
        return std::uint64_t{0};
    }

    // The lengths of the axes it opened.
    std::uint64_t Fingerprint() const override {
        return bench::Fingerprint(_lengths.data(), _lengths.size());
    }

    // The image last opened, which an open that failed leaves as it was.
    void Record(Row& result) const override {
        std::size_t item_size = 4;
        if (_image) {
            DescribeImage(result, *_image);
            item_size = ItemSize(_image->descriptor().stored_type);
        }
        result.logical_bytes = result.elements * item_size;
    }

private:
    ContextOptions _context;
    std::string _dataset;
    std::string _image_id;
    std::optional<Image> _image;
    std::vector<double> _lengths;
};

class Open final : public Workload {
public:
    explicit Open(const RunOptions& options) : Workload(Mode::open, options) {}

    std::vector<Operation> Plan(const CubeAxes& axes, const PlanSeed& at, unsigned ops) const override {
        (void)axes;
        (void)at;
        return std::vector<Operation>(ops);
    }

    std::string Describe(const Operation& operation) const override {
        (void)operation;
        return "";
    }

    Result<Planned> Begin(unsigned trial, unsigned process_index, Row& row) const override {
        (void)row;
        const PlanSeed at{options().seed, trial, options().processes, process_index};
        return Planned{Plan(CubeAxes{}, at, options().OpsFor(mode())), NewRunner()};
    }

    // Each operation opens a context and an image of its own, so the ones it is handed go unused.
    Result<std::unique_ptr<Runner>> MakeRunner(const Context& context, Image image) const override {
        (void)context;
        (void)image;
        return NewRunner();
    }

private:
    std::unique_ptr<Runner> NewRunner() const {
        return std::make_unique<OpenRunner>(options().context, options().dataset, options().image_id);
    }
};

}  // namespace

std::unique_ptr<Workload> OpenWorkload(const RunOptions& options) {
    return std::make_unique<Open>(options);
}

}  // namespace carta::zarr::bench
