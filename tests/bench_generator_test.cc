/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What tools/zarr-bench/generate.py writes, read back through the library.
//
// The generator exists so that layouts can be compared, and a comparison means something only when
// the layouts hold the same pixels and the library reads each of them as the layout it was asked to
// be. So this pins both: three layouts of one synthetic cube read back identical, each reports the
// chunk and shard shape its manifest says it was written with, and a rewritten fixture reads back as
// the window of the original it was cut from. tests/run_bench_generator.cmake writes the datasets.

#include <carta-zarr/carta_zarr.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::AxisRole;
using carta::zarr::testing::Require;

const std::string kGenerated = CARTA_ZARR_GENERATED_DIR;

carta::zarr::Context SharedContext() {
    static const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    return context.value();
}

carta::zarr::Image OpenDefault(const std::string& location) {
    const auto probe = carta::zarr::ProbeSchema(location, carta::zarr::kXradioImageSchema);
    Require(probe.has_value() && probe->kind == carta::zarr::SchemaMatchKind::match,
            location + " is not an XRADIO image dataset to the schema probe");
    const auto dataset = carta::zarr::Dataset::Open(SharedContext(), location);
    Require(dataset.has_value(), location + " did not open: " + (dataset ? "" : dataset.error().message));
    const auto& id = dataset->descriptor().default_image_id;
    Require(id.has_value(), location + " lists no image that opens");
    auto image = dataset->OpenImage(*id);
    Require(image.has_value(), location + ": " + *id + " did not open");
    return image.value();
}

std::vector<float> ReadAll(const carta::zarr::Image& image, bool apply_pixel_mask = true) {
    carta::zarr::ReadRequest request;
    std::size_t elements = 1;
    for (const auto& axis : image.descriptor().axes) {
        request.axes.push_back({0, axis.length, 1});
        elements *= axis.length;
    }
    std::vector<float> pixels(elements);
    carta::zarr::ReadOptions options;
    options.apply_pixel_mask = apply_pixel_mask;
    const auto written = image.Read(request, {pixels.data(), pixels.size()}, options);
    Require(written.has_value() && written.value() == elements, "a whole-image read did not fill the buffer");
    return pixels;
}

// Equal as stored bits, so that NaN equals NaN and -0 does not equal 0: the generator promises the
// same values, not merely values that compare equal.
bool SameBits(const std::vector<float>& left, const std::vector<float>& right) {
    return left.size() == right.size() && std::memcmp(left.data(), right.data(), left.size() * sizeof(float)) == 0;
}

nlohmann::json Manifest(const std::string& location) {
    std::ifstream file(location + "/bench-manifest.json");
    Require(file.good(), location + " has no manifest");
    auto manifest = nlohmann::json::parse(file);
    Require(manifest.at("complete").get<bool>(), location + "'s manifest does not say it is complete");
    return manifest;
}

// The manifest's shape for each axis, by name, against the library's for the same axis.
void RequireGeometryMatches(const std::string& location, const carta::zarr::Image& image) {
    const auto manifest = Manifest(location);
    const auto& dims = manifest.at("image").at("dimension_names");
    const auto& layout = manifest.at("layout");
    const auto& geometry = image.chunk_geometry();
    const auto& axes = image.descriptor().axes;
    for (std::size_t logical = 0; logical < axes.size(); ++logical) {
        std::size_t stored = dims.size();
        for (std::size_t index = 0; index < dims.size(); ++index) {
            if (dims.at(index).get<std::string>() == axes.at(logical).name) {
                stored = index;
            }
        }
        Require(stored < dims.size(), location + ": axis " + axes.at(logical).name + " is not in the manifest");
        const auto chunk = layout.at("chunk_shape").at(stored).get<std::uint64_t>();
        const auto shard =
            layout.at("shard_shape").is_null() ? chunk : layout.at("shard_shape").at(stored).get<std::uint64_t>();
        Require(geometry.chunk_shape.at(logical) == chunk,
                location + ": " + axes.at(logical).name + " is not chunked as the manifest says");
        Require(geometry.shard_shape.at(logical) == shard,
                location + ": " + axes.at(logical).name + " is not sharded as the manifest says");
    }
}

void TestLayoutsHoldTheSamePixels() {
    const auto plain = OpenDefault(kGenerated + "/plain");
    const auto sharded = OpenDefault(kGenerated + "/sharded");
    RequireGeometryMatches(kGenerated + "/plain", plain);
    RequireGeometryMatches(kGenerated + "/sharded", sharded);
    Require(sharded.chunk_geometry().sharded && !plain.chunk_geometry().sharded,
            "the sharded layout was not read as sharded, or the plain one was");

    const auto pixels = ReadAll(plain);
    Require(SameBits(pixels, ReadAll(sharded)), "two layouts of one synthetic cube hold different pixels");

    // Not a cube of fill: there is noise, there are NaN corners outside the circle, and the centre is
    // finite. Axis 0 is l and fastest, so the first element is a corner.
    std::size_t finite = 0;
    for (const float value : pixels) {
        finite += std::isfinite(value) ? 1 : 0;
    }
    Require(std::isnan(pixels.front()) && std::isnan(pixels.back()), "the corners are not NaN");
    Require(finite > pixels.size() / 2 && finite < pixels.size(), "the cube is not mostly finite with a NaN border");
}

// A flag where the plain layout writes NaN: masked, the two read the same; unmasked, the flagged
// pixels hold noise rather than NaN.
void TestAFlagMasksWhatNaNWouldHave() {
    const auto plain = OpenDefault(kGenerated + "/plain");
    const auto flagged = OpenDefault(kGenerated + "/flagged");
    Require(flagged.descriptor().has_pixel_mask, "the flagged layout has no pixel mask");
    Require(SameBits(ReadAll(plain), ReadAll(flagged)), "a masked read of the flagged layout differs from NaN");
    const auto unmasked = ReadAll(flagged, false);
    Require(std::isfinite(unmasked.front()), "an unmasked read of a flagged pixel is not the pixel");
}

// The rewrite of the committed pixel fixture, cropped to l in [1, 4), against that window of the
// original. Both are read masked, so the flag was carried across with the pixels.
void TestARewriteIsTheWindowItWasCutFrom() {
    const auto original = OpenDefault(CARTA_ZARR_PIXEL_FIXTURE);
    const auto rewritten = OpenDefault(kGenerated + "/rewritten");
    RequireGeometryMatches(kGenerated + "/rewritten", rewritten);
    Require(rewritten.descriptor().has_pixel_mask, "the rewrite lost the fixture's flag");

    const auto& axes = original.descriptor().axes;
    carta::zarr::ReadRequest window;
    std::size_t elements = 1;
    for (const auto& axis : axes) {
        window.axes.push_back({0, axis.length, 1});
    }
    const auto l = *carta::zarr::AxisIndex(axes, AxisRole::spatial_x);
    window.axes.at(l) = {1, 3, 1};
    for (const auto& range : window.axes) {
        elements *= range.count;
    }
    std::vector<float> expected(elements);
    Require(original.Read(window, {expected.data(), expected.size()}).has_value(), "the fixture's window did not read");
    Require(SameBits(expected, ReadAll(rewritten)), "the rewrite is not the window it was cut from");
}

// The source's own layout, taken by --layout-from-source: the fixture's chunks and codec, and its
// pixels over the crop.
void TestTheSourcesOwnLayoutIsTheSources() {
    const auto original = OpenDefault(CARTA_ZARR_PIXEL_FIXTURE_WIDE);
    const auto current = OpenDefault(kGenerated + "/current");
    RequireGeometryMatches(kGenerated + "/current", current);
    Require(current.chunk_geometry().chunk_shape == original.chunk_geometry().chunk_shape,
            "the source's own layout is not chunked as the source is");
    Require(current.chunk_geometry().compressor == original.chunk_geometry().compressor &&
                Manifest(kGenerated + "/current").at("layout").at("codec") == "zstd:9",
            "the source's own layout is not compressed as the source is");

    const auto& axes = original.descriptor().axes;
    carta::zarr::ReadRequest window;
    std::size_t elements = 1;
    for (const auto& axis : axes) {
        window.axes.push_back({0, axis.length, 1});
    }
    window.axes.at(*carta::zarr::AxisIndex(axes, AxisRole::spectral)) = {1, 2, 1};
    for (const auto& range : window.axes) {
        elements *= range.count;
    }
    std::vector<float> expected(elements);
    Require(original.Read(window, {expected.data(), expected.size()}).has_value(), "the fixture's window did not read");
    Require(SameBits(expected, ReadAll(current)), "the source's own layout is not the window it was cut from");
}

}  // namespace

int main() {
    try {
        TestLayoutsHoldTheSamePixels();
        TestAFlagMasksWhatNaNWouldHave();
        TestARewriteIsTheWindowItWasCutFrom();
        TestTheSourcesOwnLayoutIsTheSources();
    } catch (const std::exception& error) {
        std::cerr << "bench generator: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
