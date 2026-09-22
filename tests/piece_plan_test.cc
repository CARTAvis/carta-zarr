/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What an ordinary read decides before it reads anything: whether to cut itself into pieces, where,
// and how big each one should be.
//
// Until this existed those decisions could only be observed through a read against a directory tree
// -- a test would count progress callbacks and reason backwards to how the read must have been cut.
// None of this needs a store, a transport or a fixture, which is the point.

#include "read/pieces.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ChunkGeometry;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadOptions;
using carta::zarr::ReadRequest;
using carta::zarr::internal::PiecePlan;
using carta::zarr::internal::PlanPieces;

using carta::zarr::testing::Require;

// An image in the logical order XRADIO reports: x, y, spectral, polarization, time.
ImageDescriptor MakeImage(std::uint64_t x, std::uint64_t y, std::uint64_t channels, bool has_mask = false) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    descriptor.has_pixel_mask = has_mask;
    descriptor.pixel_mask_id = has_mask ? "FLAG" : "";
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, x},
             {"m", AxisRole::spatial_y, y},
             {"frequency", AxisRole::spectral, channels},
             {"polarization", AxisRole::polarization, 1},
             {"time", AxisRole::time, 1}};
    for (std::size_t i = 0; i < 5; ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = axes[i].name;
        axis.role = axes[i].role;
        axis.length = axes[i].length;
        axis.storage_index = std::vector<std::size_t>{3, 4, 1, 2, 0}.at(i);
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

ChunkGeometry MakeGeometry(std::uint64_t chunk_x, std::uint64_t chunk_y, std::uint64_t chunk_z) {
    ChunkGeometry geometry;
    geometry.chunk_shape = {chunk_x, chunk_y, chunk_z, 1, 1};
    return geometry;
}

ReadRequest WholeImage(const ImageDescriptor& descriptor) {
    ReadRequest request;
    for (const auto& axis : descriptor.axes) {
        request.axes.push_back(Range{0, axis.length, 1});
    }
    return request;
}

std::uint64_t ElementCount(const ReadRequest& request) {
    std::uint64_t elements = 1;
    for (const auto& range : request.axes) {
        elements *= range.count;
    }
    return elements;
}

PiecePlan Plan(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const ReadRequest& request,
               const ReadOptions& options, bool apply_mask = false, bool watching = false) {
    return PlanPieces(descriptor, geometry, request, options, watching, ElementCount(request), apply_mask);
}

// A read nobody is watching and nobody has put a ceiling on is issued exactly as it was asked for.
// The plan still describes it, as one piece covering everything, so the loop that reads it has one
// shape rather than two.
void TestAnUnconstrainedReadIsOnePiece() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    const auto request = WholeImage(image);

    const auto plan = Plan(image, geometry, request, {});
    Require(!plan.split, "a read with no reason to be cut was cut");
    Require(plan.units == 1, "an unsplit read should be one unit");
    Require(plan.elements_per_unit == ElementCount(request), "that one unit should be the whole destination");
}

// Either reason on its own is enough. The memory ceiling is the half that used to be ignored unless
// a progress callback came with it, which meant a caller who said how much memory the read could
// have and did not care to watch was refused instead of served.
void TestEitherReasonCutsTheRead() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    const auto request = WholeImage(image);

    Require(Plan(image, geometry, request, ReadOptions{}, false, /*watching=*/true).split,
            "a watched read was not cut");

    ReadOptions bounded;
    bounded.temporary_memory_limit_bytes = 4096;
    Require(Plan(image, geometry, request, bounded).split, "a read with a memory ceiling was not cut");
}

// The cut goes on the slowest-varying axis that selects more than one element, because the
// destination is dense with axis 0 fastest and that is the only axis whose pieces extend a prefix.
void TestTheCutGoesOnTheSlowestSelectedAxis() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    // Axes 3 and 4 are degenerate here, so the spectrum at index 2 is the slowest one selected.
    const auto plan = Plan(image, geometry, WholeImage(image), ReadOptions{}, false, /*watching=*/true);
    Require(plan.split && plan.axis == 2, "the cut should have gone on the spectral axis");
    Require(plan.units == 8, "the plan should cover every channel");
    Require(plan.elements_per_unit == 64 * 64, "one channel is worth a plane of the destination");

    // Narrow the spectrum to one channel and the only axis left with more than one element is m.
    auto one_channel = WholeImage(image);
    one_channel.axes.at(2) = Range{0, 1, 1};
    const auto narrowed = Plan(image, geometry, one_channel, ReadOptions{}, false, /*watching=*/true);
    Require(narrowed.split && narrowed.axis == 1, "with one channel the cut should move to m");
    Require(narrowed.elements_per_unit == 64, "one m is worth a row of the destination");
}

// A request that selects a single element of every axis has nowhere to be cut, so it is one piece
// whatever the caller asked for. A read that cannot be cut is where buffer_too_small comes from.
void TestAReadWithNowhereToCutIsOnePiece() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    ReadRequest single;
    single.axes.assign(5, Range{0, 1, 1});

    ReadOptions bounded;
    bounded.temporary_memory_limit_bytes = 1;
    Require(!Plan(image, geometry, single, bounded, false, /*watching=*/true).split,
            "a single-element read has nowhere to be cut");
}

// A piece is measured in chunks, not in elements: the budget buys whole chunks along the cut axis
// because asking for part of one decodes all of it anyway.
void TestATighterCeilingBuysFewerChunks() {
    const auto image = MakeImage(64, 64, 32);
    const auto geometry = MakeGeometry(64, 64, 4);
    const auto request = WholeImage(image);

    ReadOptions generous;
    generous.temporary_memory_limit_bytes = 64 * 64 * 4 * sizeof(float) * 8;
    ReadOptions tight;
    tight.temporary_memory_limit_bytes = 64 * 64 * 4 * sizeof(float);

    const auto wide = Plan(image, geometry, request, generous);
    const auto narrow = Plan(image, geometry, request, tight);
    Require(wide.chunk == 4 && narrow.chunk == 4, "the plan should report the chunk extent along the cut axis");
    Require(narrow.units_per_piece < wide.units_per_piece, "a tighter ceiling should buy a smaller piece");
    Require(narrow.units_per_piece >= 1, "a piece is never smaller than one unit");
}

// The flag is decoded beside the pixels, so a read that will apply it costs more per chunk and the
// same ceiling has to buy less. Counting it on one side of the division and not the other would
// size pieces against a cost the read does not have.
void TestApplyingTheFlagCostsTheBudget() {
    const auto image = MakeImage(64, 64, 32, /*has_mask=*/true);
    const auto geometry = MakeGeometry(64, 64, 1);
    const auto request = WholeImage(image);

    ReadOptions options;
    options.temporary_memory_limit_bytes = 64 * 64 * sizeof(float) * 6;

    const auto without = Plan(image, geometry, request, options, /*apply_mask=*/false);
    const auto with = Plan(image, geometry, request, options, /*apply_mask=*/true);
    Require(with.units_per_piece <= without.units_per_piece,
            "a read that also decodes the flag cannot afford as much of the spectrum per piece");
    Require(with.units_per_piece < without.units_per_piece,
            "with this budget the flag should cost at least one channel of the piece");
}

}  // namespace

int main() {
    try {
        TestAnUnconstrainedReadIsOnePiece();
        TestEitherReasonCutsTheRead();
        TestTheCutGoesOnTheSlowestSelectedAxis();
        TestAReadWithNowhereToCutIsOnePiece();
        TestATighterCeilingBuysFewerChunks();
        TestApplyingTheFlagCostsTheBudget();
    } catch (const std::exception& error) {
        std::cerr << "piece plan test failed: " << error.what() << "\n";
        return 1;
    }
    std::cout << "carta-zarr piece plan tests passed\n";
    return 0;
}
