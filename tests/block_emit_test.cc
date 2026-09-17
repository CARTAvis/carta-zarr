/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a reduction decides about its blocks, before it reads anything.
//
// The envelope around the pass was written out twice, in full, and what the two copies had in
// common is all here: where a block is cut, which numbers the walk is handed, how much of a
// part-filled block is reported as done, and what a sink saying no means. None of it needs a pixel,
// so the walk below reads nothing at all -- it only says how far it got.

#include "reduce/block_emit.h"
#include "reduce/plane_selection.h"

#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ChunkGeometry;
using carta::zarr::ErrorCode;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadOptions;
using carta::zarr::Result;
using carta::zarr::internal::BlockEmitter;
using carta::zarr::internal::CheckedPlanes;
using carta::zarr::internal::EmitBlock;
using carta::zarr::internal::MapAxes;
using carta::zarr::internal::PassPlan;
using carta::zarr::internal::PlanPass;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageDescriptor MakeImage(std::uint64_t channels) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, 64},
             {"m", AxisRole::spatial_y, 64},
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

PassPlan Plan(const ImageDescriptor& descriptor, std::uint64_t chunk_z, const Range& spectral) {
    ChunkGeometry geometry;
    geometry.fastest_spatial_axis = AxisRole::spatial_y;
    geometry.chunk_shape = {64, 64, chunk_z, 1, 1};
    const auto map = MapAxes(descriptor);
    Require(static_cast<bool>(map), "MapAxes failed on a well-formed image");
    const auto planes = CheckedPlanes::Of(descriptor, map.value(), {spectral, 0, 0});
    Require(static_cast<bool>(planes), "the spectral range does not fit this image");
    return PlanPass(descriptor, geometry, map.value(), planes.value(), 1, ReadOptions{});
}

// One hand-over, as the sink saw it.
struct Handed {
    std::uint64_t first_channel = 0;
    std::uint64_t length = 0;
    bool complete = false;
    double completeness = 0.0;
};

// The two numbers a walk is given, as the walk saw them.
struct Walked {
    std::uint64_t begin = 0;
    std::uint64_t end = 0;
};

// A block is never cut inside a spectral chunk, however small a block the caller asks for: a decode
// serving two blocks would split one read's results across them. So a hint of one channel against a
// chunk four deep still emits blocks of four.
void TestBlocksAreCutOnChunkBoundaries() {
    const auto image = MakeImage(32);
    const auto plan = Plan(image, 4, Range{0, 32, 1});
    const BlockEmitter emitter(plan, 10, sizeof(double), 1);

    std::vector<Handed> handed;
    std::vector<Walked> walked;
    const auto outcome = emitter.Over(
        [](std::uint64_t) {},
        [&](EmitBlock& block, const auto&) -> Result<void> {
            walked.push_back({block.begin, block.end});
            return {};
        },
        [&](std::uint64_t first, std::uint64_t length, bool complete, double completeness) {
            handed.push_back({first, length, complete, completeness});
            return true;
        });
    Require(static_cast<bool>(outcome), "an uncancelled walk should succeed");

    Require(handed.size() == 8, "32 channels in blocks of a 4-deep chunk is eight blocks, not " +
                                    std::to_string(handed.size()));
    Require(walked.size() == handed.size(), "every block should have been walked exactly once");

    // The one confusion this module exists to own: what the walk is given is absolute into the
    // selection, and so is what the sink is told -- a later block does not start at zero.
    std::uint64_t expected = 0;
    for (std::size_t i = 0; i < handed.size(); ++i) {
        Require(handed.at(i).first_channel == expected,
                "block " + std::to_string(i) + " should start at " + std::to_string(expected));
        Require(handed.at(i).length == 4, "block " + std::to_string(i) + " should hold four channels");
        Require(walked.at(i).begin == expected && walked.at(i).end == expected + 4,
                "the walk of block " + std::to_string(i) + " should be given the same absolute range");
        Require(handed.at(i).complete, "a block walked without reporting is handed over once, finished");
        Require(handed.at(i).completeness == 1.0, "and a finished block is complete");
        expected += 4;
    }
    Require(expected == 32, "the blocks should tile the selection");
}

// A block whose walk takes more than one read is handed over as it fills, and what it says about
// itself is the fraction of its own chunks that are in it -- counted in the layer the emit budget
// was spent against, which is the same number for both.
void TestAPartFilledBlockReportsItsOwnChunks() {
    const auto image = MakeImage(8);
    const auto plan = Plan(image, 8, Range{0, 8, 1});
    // One layer of ten chunks, and one block: eight channels of an eight-deep chunk is one chunk
    // along the spectrum, so the block's chunks are the layer's ten.
    const BlockEmitter emitter(plan, 10, sizeof(double), 0);

    std::vector<Handed> handed;
    const auto outcome = emitter.Over(
        [](std::uint64_t) {},
        [&](EmitBlock& block, const auto& report) -> Result<void> {
            block.chunks_done = 3;
            if (auto reported = report(block.chunks_done); !reported) {
                return reported.error();
            }
            block.chunks_done = 7;
            if (auto reported = report(block.chunks_done); !reported) {
                return reported.error();
            }
            block.chunks_done = 10;
            return {};
        },
        [&](std::uint64_t first, std::uint64_t length, bool complete, double completeness) {
            handed.push_back({first, length, complete, completeness});
            return true;
        });
    Require(static_cast<bool>(outcome), "an uncancelled walk should succeed");

    Require(handed.size() == 3, "two reports and one finish is three hand-overs, not " +
                                    std::to_string(handed.size()));
    Require(!handed.at(0).complete && handed.at(0).completeness == 3.0 / 10.0, "three chunks of ten");
    Require(!handed.at(1).complete && handed.at(1).completeness == 7.0 / 10.0, "seven chunks of ten");
    Require(handed.at(2).complete && handed.at(2).completeness == 1.0,
            "a finished block says one exactly, not ten tenths");
    for (const auto& one : handed) {
        Require(one.first_channel == 0 && one.length == 8, "every hand-over describes the same block");
    }
}

// A sink saying no stops the reduction there, whether it says it to a part-filled block or to a
// finished one, and the reduction reports cancelled rather than a short answer that looks complete.
void TestASinkThatSaysNoCancels() {
    const auto image = MakeImage(32);

    // Refusing the first finished block: the second is never cut.
    {
        const auto plan = Plan(image, 4, Range{0, 32, 1});
        const BlockEmitter emitter(plan, 10, sizeof(double), 1);
        std::uint64_t blocks_walked = 0;
        std::uint64_t hand_overs = 0;
        const auto outcome = emitter.Over(
            [](std::uint64_t) {},
            [&](EmitBlock&, const auto&) -> Result<void> {
                ++blocks_walked;
                return {};
            },
            [&](std::uint64_t, std::uint64_t, bool, double) {
                ++hand_overs;
                return false;
            });
        Require(!outcome, "a sink returning false should fail the reduction");
        Require(outcome.error().code == ErrorCode::cancelled, "and it should report cancelled");
        Require(blocks_walked == 1 && hand_overs == 1, "and it should stop at the first block");
    }

    // Refusing a part-filled one: the error comes back through the walk's own report, so the rest
    // of that block is not read either.
    {
        const auto plan = Plan(image, 4, Range{0, 32, 1});
        const BlockEmitter emitter(plan, 10, sizeof(double), 1);
        bool read_on = false;
        const auto outcome = emitter.Over(
            [](std::uint64_t) {},
            [&](EmitBlock& block, const auto& report) -> Result<void> {
                block.chunks_done = 1;
                if (auto reported = report(block.chunks_done); !reported) {
                    return reported.error();
                }
                read_on = true;
                return {};
            },
            [&](std::uint64_t, std::uint64_t, bool, double) { return false; });
        Require(!outcome && outcome.error().code == ErrorCode::cancelled,
                "refusing a part-filled block should cancel too");
        Require(!read_on, "and the rest of that block should not be read");
    }
}

// The hint is the caller's and the chunk alignment is the plan's, but the emit budget is the
// library's and it wins when a channel is expensive enough. Asked for the whole selection in one
// block, a request whose channels are a megabyte each gets fewer.
void TestTheBudgetLowersTheCallersHint() {
    const auto image = MakeImage(64);
    const auto plan = Plan(image, 1, Range{0, 64, 1});

    const BlockEmitter modest(plan, 1, sizeof(double), 64);
    Require(modest.emit_channels() == 64, "a cheap channel takes the hint as given");

    const BlockEmitter expensive(plan, 1, 64U << 20U, 64);
    Require(expensive.emit_channels() < 64,
            "a channel costing 64 MiB should not be emitted 64 at a time");
    Require(expensive.emit_channels() >= 1, "and a block always holds at least one channel");
}

}  // namespace

int main() {
    try {
        TestBlocksAreCutOnChunkBoundaries();
        TestAPartFilledBlockReportsItsOwnChunks();
        TestASinkThatSaysNoCancels();
        TestTheBudgetLowersTheCallersHint();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "block emit test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
