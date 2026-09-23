/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a pixel selection is, stated directly.
//
// A selection was reached by no test until this: every claim about it was made by reading pixels
// and checking the answer, which is how a selection whose order field had been overwritten to mean
// something else -- consistent with itself, and transposed -- could have gone unnoticed by any test
// that only compared totals.
//
// This target links nothing. The question is a descriptor and a request in, a selection out.

#include "zarr/pixel_selection.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadRequest;
using carta::zarr::internal::zarr::BuildSelection;
using carta::zarr::internal::zarr::DestinationOrder;
using carta::zarr::internal::zarr::PixelSelection;

using carta::zarr::testing::Require;

// An image in the layout XRADIO writes: logical l, m, frequency, polarization, time, stored as time,
// frequency, polarization, l, m -- so m is last and the array is contiguous along it.
ImageDescriptor MakeImage() {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
        std::size_t stored;
    } axes[] = {{"l", AxisRole::spatial_x, 4, 3},
                {"m", AxisRole::spatial_y, 5, 4},
                {"frequency", AxisRole::spectral, 3, 1},
                {"polarization", AxisRole::polarization, 2, 2},
                {"time", AxisRole::time, 1, 0}};
    for (const auto& axis : axes) {
        carta::zarr::AxisDescriptor described;
        described.name = axis.name;
        described.role = axis.role;
        described.length = axis.length;
        described.storage_index = axis.stored;
        descriptor.axes.push_back(described);
    }
    return descriptor;
}

ReadRequest Whole() {
    ReadRequest request;
    request.axes = {Range{0, 4, 1}, Range{0, 5, 1}, Range{0, 3, 1}, Range{0, 2, 1}, Range{0, 1, 1}};
    return request;
}

PixelSelection Built(DestinationOrder order, const ReadRequest& request = Whole()) {
    auto built = BuildSelection(MakeImage(), request, order);
    Require(static_cast<bool>(built), "BuildSelection failed: " + (built ? std::string{} : built.error().message));
    return built.value();
}

// Logical order lays the destination out the way the library reports an image: its axis i is the
// stored dimension that logical axis i is.
void TestALogicalDestinationFollowsTheImage() {
    const auto selection = Built(DestinationOrder::logical);
    Require(selection.destination_to_stored == std::vector<std::size_t>{3, 4, 1, 2, 0},
            "a logical destination was not laid out in the image's logical order");
    // l fastest, then m (4 of l), frequency (4 x 5), polarization (4 x 5 x 3), time (all 120).
    Require(selection.DestinationStrides() == std::vector<std::uint64_t>{120, 20, 60, 1, 4},
            "a logical destination's strides did not put l fastest");
}

// Stored order lays it out as the array was written, the last stored dimension fastest, so that a
// plane arrives untransposed. The strides are the ones a pass used to work out for itself: the last
// stored dimension steps by one and each earlier one by the product of those after it.
void TestAStoredDestinationFollowsTheArray() {
    const auto selection = Built(DestinationOrder::stored);
    Require(selection.destination_to_stored == std::vector<std::size_t>{4, 3, 2, 1, 0},
            "a stored destination was not laid out in reversed stored order");

    std::vector<std::uint64_t> by_hand(selection.count.size(), 1);
    std::uint64_t running = 1;
    for (std::size_t stored = selection.count.size(); stored-- > 0;) {
        by_hand.at(stored) = running;
        running *= selection.count.at(stored);
    }
    Require(selection.DestinationStrides() == by_hand,
            "a stored destination's strides were not the ones a pass computed by hand");
    Require(selection.DestinationStrides().at(4) == 1, "the dimension the array is contiguous along was not fastest");
}

// The order says where pixels land, never which ones are read.
void TestTheOrderDoesNotChangeWhatIsRead() {
    ReadRequest request = Whole();
    request.axes.at(0) = Range{1, 2, 2};
    request.axes.at(2) = Range{0, 2, 2};
    const auto logical = Built(DestinationOrder::logical, request);
    const auto stored = Built(DestinationOrder::stored, request);
    Require(logical.start == stored.start && logical.count == stored.count && logical.stride == stored.stride,
            "the destination's order changed which pixels the selection reads");
    // In stored order, whatever the request.
    Require(logical.start == std::vector<std::uint64_t>{0, 0, 0, 1, 0} &&
                logical.count == std::vector<std::uint64_t>{1, 2, 2, 2, 5} &&
                logical.stride == std::vector<std::uint64_t>{1, 2, 1, 2, 1},
            "the request was not translated into the array's stored order");
}

}  // namespace

int main() {
    try {
        TestALogicalDestinationFollowsTheImage();
        TestAStoredDestinationFollowsTheArray();
        TestTheOrderDoesNotChangeWhatIsRead();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "pixel selection test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
