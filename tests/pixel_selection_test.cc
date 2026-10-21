/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// What a pixel selection is, stated directly.
//
// A selection is asserted here directly rather than through the pixels it reads, so that one
// consistent with itself and transposed cannot pass a test that only compares totals.
//
// This target links nothing. The question is a descriptor and a request in, a selection out.

#include "zarr/pixel_selection.h"

#include "support/check.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <tuple>
#include <vector>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadRequest;
using carta::zarr::internal::zarr::BuildSelection;
using carta::zarr::internal::zarr::BuildSelectionWithin;
using carta::zarr::internal::zarr::DestinationOrder;
using carta::zarr::internal::zarr::PixelSelection;

using carta::zarr::testing::Require;

// An image in the layout XRADIO writes: logical l, m, frequency, polarization, time, stored as time,
// frequency, polarization, l, m -- so m is last and the array is contiguous along it.
ImageDescriptor MakeImage(std::uint64_t l = 4, std::uint64_t m = 5, std::uint64_t frequency = 3) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
        std::size_t stored;
    } axes[] = {{"l", AxisRole::spatial_x, l, 3},
                {"m", AxisRole::spatial_y, m, 4},
                {"frequency", AxisRole::spectral, frequency, 1},
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
    Require(selection.destination_to_stored() == std::vector<std::size_t>{3, 4, 1, 2, 0},
            "a logical destination was not laid out in the image's logical order");
    // l fastest, then m (4 of l), frequency (4 x 5), polarization (4 x 5 x 3), time (all 120).
    Require(selection.DestinationStrides() == std::vector<std::uint64_t>{120, 20, 60, 1, 4},
            "a logical destination's strides did not put l fastest");
}

// Stored order lays it out as the array was written, the last stored dimension fastest, so that a
// plane arrives untransposed: the last stored dimension steps by one and each earlier one by the
// product of those after it.
void TestAStoredDestinationFollowsTheArray() {
    const auto selection = Built(DestinationOrder::stored);
    Require(selection.destination_to_stored() == std::vector<std::size_t>{4, 3, 2, 1, 0},
            "a stored destination was not laid out in reversed stored order");

    std::vector<std::uint64_t> by_hand(selection.count().size(), 1);
    std::uint64_t running = 1;
    for (std::size_t stored = selection.count().size(); stored-- > 0;) {
        by_hand.at(stored) = running;
        running *= selection.count().at(stored);
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
    Require(
        logical.start() == stored.start() && logical.count() == stored.count() && logical.stride() == stored.stride(),
        "the destination's order changed which pixels the selection reads");
    // In stored order, whatever the request.
    Require(logical.start() == std::vector<std::uint64_t>{0, 0, 0, 1, 0} &&
                logical.count() == std::vector<std::uint64_t>{1, 2, 2, 2, 5} &&
                logical.stride() == std::vector<std::uint64_t>{1, 2, 1, 2, 1},
            "the request was not translated into the array's stored order");
}

// A part of a larger read lands where it is in that read's dense logical-order destination: each step
// moves as far as it does there, the first element is where the part begins, and the destination has
// to reach as far as the part's last element. A dense selection starts at zero and reaches its count.
void TestAPartLandsWhereItIsInTheWhole() {
    const auto dense = Built(DestinationOrder::logical);
    Require(dense.destination_offset() == 0 && dense.destination_extent() == dense.elements(),
            "a dense selection did not start at zero and reach its own count");

    ReadRequest part = Whole();
    part.axes.at(0) = Range{1, 2, 1};
    part.axes.at(1) = Range{2, 3, 1};
    part.axes.at(2) = Range{1, 1, 1};
    const auto within = BuildSelectionWithin(MakeImage(), part, Whole());
    Require(static_cast<bool>(within), "a part of the whole read was refused");
    // The whole's strides, l 1, m 4, frequency 20, polarization 60, time 120, by stored dimension.
    Require(within.value().DestinationStrides() == std::vector<std::uint64_t>{120, 20, 60, 1, 4},
            "a part did not step as the whole read's destination does");
    Require(within.value().destination_offset() == (1 * 1) + (2 * 4) + (1 * 20),
            "a part did not start where it is in the whole");
    Require(within.value().destination_extent() == 29 + (1 * 1) + (2 * 4) + (1 * 60) + 1,
            "a part did not reach as far as its last element");
    Require(within.value().elements() == 2 * 3 * 2, "a part did not count only its own elements");

    // A strided whole: the part's elements must be the whole's.
    ReadRequest strided = Whole();
    strided.axes.at(0) = Range{0, 2, 2};
    ReadRequest on = strided;
    on.axes.at(0) = Range{2, 1, 1};
    Require(static_cast<bool>(BuildSelectionWithin(MakeImage(), on, strided)),
            "an element of a strided whole was refused");
    ReadRequest between = strided;
    between.axes.at(0) = Range{1, 1, 1};
    ReadRequest finer = strided;
    finer.axes.at(0) = Range{0, 2, 1};
    ReadRequest beyond = Whole();
    beyond.axes.at(1) = Range{3, 2, 1};
    ReadRequest short_whole = Whole();
    short_whole.axes.at(1) = Range{0, 4, 1};
    for (const auto& [outside, of, what] : {std::tuple{between, strided, "between two of the whole's elements"},
                                            std::tuple{finer, strided, "stepping finer than the whole"},
                                            std::tuple{beyond, short_whole, "past the end of the whole"}}) {
        const auto refused = BuildSelectionWithin(MakeImage(), outside, of);
        Require(!refused && refused.error().code == carta::zarr::ErrorCode::invalid_argument,
                std::string("a part ") + what + " was not refused");
    }
}

// Counted once, when the selection is built, and never zero: a selection that would produce nothing
// is refused rather than made.
void TestASelectionKnowsHowManyElementsItProduces() {
    Require(Built(DestinationOrder::logical).elements() == 4 * 5 * 3 * 2 * 1,
            "the whole image did not count as 120 elements");
    ReadRequest request = Whole();
    request.axes.at(0) = Range{1, 2, 2};
    Require(Built(DestinationOrder::stored, request).elements() == 2 * 5 * 3 * 2,
            "a strided selection did not count what it selects");
}

// Every way a request can fail to be a selection, each refused where the selection is made -- which
// is what lets a reader of pixels take one on trust.
void TestAMalformedRequestNeverBecomesASelection() {
    const auto refused = [](const ImageDescriptor& image, const ReadRequest& request, const std::string& what) {
        const auto built = BuildSelection(image, request, DestinationOrder::logical);
        Require(!built && built.error().code == carta::zarr::ErrorCode::invalid_slice, what + " was accepted");
    };

    ReadRequest short_request = Whole();
    short_request.axes.pop_back();
    refused(MakeImage(), short_request, "a request with fewer axes than the image");

    ReadRequest zero_stride = Whole();
    zero_stride.axes.at(0).stride = 0;
    refused(MakeImage(), zero_stride, "a zero stride");

    ReadRequest zero_count = Whole();
    zero_count.axes.at(1).count = 0;
    refused(MakeImage(), zero_count, "an axis selecting nothing");

    ReadRequest late_start = Whole();
    late_start.axes.at(0) = Range{4, 1, 1};
    refused(MakeImage(), late_start, "a start past the end of its axis");

    ReadRequest long_run = Whole();
    long_run.axes.at(0) = Range{1, 4, 1};
    refused(MakeImage(), long_run, "a range running past the end of its axis");

    // A count of 2^32 + 1 and a stride of 2^32: the last index multiplied out is 2^64, which wraps to
    // zero and looks as though it fits in any axis. Compared by division instead, it does not.
    ReadRequest wrapping = Whole();
    wrapping.axes.at(0) = Range{0, (std::uint64_t{1} << 32) + 1, std::uint64_t{1} << 32};
    refused(MakeImage(std::uint64_t{1} << 40), wrapping, "a span that only fits when it overflows");

    // Every axis in range, and more elements than a count can hold. A count of zero here would be
    // reported by a read as a destination too small for the request.
    const auto huge = MakeImage(std::uint64_t{1} << 32, std::uint64_t{1} << 32, std::uint64_t{1} << 32);
    ReadRequest everything = Whole();
    everything.axes.at(0) = Range{0, std::uint64_t{1} << 32, 1};
    everything.axes.at(1) = Range{0, std::uint64_t{1} << 32, 1};
    everything.axes.at(2) = Range{0, std::uint64_t{1} << 32, 1};
    const auto counted = BuildSelection(huge, everything, DestinationOrder::logical);
    Require(!counted && counted.error().code == carta::zarr::ErrorCode::invalid_slice &&
                counted.error().message.find("more elements than can be counted") != std::string::npos,
            "a selection too large to count was not refused as one");
}

}  // namespace

int main() {
    try {
        TestALogicalDestinationFollowsTheImage();
        TestAStoredDestinationFollowsTheArray();
        TestTheOrderDoesNotChangeWhatIsRead();
        TestAPartLandsWhereItIsInTheWhole();
        TestASelectionKnowsHowManyElementsItProduces();
        TestAMalformedRequestNeverBecomesASelection();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "pixel selection test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
