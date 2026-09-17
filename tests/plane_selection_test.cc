/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What a reduction is allowed to ask for, asked of nothing but a descriptor.
//
// This target compiles no implementation source at all. That is the argument for the module as much
// as the assertions are: every one of these used to be reachable only through ComputeHistogram or
// ReduceSpectral over a fixture, so the same refusal was asserted twice, in two suites that each
// needed a store on disk to reach it -- and the spectral third of it was asserted somewhere else
// again.

#include "reduce/plane_selection.h"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::AxisRole;
using carta::zarr::ErrorCode;
using carta::zarr::ImageDescriptor;
using carta::zarr::PlaneSelection;
using carta::zarr::Range;
using carta::zarr::internal::AxisMap;
using carta::zarr::internal::CheckedPlanes;
using carta::zarr::internal::MapAxes;

using carta::zarr::testing::Require;

// An image in the logical order XRADIO reports. A length of zero leaves the axis out entirely,
// which is the case that separates "the image has no such axis" from "index 0 of an axis of one".
ImageDescriptor MakeImage(std::uint64_t channels, std::uint64_t polarizations, std::uint64_t times) {
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
             {"polarization", AxisRole::polarization, polarizations},
             {"time", AxisRole::time, times}};
    for (const auto& wanted : axes) {
        if (wanted.length == 0) {
            continue;
        }
        carta::zarr::AxisDescriptor axis;
        axis.name = wanted.name;
        axis.role = wanted.role;
        axis.length = wanted.length;
        axis.storage_index = descriptor.axes.size();
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

AxisMap MapOf(const ImageDescriptor& descriptor) {
    auto map = MapAxes(descriptor);
    Require(static_cast<bool>(map), "MapAxes failed on a well-formed image");
    return map.value();
}

// Asks the question the way a reduction asks it, so that each test below is one line.
carta::zarr::Result<CheckedPlanes> Check(const ImageDescriptor& descriptor, const PlaneSelection& planes) {
    return CheckedPlanes::Of(descriptor, MapOf(descriptor), planes);
}

void Accepts(const ImageDescriptor& descriptor, const PlaneSelection& planes, const std::string& what) {
    const auto checked = Check(descriptor, planes);
    Require(static_cast<bool>(checked), what + " should have been accepted");
}

void Refuses(const ImageDescriptor& descriptor, const PlaneSelection& planes, const std::string& what) {
    const auto checked = Check(descriptor, planes);
    Require(!checked, what + " should have been refused");
    Require(checked.error().code == ErrorCode::invalid_argument,
            what + " should be refused as invalid_argument");
}

// A selection that fits comes back saying what it was asked for, and saying it once: a caller that
// took three loose numbers had to carry all three to everywhere that needed any of them.
void TestASelectionThatFitsIsReportedBack() {
    const auto image = MakeImage(32, 4, 2);
    const PlaneSelection planes{Range{2, 5, 3}, 3, 1};

    const auto checked = Check(image, planes);
    Require(static_cast<bool>(checked), "a selection inside the image should have been accepted");
    Require(checked.value().spectral().start == 2 && checked.value().spectral().count == 5 &&
                checked.value().spectral().stride == 3,
            "the spectral range should come back unchanged");
    Require(checked.value().selection().polarization == 3, "the polarization should come back unchanged");
    Require(checked.value().selection().time == 1, "the time should come back unchanged");
    Require(checked.value().count() == 5, "count() is how many planes were selected");
}

// The last selected channel is start + (count - 1) * stride, and it has to exist. The boundary is
// worth both sides: one past it is the off-by-one, and exactly on it is the range a caller asking
// for every nth channel of a cube actually writes.
void TestTheSpectralRangeHasToFit() {
    const auto image = MakeImage(32, 4, 2);

    Accepts(image, {Range{0, 32, 1}, 0, 0}, "every channel");
    Accepts(image, {Range{31, 1, 1}, 0, 0}, "the last channel alone");
    Accepts(image, {Range{1, 11, 3}, 0, 0}, "a stride landing exactly on the last channel");

    Refuses(image, {Range{0, 33, 1}, 0, 0}, "one channel more than the image has");
    Refuses(image, {Range{32, 1, 1}, 0, 0}, "a start past the last channel");
    Refuses(image, {Range{1, 12, 3}, 0, 0}, "a stride one step past the last channel");
    Refuses(image, {Range{0, 0, 1}, 0, 0}, "a range of no channels");
    Refuses(image, {Range{0, 4, 0}, 0, 0}, "a stride of zero");
}

// cad6c8a, in the one place it can now be fixed. Multiplying the span out -- start + (count - 1) *
// stride -- wraps for these, and a wrapped span lands back inside the image and passes.
void TestASpanThatWouldWrapIsRefused() {
    const auto image = MakeImage(32, 4, 2);
    constexpr std::uint64_t kHuge = std::uint64_t{1} << 32U;

    Refuses(image, {Range{0, kHuge + 1, kHuge}, 0, 0}, "a span that wraps to zero");
    Refuses(image, {Range{0, 3, std::uint64_t{0} - 1}, 0, 0}, "a span that wraps past the maximum");
}

// Zero on an axis the image does not have is the image itself; any other index there names a plane
// that does not exist, which is a different mistake from asking for one that is out of range.
void TestPolarizationAndTimeAreCheckedAgainstTheirAxes() {
    const auto full = MakeImage(32, 4, 2);
    Accepts(full, {Range{0, 32, 1}, 3, 1}, "the last polarization and the last time");
    Refuses(full, {Range{0, 32, 1}, 4, 0}, "a polarization past the axis");
    Refuses(full, {Range{0, 32, 1}, 0, 2}, "a time past the axis");

    const auto no_polarization = MakeImage(32, 0, 2);
    Accepts(no_polarization, {Range{0, 32, 1}, 0, 0}, "polarization 0 of an image with no polarization axis");
    Refuses(no_polarization, {Range{0, 32, 1}, 1, 0}, "polarization 1 of an image with no polarization axis");

    const auto no_time = MakeImage(32, 4, 0);
    Accepts(no_time, {Range{0, 32, 1}, 2, 0}, "time 0 of an image with no time axis");
    Refuses(no_time, {Range{0, 32, 1}, 2, 1}, "time 1 of an image with no time axis");
}

}  // namespace

int main() {
    try {
        TestASelectionThatFitsIsReportedBack();
        TestTheSpectralRangeHasToFit();
        TestASpanThatWouldWrapIsRefused();
        TestPolarizationAndTimeAreCheckedAgainstTheirAxes();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "plane selection test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
