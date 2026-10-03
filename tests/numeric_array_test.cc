/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// NumericArray is where a Zarr array's stored order stops being the caller's problem. The cases that
// matter are the ones a schema profile used to get wrong silently: a dimension it does not have, an
// index past the end of one it does, and a buffer of another length than the shape claims.

#include "zarr/numeric_array.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::internal::zarr::ArrayMetadata;
using carta::zarr::internal::zarr::NumericArray;

using carta::zarr::testing::Require;

// (time, frequency, polarization) = (2, 3, 2), values numbered in C order.
ArrayMetadata BeamShapedMetadata() {
    ArrayMetadata metadata;
    metadata.shape = {2, 3, 2};
    metadata.dimension_names = {"time", "frequency", "polarization"};
    return metadata;
}

std::vector<double> CountingValues(std::size_t count) {
    std::vector<double> values;
    values.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        values.push_back(static_cast<double>(index));
    }
    return values;
}

NumericArray Counted(ArrayMetadata metadata, std::size_t count) {
    auto array = NumericArray::Make("BEAM", std::move(metadata), CountingValues(count));
    Require(static_cast<bool>(array), "values as many as the shape declares were refused");
    return std::move(array.value());
}

void TestAddressesByName() {
    const auto view = Counted(BeamShapedMetadata(), 12);

    // C order: polarization varies fastest, then frequency, then time.
    Require(view.At({{"time", 0}, {"frequency", 0}, {"polarization", 0}}).value() == 0.0, "origin was misaddressed");
    Require(view.At({{"time", 0}, {"frequency", 0}, {"polarization", 1}}).value() == 1.0,
            "the fastest dimension was not the last one");
    Require(view.At({{"time", 0}, {"frequency", 1}, {"polarization", 0}}).value() == 2.0,
            "the middle dimension used the wrong stride");
    Require(view.At({{"time", 1}, {"frequency", 0}, {"polarization", 0}}).value() == 6.0,
            "the slowest dimension used the wrong stride");
    Require(view.At({{"time", 1}, {"frequency", 2}, {"polarization", 1}}).value() == 11.0,
            "the last element was misaddressed");
}

// Naming the dimensions in a different order must not change the element addressed. This is the
// property that lets a profile stop tracking a store's stored order.
void TestOrderOfNamesDoesNotMatter() {
    const auto view = Counted(BeamShapedMetadata(), 12);

    Require(view.At({{"time", 1}, {"frequency", 2}, {"polarization", 1}}).value() ==
                view.At({{"polarization", 1}, {"time", 1}, {"frequency", 2}}).value(),
            "the order the dimensions were named changed the element addressed");
}

void TestUnnamedDimensionsAreZero() {
    const auto view = Counted(BeamShapedMetadata(), 12);

    Require(
        view.At({{"frequency", 1}}).value() == view.At({{"time", 0}, {"frequency", 1}, {"polarization", 0}}).value(),
        "an unnamed dimension was not taken at index 0");
}

void TestUnknownDimensionIsAnError() {
    const auto view = Counted(BeamShapedMetadata(), 12);

    const auto missing = view.At({{"beam_params_label", 0}});
    Require(!missing && missing.error().code == ErrorCode::invalid_slice,
            "a dimension the array does not have was not reported as an invalid slice");
}

void TestIndexPastTheEndIsAnError() {
    const auto view = Counted(BeamShapedMetadata(), 12);

    const auto past = view.At({{"frequency", 3}});
    Require(!past && past.error().code == ErrorCode::invalid_slice,
            "an index past the end of a dimension was not reported as an invalid slice");
}

// The case the old hand-rolled indexing turned into a silent 0.0: metadata promising more elements
// than the store actually holds. It is refused when the values are bound to their metadata rather
// than when an element past the end is asked for, so an array that exists can always be read. A
// buffer longer than the shape is as wrong, and used to be addressed as if it were not.
void TestValuesAreAsManyAsTheShapeDeclares() {
    for (const std::size_t count : {std::size_t{4}, std::size_t{13}, std::size_t{0}}) {
        const auto refused = NumericArray::Make("BEAM", BeamShapedMetadata(), CountingValues(count));
        Require(!refused && refused.error().code == ErrorCode::invalid_metadata,
                std::to_string(count) + " values were bound to a shape declaring 12");
        Require(refused.error().node_path == "BEAM", "the error did not name the array");
    }

    ArrayMetadata empty;
    empty.shape = {2, 0};
    empty.dimension_names = {"time", "frequency"};
    Require(static_cast<bool>(NumericArray::Make("BEAM", empty, {})), "an array with an empty dimension was refused");

    // A shape whose product wraps to the buffer's length is not one that holds it.
    ArrayMetadata wrapping;
    wrapping.shape = {std::uint64_t{1} << 32U, std::uint64_t{1} << 32U, 3};
    Require(!NumericArray::Make("BEAM", wrapping, CountingValues(0)),
            "a shape declaring more values than can be counted was bound to an empty buffer");
}

void TestMissingDimensionNamesIsAnError() {
    ArrayMetadata metadata;
    metadata.shape = {2, 3};
    const auto view = Counted(metadata, 6);

    const auto unnamed = view.At({{"time", 0}});
    Require(!unnamed && unnamed.error().code == ErrorCode::invalid_metadata,
            "an array without dimension names was addressed by name anyway");
}

}  // namespace

int main() {
    try {
        TestAddressesByName();
        TestOrderOfNamesDoesNotMatter();
        TestUnnamedDimensionsAreZero();
        TestUnknownDimensionIsAnError();
        TestIndexPastTheEndIsAnError();
        TestValuesAreAsManyAsTheShapeDeclares();
        TestMissingDimensionNamesIsAnError();
        std::cout << "carta-zarr numeric array tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr numeric array tests failed: " << error.what() << '\n';
        return 1;
    }
}
