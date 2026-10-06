/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// What a variable has to be before it can mask an image's pixels.
//
// Every case here is two pieces of metadata and an answer, so this target links nothing at all:
// RequireUsableFlag and IsFlag are inline, and neither reaches a Store. Choosing a flag when the
// image declares none does reach one -- that half lives in tests/schema_profile_test.cc, against
// the in-memory transport.
//
// These four refusals were checked by opening four stores on disk. They are the same four.

#include "schema/xradio/flag.h"

#include "support/check.h"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::internal::xradio::IsFlag;
using carta::zarr::internal::xradio::RequireUsableFlag;
using ArrayMetadata = carta::zarr::internal::zarr::ArrayMetadata;

using carta::zarr::testing::Require;

const std::vector<std::string> kSkyDimensions{"time", "frequency", "polarization", "l", "m"};
const std::vector<std::uint64_t> kSkyShape{1, 3, 2, 4, 5};

ArrayMetadata Image() {
    ArrayMetadata image;
    image.shape = kSkyShape;
    image.dimension_names = kSkyDimensions;
    image.data_type = "float32";
    image.attributes = nlohmann::json::object();
    return image;
}

ArrayMetadata Flag(std::string data_type = "bool", std::vector<std::uint64_t> shape = kSkyShape,
                   std::vector<std::string> dimensions = kSkyDimensions, bool marked = true) {
    ArrayMetadata flag;
    flag.shape = std::move(shape);
    flag.dimension_names = std::move(dimensions);
    flag.data_type = std::move(data_type);
    flag.attributes = marked ? nlohmann::json{{"type", "flag"}} : nlohmann::json::object();
    return flag;
}

void TestAWellFormedFlagIsAccepted() {
    const auto usable = RequireUsableFlag(Flag(), Image(), "MASK_0");
    Require(static_cast<bool>(usable), "a boolean flag over the image's own dimensions was refused");
}

// A variable is a flag because of its type attribute, never because of its name. MASK_0 is a flag
// and MASK_DECONVOLVE is an ordinary image, and the two are told apart here and nowhere else.
void TestClassificationIsNotByName() {
    Require(IsFlag(Flag()), "a variable marked type: flag was not recognised as one");
    Require(!IsFlag(Flag("bool", kSkyShape, kSkyDimensions, false)),
            "an unmarked boolean variable was recognised as a flag");
    Require(!IsFlag(Image()), "an image was recognised as a flag");

    const auto unmarked = RequireUsableFlag(Flag("bool", kSkyShape, kSkyDimensions, false), Image(), "MASK_0");
    Require(!unmarked && unmarked.error().code == ErrorCode::invalid_metadata,
            "a boolean variable that is not a flag was accepted as a pixel mask");
}

void TestANumericVariableIsNotAMask() {
    const auto numeric = RequireUsableFlag(Flag("float32"), Image(), "MASK_0");
    Require(!numeric && numeric.error().code == ErrorCode::unsupported_data_type,
            "a numeric variable was accepted as a pixel mask");
}

// The read applies the mask element by element against the selection it read pixels with, so a flag
// that is the image's own shape in another order is not a mask of it. Both halves matter: a
// reshaped flag and a transposed one fail for the same reason, and both are checked without a
// store.
void TestTheFlagMustMatchTheImagesOwnDimensions() {
    const auto reshaped = RequireUsableFlag(Flag("bool", {1, 3, 2, 4, 4}), Image(), "MASK_0");
    Require(!reshaped && reshaped.error().code == ErrorCode::invalid_metadata,
            "a flag whose shape differs from the image was accepted as a pixel mask");

    const auto transposed = RequireUsableFlag(
        Flag("bool", {1, 3, 2, 5, 4}, {"time", "frequency", "polarization", "m", "l"}), Image(), "MASK_0");
    Require(!transposed && transposed.error().code == ErrorCode::invalid_metadata,
            "a flag whose dimension order differs from the image was accepted as a pixel mask");
}

// The refusal names the variable, because a consumer shows it to whoever picked the file.
// xarray stores a bool variable as int8 and records `dtype: "bool"` in its attributes, so that is
// what XRADIO's own writer leaves on disk for every flag it converts. An int8 without that record
// is a number, not a flag, and so is any other integer carrying it: xarray writes only int8.
void TestXarraysEncodingOfABoolIsAFlag() {
    auto encoded = Flag("int8");
    encoded.attributes["dtype"] = "bool";
    Require(static_cast<bool>(RequireUsableFlag(encoded, Image(), "FLAG_SKY")),
            "an int8 flag recorded as bool, which is how XRADIO writes one, was refused");

    const auto bare = RequireUsableFlag(Flag("int8"), Image(), "FLAG_SKY");
    Require(!bare && bare.error().code == ErrorCode::unsupported_data_type,
            "an int8 flag with no record that it holds booleans was accepted");

    auto widened = Flag("uint16");
    widened.attributes["dtype"] = "bool";
    const auto other = RequireUsableFlag(widened, Image(), "FLAG_SKY");
    Require(!other && other.error().code == ErrorCode::unsupported_data_type,
            "an integer other than int8 was taken as xarray's encoding of a bool");
}

void TestTheRefusalNamesTheVariable() {
    const auto numeric = RequireUsableFlag(Flag("float32"), Image(), "MASK_7");
    Require(!numeric && numeric.error().node_path == "MASK_7", "the refusal did not name the flag variable");
}

}  // namespace

int main() {
    try {
        TestAWellFormedFlagIsAccepted();
        TestClassificationIsNotByName();
        TestANumericVariableIsNotAMask();
        TestTheFlagMustMatchTheImagesOwnDimensions();
        TestXarraysEncodingOfABoolIsAFlag();
        TestTheRefusalNamesTheVariable();
        std::cout << "carta-zarr flag tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr flag tests failed: " << error.what() << '\n';
        return 1;
    }
}
