/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "zarr/string_array.h"

#include "zarr/array_metadata.h"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::ErrorCode;
namespace zarr_metadata = carta::zarr::internal::zarr;

const std::filesystem::path kFixtureDir{CARTA_ZARR_STRING_FIXTURE_DIR};

using carta::zarr::testing::Require;

carta::zarr::Result<std::vector<std::string>> ReadFixture(const std::string& name) {
    const std::filesystem::path array_dir = kFixtureDir / name;
    std::ifstream metadata_file(array_dir / "zarr.json");
    Require(metadata_file.is_open(), "Missing fixture " + array_dir.string());
    const nlohmann::json metadata = nlohmann::json::parse(metadata_file);

    auto array_metadata = zarr_metadata::ParseArrayMetadata(metadata, name);
    Require(static_cast<bool>(array_metadata), "Fixture " + name + " has unreadable array metadata");
    return zarr_metadata::ReadFixedLengthUtf32StringArray(array_dir, array_metadata.value(), name);
}

void ExpectValues(const std::string& name, const std::vector<std::string>& expected) {
    auto result = ReadFixture(name);
    Require(static_cast<bool>(result),
            "Fixture " + name + " failed to decode: " + (result ? std::string{} : result.error().message));
    Require(result.value() == expected, "Fixture " + name + " decoded unexpected values");
}

void ExpectError(const std::string& name, ErrorCode expected) {
    auto result = ReadFixture(name);
    Require(!result, "Fixture " + name + " decoded successfully but should have failed");
    Require(result.error().code == expected, "Fixture " + name + " reported " +
                                                 zarr_metadata::ErrorCodeName(result.error().code) + " instead of " +
                                                 zarr_metadata::ErrorCodeName(expected));
    Require(result.error().node_path == name, "Fixture " + name + " did not report its node path");
}

// Codec chains and chunk key encodings that XRADIO can emit must all decode to the same values.
void TestSupportedLayouts() {
    const std::vector<std::string> expected{"A", "BC"};
    ExpectValues("zstd_overhang", expected);
    ExpectValues("gzip", expected);
    ExpectValues("blosc", expected);
    ExpectValues("crc_before_after_zstd", expected);
    ExpectValues("v2_key", expected);
    ExpectValues("dot_key", expected);
}

// The bytes codec endianness must be honoured, including for multi-byte code points.
void TestBigEndian() {
    ExpectValues("big_endian", std::vector<std::string>{"\xCE\xA9", "\xF0\x9F\x99\x82"});
}

// A chunk that was never written falls back to the empty fill value.
void TestMissingChunk() {
    ExpectValues("missing_chunk", std::vector<std::string>{"", ""});
}

// A label array split over several chunks reads as one, whichever way its chunk keys are spelled and
// whether or not its last chunk overhangs the array; a chunk never written takes the fill value
// while those around it are still read.
void TestMultipleChunks() {
    const std::vector<std::string> expected{"I", "Q", "U"};
    ExpectValues("multi_chunk", expected);
    ExpectValues("multi_chunk_overhang", expected);
    ExpectValues("multi_chunk_v2_key", expected);
    ExpectValues("multi_chunk_missing", std::vector<std::string>{"I", "", "U"});
}

// A chunk never written reads as the fill value the array declares, as zarr-python reads it, rather
// than as an empty string whatever was declared.
void TestDeclaredFillValue() {
    ExpectValues("missing_chunk_fill", std::vector<std::string>{"I", "I"});
    ExpectValues("multi_chunk_missing_fill", std::vector<std::string>{"I", "V", "U"});
}

// Corrupted chunks must be reported, never decoded past the end of the buffer.
void TestCorruptChunks() {
    ExpectError("crc_mismatch", ErrorCode::decode_error);
    ExpectError("truncated", ErrorCode::decode_error);
    ExpectError("invalid_unicode", ErrorCode::decode_error);
}

}  // namespace

int main() {
    try {
        Require(std::filesystem::is_directory(kFixtureDir), "Missing string fixture directory " + kFixtureDir.string());
        TestSupportedLayouts();
        TestBigEndian();
        TestMissingChunk();
        TestMultipleChunks();
        TestDeclaredFillValue();
        TestCorruptChunks();
        std::cout << "carta-zarr string array tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr string array tests failed: " << error.what() << '\n';
        return 1;
    }
}
