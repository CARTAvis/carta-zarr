/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "zarr/string_array.h"

#include "zarr/array_metadata.h"

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <new>
#include <string>
#include <vector>

#include "support/check.h"

// One allocation on the calling thread can be made to fail once, so that a decode can be run with
// each of its allocations failing in turn.
namespace {

thread_local bool counting = false;
thread_local std::size_t allocations = 0;
thread_local std::size_t failing = 0;

void* Allocate(std::size_t size) {
    if (counting) {
        ++allocations;
        if (allocations == failing) {
            failing = 0;
            throw std::bad_alloc();
        }
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) {
        return memory;
    }
    throw std::bad_alloc();
}

}  // namespace

void* operator new(std::size_t size) {
    return Allocate(size);
}
void* operator new[](std::size_t size) {
    return Allocate(size);
}
void operator delete(void* memory) noexcept {
    std::free(memory);
}
void operator delete[](void* memory) noexcept {
    std::free(memory);
}
void operator delete(void* memory, std::size_t) noexcept {
    std::free(memory);
}
void operator delete[](void* memory, std::size_t) noexcept {
    std::free(memory);
}

namespace {

using carta::zarr::ErrorCode;
namespace zarr_metadata = carta::zarr::internal::zarr;

const std::filesystem::path kFixtureDir{CARTA_ZARR_STRING_FIXTURE_DIR};

using carta::zarr::testing::Require;

zarr_metadata::ArrayMetadata FixtureMetadata(const std::string& name) {
    std::ifstream metadata_file(kFixtureDir / name / "zarr.json");
    Require(metadata_file.is_open(), "Missing fixture " + (kFixtureDir / name).string());
    const nlohmann::json metadata = nlohmann::json::parse(metadata_file);

    auto array_metadata = zarr_metadata::ParseArrayMetadata(metadata, name);
    Require(static_cast<bool>(array_metadata), "Fixture " + name + " has unreadable array metadata");
    return array_metadata.value();
}

carta::zarr::Result<std::vector<std::string>> ReadFixture(const std::string& name) {
    return zarr_metadata::ReadFixedLengthUtf32StringArray(kFixtureDir / name, FixtureMetadata(name), name);
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
    ExpectError("oversized", ErrorCode::decode_error);
    ExpectError("invalid_unicode", ErrorCode::decode_error);
}

// An allocation that fails is not corrupt data. The store remembers what a read of labels came to
// for as long as it lives, and asks again only after a read that could not be made; a decode_error
// is an answer about the array, so a shortage of memory reported as one kept good labels unreadable
// until the store was closed. With each of the decode's allocations failing in turn, what comes out
// is the std::bad_alloc itself -- which the store reports as a read it could not make -- or, where the
// failure was survived, the values.
void TestAnAllocationThatFailsIsNotCorruptData() {
    for (const auto* name : {"zstd_overhang", "gzip", "blosc", "multi_chunk"}) {
        const auto metadata = FixtureMetadata(name);
        const auto decode = [&]() {
            return zarr_metadata::ReadFixedLengthUtf32StringArray(kFixtureDir / name, metadata, name);
        };
        allocations = 0;
        counting = true;
        Require(static_cast<bool>(decode()), std::string("Fixture ") + name + " failed to decode");
        counting = false;
        const auto total = allocations;
        Require(total > 0, std::string("Decoding ") + name + " allocated nothing to fail");
        for (std::size_t which = 1; which <= total; ++which) {
            allocations = 0;
            failing = which;
            counting = true;
            std::string outcome;
            try {
                const auto result = decode();
                outcome = result ? "values" : std::string("a ") + zarr_metadata::ErrorCodeName(result.error().code);
            } catch (const std::bad_alloc&) {
                outcome = "a std::bad_alloc";
            }
            counting = false;
            failing = 0;
            Require(outcome == "values" || outcome == "a std::bad_alloc", std::string("Decoding ") + name + " with allocation " +
                                             std::to_string(which) + " of " + std::to_string(total) +
                                             " failing gave " + outcome);
        }
    }
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
        TestAnAllocationThatFailsIsNotCorruptData();
        std::cout << "carta-zarr string array tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr string array tests failed: " << error.what() << '\n';
        return 1;
    }
}
