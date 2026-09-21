/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What an array's zarr.json says about how it is stored.
//
// These were Store methods reaching into the raw node document, so the only way to ask them was to
// write a store to a temporary directory, open a Context, open a Dataset and open an Image -- six
// files and four objects to check four facts about one codec chain. Two of the refusals below could
// not be reached from any test at all, because a store that reaches them is one no fixture would be
// written to contain.
//
// This target links nothing. A document goes in and a layout comes out.

#include "zarr/array_metadata.h"

#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::StorageLayout;
using carta::zarr::internal::zarr::ArrayMetadata;
using carta::zarr::internal::zarr::ParseArrayMetadata;
using carta::zarr::internal::zarr::ParseStorageLayout;

using carta::zarr::testing::Require;

// One array document, with whatever codec chain the case is about.
std::string Document(const std::string& codecs, const std::string& chunk_shape = "[4,4]") {
    return R"({"shape":[8,8],"data_type":"float32","chunk_grid":{"name":"regular","configuration":)"
           R"({"chunk_shape":)" + chunk_shape + R"(}},"attributes":{},"dimension_names":["l","m"],)" +
           (codecs.empty() ? std::string{} : "\"codecs\":" + codecs + ",") +
           R"("zarr_format":3,"node_type":"array"})";
}

ArrayMetadata Parsed(const std::string& document) {
    auto metadata = ParseArrayMetadata(nlohmann::json::parse(document), "SKY");
    Require(static_cast<bool>(metadata),
            "the document did not parse: " + (metadata ? std::string{} : metadata.error().message));
    return std::move(metadata.value());
}

StorageLayout LayoutOf(const std::string& codecs, const std::string& chunk_shape = "[4,4]") {
    auto layout = ParseStorageLayout(Parsed(Document(codecs, chunk_shape)), "SKY");
    Require(static_cast<bool>(layout), "the layout was refused: " + (layout ? std::string{} : layout.error().message));
    return std::move(layout.value());
}

carta::zarr::Error RefusedLayout(const std::string& codecs) {
    auto layout = ParseStorageLayout(Parsed(Document(codecs)), "SKY");
    Require(!layout, "the layout was accepted when it should have been refused");
    return layout.error();
}

void TestTheCodecChainAndChunkKeyEncodingAreCarried() {
    const auto bare = Parsed(Document(""));
    Require(bare.codecs.is_array() && bare.codecs.empty(),
            "an array with no codecs should carry an empty chain, not a null");
    Require(bare.chunk_key_encoding.is_object() && bare.chunk_key_encoding.empty(),
            "an array with no chunk key encoding should carry an empty object");

    const auto chained = Parsed(Document(R"([{"name":"bytes"},{"name":"zstd"}])"));
    Require(chained.codecs.size() == 2, "the chain was not carried whole");
    Require(chained.codecs.at(1).at("name") == "zstd", "the chain was carried out of order");
}

void TestAnUnshardedArrayReportsItsOwnChunksAndCompressor() {
    const auto raw = LayoutOf("");
    Require(!raw.sharded, "an array with no sharding codec is not sharded");
    Require(raw.chunk_shape == std::vector<std::uint64_t>{4, 4}, "the chunk grid's shape is the chunk shape");
    Require(raw.compressor.empty(), "a chain with no compressor stores raw bytes");

    const auto zstd = LayoutOf(R"([{"name":"bytes"},{"name":"zstd"},{"name":"crc32c"}])");
    Require(zstd.compressor == "zstd", "the bytes-to-bytes compressor was not reported, got '" + zstd.compressor + "'");

    // crc32c is a checksum and bytes is array-to-bytes; neither is a compressor.
    const auto checksummed = LayoutOf(R"([{"name":"bytes"},{"name":"crc32c"}])");
    Require(checksummed.compressor.empty(), "a checksum was reported as a compressor");
}

void TestASharedArrayReportsTheChunksInsideTheShard() {
    // The chunk_grid describes the shard. What is decoded is the sharding codec's own chunk_shape,
    // and the compressor that applies to it is the one inside the sharding codec -- not any that
    // happens to sit beside the sharding codec in the outer chain.
    const auto layout = LayoutOf(
        R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[2,1],)"
        R"("codecs":[{"name":"bytes"},{"name":"blosc"}]}}])",
        "[8,8]");

    Require(layout.sharded, "a sharding_indexed codec was not reported as sharding");
    Require(layout.shard_shape == std::vector<std::uint64_t>{8, 8}, "the shard shape is the chunk grid's");
    Require(layout.chunk_shape == std::vector<std::uint64_t>{2, 1},
            "the inner chunk shape was not taken from the sharding codec");
    Require(layout.compressor == "blosc", "the compressor inside the sharding codec was not reported");
}

void TestTheOuterCompressorIsNotUsedForShardedChunks() {
    // A sharding codec whose configuration names no codecs leaves the compressor unset, and the
    // outer chain is then consulted -- which is what `layout.compressor.empty()` guards.
    const auto layout = LayoutOf(
        R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[2,1]}},{"name":"zstd"}])", "[8,8]");
    Require(layout.sharded && layout.chunk_shape == std::vector<std::uint64_t>{2, 1},
            "the inner chunk shape was lost");
    Require(layout.compressor == "zstd", "with no codecs inside the shard, the outer chain answers");
}

void TestAShardingCodecWithAnUnusableChunkShapeIsRefused() {
    const auto zeroed = RefusedLayout(
        R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[0,1]}}])");
    Require(zeroed.code == ErrorCode::invalid_metadata,
            "a zero inner chunk dimension should be invalid metadata");

    const auto wrong_rank = RefusedLayout(
        R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[2]}}])");
    Require(wrong_rank.code == ErrorCode::invalid_metadata,
            "an inner chunk shape of the wrong rank should be invalid metadata");

    // A sharding codec carrying no chunk_shape at all leaves the inner shape empty, which is the
    // same failure: rank zero against a rank-two shard.
    const auto absent = RefusedLayout(R"([{"name":"sharding_indexed"}])");
    Require(absent.code == ErrorCode::invalid_metadata,
            "a sharding codec with no chunk_shape should be invalid metadata");
}

}  // namespace

int main() {
    try {
        TestTheCodecChainAndChunkKeyEncodingAreCarried();
        TestAnUnshardedArrayReportsItsOwnChunksAndCompressor();
        TestASharedArrayReportsTheChunksInsideTheShard();
        TestTheOuterCompressorIsNotUsedForShardedChunks();
        TestAShardingCodecWithAnUnusableChunkShapeIsRefused();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "array metadata test failed: %s\n", error.what());
        return 1;
    }
    return 0;
}
