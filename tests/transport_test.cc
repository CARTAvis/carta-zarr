/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// What the filesystem transport says is in a store.
//
// Nothing pinned this. Store::Inventory sorts and de-duplicates whatever comes back, so a
// transport that named a node twice, or walked into a million chunk files to find nothing, looked
// exactly like one that did neither. These assert the listing as the transport hands it over,
// before that sort has a chance to tidy it.

#include "zarr/transport.h"

#include <algorithm>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "support/check.h"

namespace {

using carta::zarr::internal::OpenFilesystemTransport;
using carta::zarr::testing::Require;

void WriteNode(const std::filesystem::path& directory, const std::string& node_type) {
    std::filesystem::create_directories(directory);
    std::ofstream out(directory / "zarr.json");
    out << R"({"zarr_format": 3, "node_type": ")" << node_type << R"("})";
}

void WriteFile(const std::filesystem::path& file, const std::string& contents) {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream out(file);
    out << contents;
}

std::vector<std::string> Listing(const std::filesystem::path& root) {
    auto transport = OpenFilesystemTransport(root.string());
    Require(static_cast<bool>(transport), "the fixture did not open as a transport");
    auto nodes = transport.value()->ListNodes();
    Require(static_cast<bool>(nodes), "the fixture did not list");
    auto listed = std::move(nodes.value());
    std::sort(listed.begin(), listed.end());
    return listed;
}

std::string Joined(const std::vector<std::string>& nodes) {
    std::string joined;
    for (const auto& node : nodes) {
        if (!joined.empty()) {
            joined += ", ";
        }
        joined += node;
    }
    return "[" + joined + "]";
}

void RequireListing(const std::filesystem::path& root, const std::vector<std::string>& expected,
                    const std::string& what) {
    const auto listed = Listing(root);
    Require(listed == expected, what + ": listed " + Joined(listed) + ", expected " + Joined(expected));
}

// The root carries metadata like any other node, and is the one node a listing leaves out: a
// caller asked for what is in the store, not for the store.
void TestTheRootIsNotOneOfItsOwnNodes(const std::filesystem::path& root) {
    WriteNode(root, "group");
    RequireListing(root, {}, "a store holding nothing but its root");
}

// Every node once. A group's own zarr.json is a file inside a directory that was already listed, so
// a walk that counts files as well as directories reaches the same node twice.
void TestEveryNodeIsListedExactlyOnce(const std::filesystem::path& root) {
    WriteNode(root, "group");
    WriteNode(root / "GROUP", "group");
    WriteNode(root / "GROUP" / "INNER", "group");
    WriteNode(root / "SKY", "array");
    RequireListing(root, {"GROUP", "GROUP/INNER", "SKY"}, "a store with a nested group and an array");
}

// An array's chunks are its descendants and are not nodes. A real dataset has millions of them, so
// this is the difference between a listing and a full tree walk.
//
// The chunk files alone would not pin that. A chunk is neither a directory nor a zarr.json, so a
// walk that descended into the array would pass over it anyway and still answer {SKY}. What makes
// the pruning observable is a directory below the array that would otherwise qualify as a node:
// listing it is the only difference between stopping at the array and walking through it.
void TestChunksAreNotWalked(const std::filesystem::path& root) {
    WriteNode(root, "group");
    WriteNode(root / "SKY", "array");
    WriteFile(root / "SKY" / "c" / "0" / "0", "chunk bytes");
    WriteFile(root / "SKY" / "c" / "0" / "1", "chunk bytes");
    WriteNode(root / "SKY" / "BELOW", "group");
    RequireListing(root, {"SKY"}, "an array with chunks below it");
}

// A directory carrying no zarr.json is not a node, and does not stop the walk reaching the ones
// below it.
void TestADirectoryWithoutMetadataIsNotANode(const std::filesystem::path& root) {
    WriteNode(root, "group");
    std::filesystem::create_directories(root / "not_a_node");
    WriteNode(root / "not_a_node" / "DEEPER", "array");
    RequireListing(root, {"not_a_node/DEEPER"}, "a plain directory between the root and a node");
}

// An array linked in from elsewhere -- where the bytes of a large dataset as often as not are -- is
// named by where it sits in the store, not by where the link resolves: resolved, it came out as
// ../elsewhere/SKY, which the store refuses as a node name, and the dataset opened only when
// consolidated. The same goes for a store reached through a link of its own.
void TestALinkedArrayIsNamedWhereItSits(const std::filesystem::path& root) {
    const auto store = root / "store";
    WriteNode(store, "group");
    WriteNode(root / "elsewhere" / "SKY", "array");
    WriteNode(root / "elsewhere" / "GROUP", "group");
    std::filesystem::create_directory_symlink(root / "elsewhere" / "SKY", store / "SKY");
    std::filesystem::create_directory_symlink(root / "elsewhere" / "GROUP", store / "GROUP");
    RequireListing(store, {"GROUP", "SKY"}, "a store whose array and group are links");

    std::filesystem::create_directory_symlink(store, root / "alias");
    RequireListing(root / "alias", {"GROUP", "SKY"}, "the same store reached through a link");
}

}  // namespace

int main() {
    const auto root = std::filesystem::temp_directory_path() / "carta-zarr-transport-test";
    try {
        std::filesystem::remove_all(root);
        TestTheRootIsNotOneOfItsOwnNodes(root / "root-only");
        TestEveryNodeIsListedExactlyOnce(root / "nested");
        TestChunksAreNotWalked(root / "chunks");
        TestADirectoryWithoutMetadataIsNotANode(root / "plain-directory");
        TestALinkedArrayIsNamedWhereItSits(root / "linked");
        std::filesystem::remove_all(root);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "transport test failed: %s\n", error.what());
        std::filesystem::remove_all(root);
        return 1;
    }
    return 0;
}
