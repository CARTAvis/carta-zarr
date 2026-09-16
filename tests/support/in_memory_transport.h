/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_TESTS_SUPPORT_IN_MEMORY_TRANSPORT_H_
#define CARTA_ZARR_TESTS_SUPPORT_IN_MEMORY_TRANSPORT_H_

#include "zarr/transport.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace carta::zarr::testing {

/**
 * A Transport backed by a map of node name to zarr.json text; an empty node names the root.
 *
 * This is the second adapter at the transport seam, and it is what lets the schema profile's
 * structural decisions be exercised without a directory tree. It holds no array data at all, so
 * coordinate value reads report unsupported_transport: this transport serves probing and discovery,
 * never descriptors.
 *
 * Hand-written store metadata is for negative and structural cases only -- missing fields, wrong
 * types, mismatched axes, consolidated metadata that disagrees with the node listing. Positive
 * conformance is settled by generator fixtures on disk, in the schema probe tests.
 *
 * It records which nodes were asked for, because "the store did not have to read this" is a
 * requirement in its own right: consolidated metadata exists to save those reads, and the only way
 * to tell a store that used it from one that ignored it is to ask what was read.
 */
class InMemoryTransport final : public internal::Transport {
public:
    explicit InMemoryTransport(std::map<std::string, std::string> nodes) : _nodes(std::move(nodes)) {}

    Result<std::string> ReadNodeBytes(std::string_view node) const override {
        _nodes_read.insert(std::string(node));
        const auto found = _nodes.find(std::string(node));
        if (found == _nodes.end()) {
            return Error{ErrorCode::not_found, "Zarr node is missing zarr.json", std::string(node)};
        }
        return found->second;
    }

    Result<std::vector<std::string>> ListNodes() const override {
        ++_listings;
        std::vector<std::string> nodes;
        for (const auto& entry : _nodes) {
            if (!entry.first.empty()) {
                nodes.push_back(entry.first);
            }
        }
        return nodes;
    }

    Result<std::filesystem::path> ArrayDirectory(std::string_view node) const override {
        return Error{ErrorCode::unsupported_transport, "An in-memory transport holds no array data", std::string(node)};
    }

    // The nodes whose bytes were asked for, and how many times the hierarchy was listed. The root
    // is a node like any other, and reading it is what every open does first.
    const std::set<std::string>& nodes_read() const {
        return _nodes_read;
    }
    std::size_t listings() const {
        return _listings;
    }

private:
    std::map<std::string, std::string> _nodes;
    mutable std::set<std::string> _nodes_read;
    mutable std::size_t _listings = 0;
};

inline std::shared_ptr<const InMemoryTransport> MakeInMemoryTransport(std::map<std::string, std::string> nodes) {
    return std::make_shared<const InMemoryTransport>(std::move(nodes));
}

}  // namespace carta::zarr::testing

#endif  // CARTA_ZARR_TESTS_SUPPORT_IN_MEMORY_TRANSPORT_H_
