/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_TRANSPORT_H_
#define CARTA_ZARR_SRC_ZARR_TRANSPORT_H_

#include "carta-zarr/result.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal {

// The seam between a Zarr hierarchy and where its bytes live. A Transport supplies raw access: it
// hands a node's document up verbatim and never decides what a node means. Keeping interpretation
// above the seam is deliberate -- it is what stops a test transport and the filesystem transport
// from disagreeing about the very metadata the tests exist to pin down.
//
// One transport does parse, and it is worth naming rather than denying: FilesystemTransport reads
// `node_type` out of a node's document while listing, to know whether to descend into it. The
// alternative is walking the chunk files, of which a real dataset has millions. It parses to decide
// where to walk and never to decide what anything means, and the document it hands up is the bytes
// it found; a document it cannot parse is passed on for Store to diagnose.
//
// Node names reaching a Transport have already been validated and normalized by Store -- including
// through Store::ResolveArrayDirectory, which did not always do so -- so every Transport is held to
// one rule about what a legal node path is. A Transport may check again; it must not have a second
// opinion.
class Transport {
public:
    Transport() = default;
    Transport(const Transport&) = delete;
    Transport& operator=(const Transport&) = delete;
    Transport(Transport&&) = delete;
    Transport& operator=(Transport&&) = delete;
    virtual ~Transport() = default;

    // Read one node's zarr.json verbatim. An empty node names the root. Reports not_found when the
    // node carries no metadata, and io_error for every other failure.
    virtual Result<std::string> ReadNodeBytes(std::string_view node) const = 0;

    // Every node in the hierarchy, excluding the root, in any order. Only consulted when the store
    // carries no consolidated metadata.
    virtual Result<std::vector<std::string>> ListNodes() const = 0;

    // Where a node's array lives, as a filesystem path. This is the one place the seam is not
    // about bytes: a transport holding no filesystem data reports unsupported_transport rather than
    // pretending TensorStore or the string decoder can reach it, which means such a transport can
    // serve a probe and a discovery but cannot open an image. ADR 0004 records what that costs and
    // why the seam still has this shape.
    //
    // The path is absolute and does not depend on the process's working directory: an image is read
    // long after it is opened, and a store that answered relatively would hand out locations that
    // stop meaning the same thing. Nothing above or below this normalizes it again.
    virtual Result<std::filesystem::path> ArrayDirectory(std::string_view node) const = 0;

    // How many bytes this store occupies where it lives, counting everything: chunks, shards and
    // metadata, not only what the arrays logically hold.
    //
    // Here rather than above the seam because it is the one question about a store that only the
    // place the bytes live can answer. A facade that walked a directory itself would be holding a
    // second opinion about what a location means, and a weaker one -- it would not have the
    // absolute path OpenFilesystemTransport resolved, so a process that changed directory between
    // opening the store and asking its size would measure somewhere else.
    //
    // Reports `cancelled` when the deadline passes mid-walk, and `unsupported_transport` when the
    // transport has no stored bytes to count. Pure virtual for the same reason ArrayDirectory is: a
    // transport that cannot answer says so rather than serving a number that looks measured.
    //
    // Every failure here means the caller gets the declared size instead, so this is allowed to
    // give up. It is a size, not a read.
    virtual Result<std::uint64_t> StoredSizeBytes(std::chrono::steady_clock::time_point deadline) const = 0;
};

using TransportPtr = std::shared_ptr<const Transport>;

// Normalize a location (a bare path, or file://) and check that it names a directory.
Result<TransportPtr> OpenFilesystemTransport(std::string_view location);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_ZARR_TRANSPORT_H_
