/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_PROFILE_H_
#define CARTA_ZARR_SRC_SCHEMA_PROFILE_H_

#include "../store.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal {

// Which variables of a store a schema profile will open, and what it had to say about the ones it
// would not. Defined here rather than in store.h because it is a statement about images, and a
// Store deals in nodes and arrays; it used to sit down there only so that a cache down there could
// be keyed by profile.
struct ImageDiscovery {
    std::vector<ImageEntry> images;
    std::optional<std::string> default_image_id;
    std::vector<Diagnostic> diagnostics;
};

// Everything a profile has to say about a store it is meeting for the first time. The two halves
// come back together because they are one enumeration: deciding whether a store matches means
// finding out what is in it, and answering the two questions separately meant enumerating twice or
// caching the answer somewhere both could reach.
//
// `discovery` is what `probe` was decided from, so it is filled whatever the match kind says.
struct SchemaInspection {
    SchemaProbeResult probe;
    ImageDiscovery discovery;
};

// Whether a profile will open this image, as the error a caller should report when it will not.
//
// Shared because the facade asks it of a descriptor it already holds and a profile asks it of a
// store it has just inspected, and the two must give the same answer. A variable that was listed
// but will not open is a different answer from one that was never seen, and a caller can act on the
// difference.
Result<void> RequireOpenable(const std::vector<ImageEntry>& images, std::string_view image_id);

/**
 * A named, versioned description of how an image dataset is laid out, bound to its identifier.
 *
 * A profile is looked up once and then asked questions, rather than every question carrying an
 * identifier to be resolved again. It binds the identifier only: the built-in table outlives every
 * store, so binding a store here would invent a lifetime problem that does not exist.
 *
 * "Profile" rather than "adapter" deliberately -- CONTEXT.md names this concept the schema profile,
 * and an adapter in this codebase is the thing satisfying an interface at a seam, as the filesystem
 * and in-memory transports do.
 */
class SchemaProfile {
public:
    // Reports unsupported_schema when no built-in profile carries this identifier.
    static Result<SchemaProfile> For(std::string_view schema_id);

    const SchemaId& id() const noexcept;

    // One enumeration of the store, answering both of the questions a profile is asked about it.
    Result<SchemaInspection> Inspect(const Store& store) const;
    // Halves of an Inspect, for a caller that wants only one of them. Each one inspects: there is
    // no cache behind them, so a caller wanting both should ask once.
    Result<SchemaProbeResult> Probe(const Store& store) const;
    Result<ImageDiscovery> Discover(const Store& store) const;

    // Inspects the store first, to ask whether the profile will open this image at all.
    Result<ImageDescriptor> Describe(const Store& store, std::string_view image_id) const;
    // For a caller that has already established openability against a discovery it kept -- which is
    // what Dataset does, so that opening an image does not re-enumerate the store.
    Result<ImageDescriptor> DescribeVerified(const Store& store, std::string_view image_id) const;
    // Verified in the same sense: an Image handle only exists for a variable Dataset::OpenImage
    // already accepted, so asking again would enumerate the whole store on every call.
    Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id) const;

private:
    struct Entry {
        SchemaId id;
        Result<SchemaInspection> (*inspect)(const Store&);
        Result<ImageDescriptor> (*describe)(const Store&, std::string_view);
        Result<std::vector<Beam>> (*read_beams)(const Store&, std::string_view);
    };

    static const std::vector<Entry>& BuiltIn();

    explicit SchemaProfile(const Entry& entry) : _entry(&entry) {}

    const Entry* _entry;

    friend Result<ProbeResult> ProbeStore(const Store& store);
};

// Ask every built-in profile about a store. More than one match is reported as ambiguous rather
// than resolved silently.
Result<ProbeResult> ProbeStore(const Store& store);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_SCHEMA_PROFILE_H_
