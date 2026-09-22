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

// What a probe's refusal should be reported as.
//
// A profile that rejects a store has already worked out why, often down to the attribute, and a
// consumer shows whatever comes back here to whoever picked the file. Reporting "not a supported
// dataset" instead throws that away and names a schema profile the store may have nothing to do
// with. The fallback is for the case the probe genuinely had nothing to say.
std::string RejectionMessage(const std::vector<Diagnostic>& diagnostics, std::string_view fallback);

// Whether a probe's answer is a dataset that can be opened, as the error a caller should report
// when it is not. The counterpart, one level up, of the question a profile answers about a single
// image within a dataset: that one is about an image, this one about the dataset holding it.
//
// Three ways it is not, and they are three different errors: nothing matched, something matched and
// was malformed, and a profile matched a store that has no images in it. Each used to be written
// out at the one call site that needed it, beside a rejection-message rule that lived there too --
// so what a probe's refusal means was decided by the facade, which is the profile's question.
Result<void> RequireOpenableDataset(const ProbeResult& probe, std::string_view location);

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

    // The two questions a profile is asked about a store it is meeting for the first time. Both go
    // through one enumeration, because deciding whether a store matches means finding out what is in
    // it -- but there is no cache behind them, so asking both costs two enumerations. Nothing wants
    // both today; ProbeStore, which does, reaches the enumeration itself.
    Result<SchemaProbeResult> Probe(const Store& store) const;
    Result<ImageDiscovery> Discover(const Store& store) const;

    // One of these rather than a verified and an unverified form. Describing an image establishes
    // its own precondition now, against the one variable it was handed rather than by enumerating
    // the store, so there is no longer a caller that has established it first and none that has to
    // be told to.
    Result<ImageDescriptor> Describe(const Store& store, std::string_view image_id) const;
    // Asks nothing about openability, deliberately: an Image handle exists only for a variable
    // Dataset::OpenImage already opened, so there is nothing left here to establish.
    Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id) const;

private:
    struct Entry {
        SchemaId id;
        Result<SchemaInspection> (*inspect)(const Store&);
        Result<ImageDescriptor> (*describe)(const Store&, std::string_view);
        Result<std::vector<Beam>> (*read_beams)(const Store&, std::string_view);
    };

    static const std::vector<Entry>& BuiltIn();

    // One enumeration, answering both halves at once. Private because no caller wants both: it is
    // what Probe and Discover are each one half of.
    Result<SchemaInspection> Inspect(const Store& store) const;

    explicit SchemaProfile(const Entry& entry) : _entry(&entry) {}

    const Entry* _entry;

    friend Result<ProbeResult> ProbeStore(const Store& store);
};

// Ask every built-in profile about a store. More than one match is reported as ambiguous rather
// than resolved silently.
Result<ProbeResult> ProbeStore(const Store& store);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_SCHEMA_PROFILE_H_
