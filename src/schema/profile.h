/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_PROFILE_H_
#define CARTA_ZARR_SRC_SCHEMA_PROFILE_H_

#include "../store.h"
#include "chunk_geometry.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal {

// One image as a profile describes it: what a consumer is told about it, and the geometry its reads
// and reductions are cut by. Made together because the geometry is the image's stored layout in the
// descriptor's logical order, and the layout is read only while the image is described.
struct DescribedImage {
    ImageDescriptor descriptor;
    ChunkGeometry geometry;
    // The flag's geometry, from its own layout, which may differ from the pixels': a flag has to share
    // the image's dimensions and nothing more. Empty when the image has no pixel mask.
    ChunkGeometry flag_geometry;
};

// What every built-in profile made of one store, as ProbeStore reports it: whether one of them
// describes it, and if so the dataset it describes.
enum class ProbeKind {
    zarr_without_supported_schema,
    supported_dataset,
    invalid_dataset,
};

// A probe answers one more question than a descriptor does, whether a profile describes the store at
// all, and is otherwise exactly what Dataset::descriptor() reports; Dataset::Open moves it from here.
// Inheriting rather than copying the fields means a field added to DatasetDescriptor reaches both.
//
// The inherited fields are filled to the extent the kind allows. A supported_dataset fills them all;
// an invalid_dataset names the schema it failed to be and says why in diagnostics; anything else
// carries diagnostics alone.
struct ProbeResult : DatasetDescriptor {
    ProbeKind kind = ProbeKind::zarr_without_supported_schema;
};

// Which variables of a store a schema profile will open, and what it had to say about the ones it
// would not. Defined here rather than in store.h because it is a statement about images, and a
// Store deals in nodes and arrays.
struct ImageDiscovery {
    std::vector<ImageEntry> images;
    std::optional<std::string> default_image_id;
    std::vector<Diagnostic> diagnostics;
};

// Everything a profile has to say about a store it is meeting for the first time. The two halves
// come back together because they come from one enumeration: deciding whether a store matches means
// finding out what is in it.
//
// `discovery` is what `probe` was decided from, so it is filled whatever the match kind says.
struct SchemaInspection {
    SchemaProbeResult probe;
    ImageDiscovery discovery;
};

// What a probe's refusal should be reported as.
//
// A profile that rejects a store has already worked out why, often down to the attribute, and a
// consumer shows this message to whoever picked the file. "Not a supported dataset" would throw that
// away and name a schema profile the store may have nothing to do with. The fallback is for when the
// probe had nothing to say.
std::string RejectionMessage(const std::vector<Diagnostic>& diagnostics, std::string_view fallback);

// Whether a probe's answer is a dataset that can be opened, as the error a caller should report when
// it is not. The dataset-level counterpart of RequireQualified, which asks it of one image.
//
// There are three ways it is not, each a different error: nothing matched, something matched and was
// malformed, or a profile matched a store with no images in it. The profile, not the facade, decides
// which one a refusal is and what it says.
Result<void> RequireOpenableDataset(const ProbeResult& probe, std::string_view location);

/**
 * A named, versioned description of how an image dataset is laid out, bound to its identifier.
 *
 * A profile is looked up once and then asked questions, rather than each question resolving an
 * identifier again. It binds only the identifier: the built-in table outlives every store, so
 * binding a store here would create a lifetime problem.
 *
 * "Profile" rather than "adapter": GLOSSARY.md calls this concept the schema profile, and an adapter
 * in this codebase is what implements an interface at a seam, as the filesystem and in-memory
 * transports do.
 */
class SchemaProfile {
public:
    // Reports unsupported_schema when no built-in profile carries this identifier.
    static Result<SchemaProfile> For(std::string_view schema_id);

    // The two questions a profile is asked about a store it is meeting for the first time. Each does
    // one enumeration with no cache behind it, so asking both costs two. ProbeStore, which needs
    // both, calls Inspect instead.
    Result<SchemaProbeResult> Probe(const Store& store) const;
    Result<ImageDiscovery> Discover(const Store& store) const;

    // Describing an image checks its own precondition against the one variable it was handed,
    // without enumerating the store, so no caller has to check it first.
    Result<DescribedImage> Describe(const Store& store, std::string_view image_id) const;
    // Does not check openability: an Image handle exists only for a variable Dataset::OpenImage
    // already opened.
    Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id) const;

private:
    struct Entry {
        SchemaId id;
        Result<SchemaInspection> (*inspect)(const Store&);
        Result<DescribedImage> (*describe)(const Store&, std::string_view);
        Result<std::vector<Beam>> (*read_beams)(const Store&, std::string_view);
    };

    static const std::vector<Entry>& BuiltIn();

    // One enumeration, answering both halves at once. Probe and Discover each return one half.
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
