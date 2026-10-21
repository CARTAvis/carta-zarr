/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The schema profile's structural decisions -- does this store match, which variables are images,
// which of them are openable -- read metadata and nothing else. They are exercised here through an
// in-memory transport, so a case costs a map entry rather than a directory tree, and the build links
// no TensorStore.
//
// Whole descriptors are deliberately absent: they read coordinate values, which an in-memory
// transport has none of. Descriptor behaviour is covered in tests/schema_probe_test.cc against
// fixtures on disk. The parts of describing an image that read only metadata do belong here, which
// is why choosing a flag is below: what a flag has to be is checked in tests/flag_test.cc without a
// store at all, and what happens when the store is consulted is checked here.

#include "schema/profile.h"
#include "schema/xradio/coordinates.h"
#include "schema/xradio/flag.h"
#include "schema/xradio/image.h"
#include "store.h"
#include "support/check.h"
#include "support/in_memory_transport.h"
#include "zarr/array_metadata.h"
#include "zarr/pixel_selection.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::SchemaMatchKind;
using carta::zarr::internal::ProbeKind;
using carta::zarr::internal::xradio::DetermineFlag;
using carta::zarr::testing::MakeInMemoryTransport;

using carta::zarr::testing::Require;

bool HasDiagnostic(const carta::zarr::SchemaProbeResult& probe, carta::zarr::DiagnosticCode code);

bool HasDiagnostic(const std::vector<carta::zarr::Diagnostic>& diagnostics, carta::zarr::DiagnosticCode code) {
    for (const auto& diagnostic : diagnostics) {
        if (diagnostic.code == code) {
            return true;
        }
    }
    return false;
}

bool HasDiagnostic(const carta::zarr::SchemaProbeResult& probe, carta::zarr::DiagnosticCode code) {
    return HasDiagnostic(probe.diagnostics, code);
}

std::vector<std::string> ImageIds(const std::vector<carta::zarr::ImageEntry>& images) {
    std::vector<std::string> ids;
    ids.reserve(images.size());
    for (const auto& image : images) {
        ids.push_back(image.id);
    }
    return ids;
}

std::vector<std::string> OpenableImageIds(const std::vector<carta::zarr::ImageEntry>& images) {
    std::vector<std::string> ids;
    for (const auto& image : images) {
        if (image.openable) {
            ids.push_back(image.id);
        }
    }
    return ids;
}

// A root group with `attributes`, and the coordinate system every store here carries unless it is
// left out.
std::string RootWith(const std::string& attributes, bool coordinate_system = true) {
    std::string all = attributes;
    if (coordinate_system) {
        all += std::string(all.empty() ? "" : ",") + R"("coordinate_system_info": {
      "projection": "SIN",
      "reference_direction": {"data": [1.0, 0.5]},
      "native_pole_direction": {"data": [0.0, 1.5707963267948966]},
      "pixel_coordinate_transformation_matrix": [[1.0, 0.0], [0.0, 1.0]]
    })";
    }
    return "{\"attributes\":{" + all + "},\"zarr_format\":3,\"node_type\":\"group\"}";
}

// The root of a store of schema 0.0.2, whose base data group names SKY as its sky.
std::string RootGroup(bool coordinate_system = true) {
    return RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"sky":"SKY"}})", coordinate_system);
}

std::string NumericArray(const std::string& shape, const std::string& dimensions,
                         const std::string& data_type = "float64", const std::string& attributes = "{}") {
    return "{\"shape\":" + shape + ",\"data_type\":\"" + data_type +
           "\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":" + shape +
           "}},\"attributes\":" + attributes + ",\"dimension_names\":" + dimensions +
           ",\"zarr_format\":3,\"node_type\":\"array\"}";
}

std::string SkyArray(const std::string& data_type = "float32",
                     const std::string& attributes = R"({"units":"Jy/beam","type":"sky"})",
                     const std::string& dimensions = R"(["time","frequency","polarization","l","m"])") {
    return NumericArray("[1,3,2,4,5]", dimensions, data_type, attributes);
}

std::string PolarizationArray(const std::string& data_type = R"({"name": "fixed_length_utf32",
                                    "configuration": {"length_bytes": 4}})") {
    return "{\"shape\":[2],\"data_type\":" + data_type +
           ",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":[2]}},"
           "\"attributes\":{},\"dimension_names\":[\"polarization\"],\"zarr_format\":3,\"node_type\":\"array\"}";
}

// A spectral coordinate as schema 0.0.2 writes it: with its frame and its units.
std::string FrequencyArray(const std::string& shape = "[3]",
                           const std::string& attributes = R"({"units":"Hz","frame":"LSRK"})") {
    return NumericArray(shape, R"(["frequency"])", "float64", attributes);
}

// A complete image dataset of schema 0.0.2 carries a coordinate array for every axis its image uses,
// time included.
std::map<std::string, std::string> CompleteStore() {
    return {
        {"", RootGroup()},
        {"SKY", SkyArray()},
        {"time", NumericArray("[1]", R"(["time"])")},
        {"frequency", FrequencyArray()},
        {"polarization", PolarizationArray()},
        {"l", NumericArray("[4]", R"(["l"])")},
        {"m", NumericArray("[5]", R"(["m"])")},
    };
}

// The same store with a copy of every child's metadata in the root, which is what zarr-python
// writes when a dataset is consolidated. The children keep their own metadata: consolidated
// metadata is a copy that saves reading them, carrying `must_understand: false` precisely so that a
// reader which ignores it still reads the same hierarchy. A store whose children exist only in the
// root is malformed, and no reader but this one could open it -- the array data behind those names
// is read by TensorStore, which needs each array's own metadata.
//
// `copies` replaces what the root's copy says about a node, leaving the node's own document as it
// was: the store a rewrite leaves behind when it does not consolidate again.
std::map<std::string, std::string> ConsolidatedStore(std::map<std::string, std::string> nodes = CompleteStore(),
                                                     const std::map<std::string, std::string>& copies = {}) {
    std::string metadata;
    for (const auto& child : nodes) {
        if (child.first.empty()) {
            continue;
        }
        if (!metadata.empty()) {
            metadata += ",";
        }
        const auto copy = copies.find(child.first);
        metadata += "\"" + child.first + "\":" + (copy != copies.end() ? copy->second : child.second);
    }
    nodes[""] =
        "{\"attributes\":{\"schema_version\":\"0.0.2\",\"data_groups\":{\"base\":{\"sky\":\"SKY\"}},"
        "\"coordinate_system_info\":{"
        "\"projection\":\"SIN\","
        "\"reference_direction\":{\"data\":[1.0,0.5]},"
        "\"native_pole_direction\":{\"data\":[0.0,1.5707963267948966]},"
        "\"pixel_coordinate_transformation_matrix\":[[1.0,0.0],[0.0,1.0]]}},"
        "\"zarr_format\":3,\"node_type\":\"group\","
        "\"consolidated_metadata\":{\"kind\":\"inline\",\"must_understand\":false,\"metadata\":{" +
        metadata + "}}}";
    return nodes;
}

carta::zarr::internal::SchemaProfile XradioProfile() {
    auto profile = carta::zarr::internal::SchemaProfile::For(carta::zarr::kXradioImageSchema);
    Require(static_cast<bool>(profile), "the built-in XRADIO profile was not found");
    return profile.value();
}

carta::zarr::Result<carta::zarr::internal::Store> Open(std::map<std::string, std::string> nodes) {
    return carta::zarr::internal::OpenStore(MakeInMemoryTransport(std::move(nodes)));
}

carta::zarr::SchemaProbeResult Probe(std::map<std::string, std::string> nodes) {
    auto store = Open(std::move(nodes));
    Require(static_cast<bool>(store), "OpenStore rejected the in-memory store");
    auto probe = XradioProfile().Probe(store.value());
    Require(static_cast<bool>(probe), "the profile probe reported an error");
    return probe.value();
}

void TestCompleteStoreMatches() {
    const auto probe = Probe(CompleteStore());
    Require(probe.kind == SchemaMatchKind::match, "a complete image dataset did not match");
    Require(probe.schema_version == "0.0.2", "a 0.0.2 store reported '" + probe.schema_version + "'");
    Require(probe.diagnostics.empty(), "a complete image dataset was matched with something to say");
}

void TestDatasetWithoutSkyMatches() {
    auto nodes = CompleteStore();
    nodes.erase("SKY");
    nodes[""] = RootWith(R"("schema_version":"0.0.2","data_groups":{"residual":{"sky":"SKY_RESIDUAL"}})");
    nodes["SKY_RESIDUAL"] = SkyArray();
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::match, "an image dataset with only SKY_RESIDUAL was incorrectly rejected");
}

void TestNonMatch() {
    auto store = Open({{"", RootGroup()}, {"OTHER", NumericArray("[2]", R"(["x"])")}});
    Require(static_cast<bool>(store), "a valid non-XRADIO group failed to open");
    const auto probe = carta::zarr::internal::ProbeStore(store.value());
    Require(static_cast<bool>(probe), "ProbeStore reported an error for a valid non-XRADIO group");
    Require(probe.value().kind == ProbeKind::zarr_without_supported_schema,
            "a valid non-XRADIO group was not reported as an unsupported schema");
    // Nothing in it looked like an image, so the profile has nothing to say about it either.
    Require(probe.value().diagnostics.empty() && probe.value().schema_id.empty(),
            "a store with nothing the profile recognised was reported with the profile's say-so");
}

// A store whose only images are ones this library does not open is still one no profile matched,
// and the reason is the profile's own: it recognised the images and refused them. That reason used
// to be dropped on the way out of ProbeStore, so the refusal a consumer saw was "no profile matched"
// with nothing to say which variable, or why.
void TestNothingOpenableSaysWhy() {
    using carta::zarr::internal::ProbeStore;
    using carta::zarr::internal::RequireOpenableDataset;

    const auto refused = [](std::map<std::string, std::string> nodes, carta::zarr::DiagnosticCode code,
                            const std::string& message) {
        auto store = Open(std::move(nodes));
        Require(static_cast<bool>(store), "the unopenable store failed to open");
        const auto probe = ProbeStore(store.value());
        Require(static_cast<bool>(probe), "ProbeStore reported an error for a well-formed store");
        Require(probe.value().kind == ProbeKind::zarr_without_supported_schema,
                "a store of variables this library does not open was not reported as unsupported");
        // Diagnostics alone: naming the schema would say the store is one.
        Require(probe.value().schema_id.empty(), "an unmatched store was attributed to a schema");
        Require(HasDiagnostic(probe.value().diagnostics, code) && probe.value().diagnostics.front().node_path == "SKY",
                "an unmatched store lost the reason its images were refused");

        const auto openable = RequireOpenableDataset(probe.value(), "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::unsupported_schema,
                "a store of variables this library does not open was not refused as an unsupported schema");
        Require(openable.error().message == message,
                "the refusal did not carry the profile's reason: " + openable.error().message);
    };

    auto complex = CompleteStore();
    complex["SKY"] = SkyArray("complex64");
    refused(complex, carta::zarr::DiagnosticCode::unsupported_data_type,
            "Complex sky-plane variables are not openable");

    // An image of booleans holds a flag's values, and is said to, not called complex.
    auto boolean = CompleteStore();
    boolean["SKY"] = SkyArray("bool");
    refused(boolean, carta::zarr::DiagnosticCode::unsupported_data_type,
            "Boolean sky-plane variables are not openable; a flag is typed 'flag'");

    // An int8 is a number, whatever its attributes record: typed as a sky image, it is one, and opens.
    auto int8_sky = CompleteStore();
    int8_sky["SKY"] = SkyArray("int8", R"({"type":"sky","dtype":"bool"})");
    auto int8_store = Open(int8_sky);
    Require(static_cast<bool>(int8_store), "the in-memory store failed to open");
    const auto int8_listing = XradioProfile().Discover(int8_store.value());
    Require(static_cast<bool>(int8_listing) &&
                OpenableImageIds(int8_listing.value().images) == std::vector<std::string>{"SKY"},
            "an int8 sky image recorded as bool was not listed openable as the number it is");
    Require(Probe(int8_sky).kind == SchemaMatchKind::match, "a store of an int8 sky image did not match");

    // A real image typed for the aperture plane opens (ADR 0024); a complex one does not, on either plane.
    // Named SKY only so that the reason is found where the others' is; typed and grouped as what it is.
    auto aperture = CompleteStore();
    aperture[""] = RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"aperture":"SKY"}})");
    aperture["SKY"] = SkyArray("complex64", R"({"type":"aperture"})", R"(["time","frequency","polarization","u","v"])");
    refused(aperture, carta::zarr::DiagnosticCode::unsupported_data_type,
            "Complex aperture-plane variables are not openable");
}

void TestMissingCoordinateIsInvalid() {
    auto nodes = CompleteStore();
    nodes.erase("time");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid,
            "an image dataset missing its time coordinate was not reported as invalid");
    Require(HasDiagnostic(probe.diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
            "the missing coordinate produced no diagnostic");
}

void TestCoordinateShapeMismatchIsInvalid() {
    auto nodes = CompleteStore();
    nodes["frequency"] = FrequencyArray("[7]");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a coordinate whose length disagrees with the image was accepted");
}

void TestCoordinateDataTypeIsChecked() {
    auto nodes = CompleteStore();
    // Polarization labels are strings; a numeric polarization coordinate is malformed.
    nodes["polarization"] = PolarizationArray("\"float64\"");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a numeric polarization coordinate was accepted");
    Require(HasDiagnostic(probe.diagnostics, carta::zarr::DiagnosticCode::unsupported_data_type),
            "the polarization data type produced no diagnostic");

    // And a number has to be a real one. Labels are also held to being decodable, which refuses a
    // numeric polarization on its own; nothing but the type check refuses a complex frequency.
    auto complex = CompleteStore();
    complex["frequency"] = NumericArray("[3]", R"(["frequency"])", "complex64", R"({"units":"Hz","frame":"LSRK"})");
    const auto complex_probe = Probe(complex);
    Require(complex_probe.kind == SchemaMatchKind::invalid, "a complex frequency coordinate was accepted");
    Require(HasDiagnostic(complex_probe.diagnostics, carta::zarr::DiagnosticCode::unsupported_data_type),
            "the frequency data type produced no diagnostic");
}

// Labels this library cannot decode are refused when the dataset is probed, not when an image is
// opened: the probe is what a listing of openable images rests on, so a polarization coordinate in
// a codec the string reader does not run would otherwise list SKY as openable and fail it on open.
void TestUndecodableLabelsAreInvalid() {
    auto nodes = CompleteStore();
    nodes["polarization"] =
        R"({"shape":[2],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":4}},)"
        R"("chunk_grid":{"name":"regular","configuration":{"chunk_shape":[2]}},)"
        R"("codecs":[{"name":"bytes","configuration":{"endian":"little"}},{"name":"bz2"}],)"
        R"("attributes":{},"dimension_names":["polarization"],"zarr_format":3,"node_type":"array"})";
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid,
            "polarization labels in a codec this library cannot run were accepted");
    Require(HasDiagnostic(probe.diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
            "the undecodable polarization labels produced no diagnostic");
}

void TestMalformedCoordinateSystemIsInvalid() {
    auto nodes = CompleteStore();
    nodes[""] = R"({"attributes":{"coordinate_system_info":{}},"zarr_format":3,"node_type":"group"})";
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "malformed coordinate_system_info was accepted");
}

// Every XRADIO writes coordinate_system_info, so a store without it is malformed too, and refused
// rather than opened with a direction of "" at (0, 0).
void TestAMissingCoordinateSystemIsInvalid() {
    auto nodes = CompleteStore();
    nodes[""] = RootGroup(false);
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a store with no coordinate_system_info was accepted");
    Require(!probe.diagnostics.empty() && probe.diagnostics.front().node_path == "/attributes/coordinate_system_info",
            "the refusal should name the attribute that is missing");
}

// A variable typed as a sky image is held to the sky plane's axes. One without them is not quietly
// passed over as something else: the store says SKY is an image and gives it axes it cannot be read
// by, which is the store being wrong about itself.
void TestAnImageWithoutItsPlanesAxesIsInvalid() {
    const auto refused = [](std::map<std::string, std::string> nodes, const std::string& what) {
        const auto probe = Probe(nodes);
        Require(probe.kind == SchemaMatchKind::invalid, what + " was not refused as invalid");
        Require(!probe.diagnostics.empty() && probe.diagnostics.front().node_path == "SKY" &&
                    probe.diagnostics.front().code == carta::zarr::DiagnosticCode::invalid_metadata,
                what + " was refused without saying SKY is why");
    };
    auto missing_axis = CompleteStore();
    missing_axis["SKY"] =
        SkyArray("float32", R"({"units":"Jy/beam","type":"sky"})", R"(["time","frequency","polarization","l","x"])");
    refused(missing_axis, "a SKY with x in place of m");

    auto wrong_rank = CompleteStore();
    wrong_rank["SKY"] = NumericArray("[1,3,2,4]", R"(["time","frequency","polarization","l"])", "float32",
                                     R"({"units":"Jy/beam","type":"sky"})");
    refused(wrong_rank, "a four-dimensional SKY");
}

// Discovery reports every variable typed as an image, and says why the ones it cannot open are
// closed. Flags, the optional (l, m) coordinates and an untyped variable over a whole plane's axes
// are never images at all: what an image is, is its type (ADR 0001).
void TestDiscoveryClassifiesVariables() {
    auto nodes = CompleteStore();
    nodes["SKY_MODEL"] = SkyArray();
    nodes["COMPLEX"] = SkyArray("complex64");
    nodes["APERTURE"] = SkyArray("float32", R"({"type":"aperture"})", R"(["time","frequency","polarization","u","v"])");
    nodes["VISIBILITY"] =
        SkyArray("complex64", R"({"type":"visibility"})", R"(["time","frequency","polarization","u","v"])");
    nodes["MASK_0"] = SkyArray("bool", R"({"type":"flag"})");
    nodes["UNTYPED"] = SkyArray("float32", "{}");
    nodes["PRIMARY_BEAM"] = SkyArray("float32", R"({"type":"primary_beam"})");
    nodes["UV_SAMPLING"] =
        SkyArray("float32", R"({"type":"uv_sampling"})", R"(["time","frequency","polarization","u","v"])");
    // Typed as an image of the other plane: listed, and malformed.
    nodes["SKY_ON_UV"] = SkyArray("float32", R"({"type":"sky"})", R"(["time","frequency","polarization","u","v"])");
    nodes["APERTURE_ON_LM"] = SkyArray("float32", R"({"type":"aperture"})");
    nodes["right_ascension"] = NumericArray("[4,5]", R"(["l","m"])");
    nodes["u"] = NumericArray("[4]", R"(["u"])");
    nodes["v"] = NumericArray("[5]", R"(["v"])");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the discovery store failed to open");
    auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "profile discovery reported an error");

    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY", "SKY_MODEL", "PRIMARY_BEAM",
                                                                           "APERTURE", "APERTURE_ON_LM", "COMPLEX",
                                                                           "SKY_ON_UV", "UV_SAMPLING", "VISIBILITY"},
            "discovery did not enumerate the images in display order");
    Require(OpenableImageIds(discovery.value().images) ==
                std::vector<std::string>{"SKY", "SKY_MODEL", "PRIMARY_BEAM", "APERTURE", "UV_SAMPLING"},
            "discovery did not restrict the openable images to real images of either plane");
    for (const auto& id : {"SKY_ON_UV", "APERTURE_ON_LM"}) {
        const auto& images = discovery.value().images;
        const auto entry =
            std::find_if(images.begin(), images.end(), [&](const auto& image) { return image.id == id; });
        Require(
            entry != images.end() && HasDiagnostic(entry->diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
            std::string(id) + ", typed as an image of the other plane, was not listed as malformed");
    }
    Require(discovery.value().default_image_id == "SKY", "discovery did not select the default openable image");
    for (const auto& id : {"COMPLEX", "VISIBILITY"}) {
        const auto& images = discovery.value().images;
        const auto entry =
            std::find_if(images.begin(), images.end(), [&](const auto& image) { return image.id == id; });
        Require(entry != images.end() &&
                    HasDiagnostic(entry->diagnostics, carta::zarr::DiagnosticCode::unsupported_data_type),
                std::string("the complex variable ") + id + " produced no diagnostic");
    }
}

// The openable gate sits on the profile handle, ahead of any descriptor work.
void TestOpenableGate() {
    auto nodes = CompleteStore();
    nodes["COMPLEX"] = SkyArray("complex64");
    nodes["MASK_0"] = SkyArray("bool", R"({"type":"flag"})");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the openable-gate store failed to open");
    auto profile = carta::zarr::internal::SchemaProfile::For(carta::zarr::kXradioImageSchema);
    Require(static_cast<bool>(profile), "the built-in XRADIO profile was not found");

    const auto complex = profile.value().Describe(store.value(), "COMPLEX");
    Require(!complex && complex.error().code == ErrorCode::unsupported_data_type,
            "a complex variable was not reported as unopenable");

    const auto flag = profile.value().Describe(store.value(), "MASK_0");
    Require(!flag && flag.error().code == ErrorCode::not_found, "a flag variable was exposed as an image");

    const auto absent = profile.value().Describe(store.value(), "NOPE");
    Require(!absent && absent.error().code == ErrorCode::not_found, "an absent variable was not reported as missing");
}

// Coordinates belong to the dataset and every image references them by dimension name, so an image
// whose own extent disagrees with the coordinate it names cannot be described with it: it would be
// reported with three channels' worth of coordinates over seven channels of pixels. It is not
// listed as openable.
void TestAnImageDisagreeingWithACoordinateIsNotOpenable() {
    auto nodes = CompleteStore();
    // Seven channels of pixels where the dataset's frequency coordinate has three.
    nodes["MODEL"] = NumericArray("[1,7,2,4,5]", R"(["time","frequency","polarization","l","m"])", "float32",
                                  R"({"units":"Jy/beam","type":"sky"})");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the disagreeing-image store failed to open");
    const auto profile = XradioProfile();
    const auto discovery = profile.Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed on the disagreeing-image store");
    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY", "MODEL"},
            "the disagreeing image was dropped from the listing rather than listed with its reason");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "an image that will not open was listed as openable");
    Require(discovery.value().default_image_id == "SKY", "the default image was not the first openable one");

    const auto model = profile.Describe(store.value(), "MODEL");
    Require(!model && model.error().code == ErrorCode::invalid_metadata,
            "a store that is wrong about an image was refused as a library limitation");
}

// Nothing openable is two different answers, and a consumer acts on the difference. A store of
// well-formed variables this library does not open is a store it is not for; a store whose only
// image disagrees with its coordinates is one this profile recognised and found broken, and
// reporting that as "no profile matched" names nothing anyone can act on.
void TestNothingOpenableIsInvalidOnlyWhenSomethingIsMalformed() {
    auto capability = CompleteStore();
    capability["SKY"] = SkyArray("complex64");
    Require(Probe(capability).kind == SchemaMatchKind::no_match,
            "a store of variables this library does not open was reported as broken");

    auto malformed = CompleteStore();
    malformed["SKY"] = NumericArray("[1,7,2,4,5]", R"(["time","frequency","polarization","l","m"])", "float32",
                                    R"({"units":"Jy/beam","type":"sky"})");
    // Enumerated before SKY, and a limitation rather than a defect: without the reason for the
    // refusal being put first, this is the diagnostic a consumer would be shown.
    malformed["APERTURE"] = NumericArray("[1,3,2,4,5]", R"(["time","frequency","polarization","u","v"])", "complex64",
                                         R"({"type":"aperture"})");

    const auto broken = Probe(malformed);
    Require(broken.kind == SchemaMatchKind::invalid,
            "a store whose only image disagrees with a coordinate was reported as unrecognised");
    Require(broken.diagnostics.size() == 2 && broken.diagnostics.front().node_path == "SKY" &&
                broken.diagnostics.back().node_path == "APERTURE",
            "the refusal did not lead with the reason the store was refused, and keep the other");
}

// A group is a node in the hierarchy, and a store is free to hold one. It was being read as an
// array, failing, and reported as a node that could not be read -- so every store with a nested
// group carried a diagnostic about a node that was perfectly well formed, and a promotion rule that
// counted such a diagnostic would close a store over it.
void TestANestedGroupIsNotABrokenArray() {
    auto nodes = CompleteStore();
    nodes["SUBDIR"] = RootGroup(false);

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the nested-group store failed to open");
    const auto profile = XradioProfile();
    const auto discovery = profile.Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed on the nested-group store");
    Require(!HasDiagnostic(discovery.value().diagnostics, carta::zarr::DiagnosticCode::unreadable_array),
            "a valid group was diagnosed as an array that could not be read");
    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "a group was listed among the dataset's images");

    // Naming one is naming something that is not an image, which is the answer a flag and a
    // coordinate get as well.
    const auto group = profile.Describe(store.value(), "SUBDIR");
    Require(!group && group.error().code == ErrorCode::not_found,
            "a group was refused as a broken array rather than as something that is not an image");
}

// A refusal names the variable's own reason rather than answering generically.
//
// The reason an image is not openable is the same whichever way a consumer arrives at it. One rule
// decides openability; what is asserted here is the half a single rule does not guarantee on its
// own -- that the reason survives the trip from the listing to the refusal.
void TestARefusalCarriesTheVariablesOwnReason() {
    auto nodes = CompleteStore();
    nodes["COMPLEX"] = SkyArray("complex64");
    // Valid JSON, and not array metadata: node_type says array and there is no shape.
    nodes["BROKEN"] = R"({"zarr_format":3,"node_type":"array","data_type":"float32"})";

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the refusal store failed to open");
    const auto profile = XradioProfile();
    const auto discovery = profile.Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed on the refusal store");
    const auto& images = discovery.value().images;

    const auto complex = profile.Describe(store.value(), "COMPLEX");
    Require(!complex && complex.error().code == ErrorCode::unsupported_data_type,
            "a listed but unopenable variable should be refused as an unsupported data type");

    const auto listed =
        std::find_if(images.begin(), images.end(), [](const auto& image) { return image.id == "COMPLEX"; });
    Require(listed != images.end() && !listed->diagnostics.empty(),
            "the complex variable was not listed with a diagnostic to pass on");
    Require(complex.error().message == listed->diagnostics.front().message,
            "the refusal did not carry the variable's own diagnostic");

    const auto absent = profile.Describe(store.value(), "NOPE");
    Require(!absent && absent.error().code == ErrorCode::not_found,
            "a variable that is not there should be reported as missing, not as unopenable");

    // A node that is there and will not parse is neither of those, and not "not found" either,
    // which is a name the dataset does not have.
    const auto broken = profile.Describe(store.value(), "BROKEN");
    Require(!broken && broken.error().code == ErrorCode::invalid_metadata,
            "a node that is present and will not parse was reported as a missing variable");
}

// The dataset-level counterpart, and the three answers it has to keep apart: what a probe's refusal
// means is the profile's to say, and each kind is reached here without a store on disk shaped to
// produce it.
void TestOneRuleDecidesWhatDatasetIsOpenable() {
    using carta::zarr::internal::RequireOpenableDataset;

    // Nothing matched. Not an error about the file's contents: this library is not for it.
    {
        carta::zarr::internal::ProbeResult probe;
        probe.kind = ProbeKind::zarr_without_supported_schema;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::unsupported_schema,
                "an unmatched store should be refused as an unsupported schema");
    }

    // Something matched and was malformed. A different answer, and a consumer acts on it.
    {
        carta::zarr::internal::ProbeResult probe;
        probe.kind = ProbeKind::invalid_dataset;
        probe.diagnostics.push_back(
            carta::zarr::Diagnostic{carta::zarr::DiagnosticCode::invalid_metadata, "frequency is missing", "SKY"});
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::invalid_schema,
                "a malformed store should be refused as not conforming to its schema");
        Require(openable.error().message == "frequency is missing",
                "and it should carry what the probe worked out, not a generic refusal");
    }

    // A probe that had nothing to say falls back rather than reporting an empty message.
    {
        carta::zarr::internal::ProbeResult probe;
        probe.kind = ProbeKind::zarr_without_supported_schema;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && !openable.error().message.empty(), "a silent probe still needs a message");
    }

    // A profile matched a store with nothing openable in it. The profile is not refusing, so there
    // are no diagnostics; opening it would hand back a dataset a consumer can do nothing with.
    {
        carta::zarr::internal::ProbeResult probe;
        probe.kind = ProbeKind::supported_dataset;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::invalid_schema,
                "a matched store with no images should be refused");
    }

    // And the one case that opens.
    {
        carta::zarr::internal::ProbeResult probe;
        probe.kind = ProbeKind::supported_dataset;
        probe.images.push_back(carta::zarr::ImageEntry{});
        Require(static_cast<bool>(RequireOpenableDataset(probe, "/tmp/store")),
                "a matched store with an image should open");
    }
}

void TestDefaultImageSkipsUnopenablePreferredImage() {
    auto nodes = CompleteStore();
    nodes["SKY"] = SkyArray("complex64");
    nodes["RESIDUAL"] = SkyArray();

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the default-image store failed to open");
    auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "default-image discovery reported an error");
    Require(discovery.value().default_image_id == "RESIDUAL",
            "discovery selected an unopenable preferred image as the default");

    const auto sky = std::find_if(discovery.value().images.begin(), discovery.value().images.end(),
                                  [](const auto& image) { return image.id == "SKY"; });
    Require(sky != discovery.value().images.end() && !sky->openable &&
                HasDiagnostic(sky->diagnostics, carta::zarr::DiagnosticCode::unsupported_data_type),
            "the unopenable image did not carry its capability diagnostic");
}

// What consolidated metadata is for: the root's copy answers for every child, so discovery reads
// one node instead of one per variable. That saving is the whole reason the block exists -- on a
// store reached over a network each of those reads is a round trip -- so it is asserted as reads
// not taken, which is the only thing that tells a store that used the copy from one that ignored
// it. The children are present throughout: this is a cache, and the test would pass either way if
// it only checked the answer.
void TestConsolidatedMetadataDiscovery() {
    auto nodes = ConsolidatedStore();
    Require(nodes.size() > 1, "the consolidated store must carry its child nodes as well as the root copy");

    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::match, "a store using consolidated metadata was not matched");

    auto transport = MakeInMemoryTransport(nodes);
    auto store = carta::zarr::internal::OpenStore(transport);
    Require(static_cast<bool>(store), "the consolidated store failed to open");
    auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery over consolidated metadata reported an error");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "discovery did not find SKY through consolidated metadata");
    Require(transport->nodes_read() == std::set<std::string>{""},
            "discovery read a child node that the root's consolidated copy already answered for");
    Require(transport->listings() == 0,
            "discovery listed the hierarchy although consolidated metadata names every node");

    // The control: the same store without the root's copy reads every child instead, and reaches
    // the same answer. Both are supported; only one of them costs a read per variable.
    auto plain = MakeInMemoryTransport(CompleteStore());
    auto plain_store = carta::zarr::internal::OpenStore(plain);
    Require(static_cast<bool>(plain_store), "the unconsolidated store failed to open");
    auto plain_discovery = XradioProfile().Discover(plain_store.value());
    Require(static_cast<bool>(plain_discovery), "discovery without consolidated metadata reported an error");
    Require(OpenableImageIds(plain_discovery.value().images) == OpenableImageIds(discovery.value().images),
            "the two metadata layouts did not describe the same images");
    Require(plain->nodes_read().size() > 1, "a store without consolidated metadata has to read its children");
}

// A consolidated key is filed under the canonical name of the node it describes. A block is entitled
// to spell a node "./SKY" -- it names the node SKY -- and filed under that spelling the copy was
// never found by a read of "SKY": the node was read from the transport regardless, cached twice,
// and listed as "./SKY". That is an image id no reader asks for, and since openability is decided
// against the listing, SKY could not be opened by its own name.
void TestAConsolidatedKeyIsFiledUnderItsCanonicalName() {
    auto nodes = ConsolidatedStore();
    auto& root = nodes[""];
    const auto key = root.find("\"SKY\":");
    Require(key != std::string::npos, "the consolidated fixture has no SKY entry to respell");
    root.replace(key, 6, "\"./SKY\":");

    auto transport = MakeInMemoryTransport(nodes);
    auto store = carta::zarr::internal::OpenStore(transport);
    Require(static_cast<bool>(store), "a consolidated key spelled with ./ was refused");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed over a respelled consolidated key");
    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "the image was listed under the spelling the block used rather than under its name");
    Require(transport->nodes_read() == std::set<std::string>{""},
            "the node was read again although the consolidated copy described it");
}

// Two ways a consolidated block can be wrong about itself, and both refuse the store. For a store
// that consolidated its metadata the block is the whole listing, so a key passed over would be a
// variable missing from the dataset without a word said about it.
void TestAConsolidatedBlockThatMisnamesItsNodesIsRefused() {
    // A key that no node path can be.
    auto escaping = ConsolidatedStore();
    auto& escaping_root = escaping[""];
    escaping_root.replace(escaping_root.find("\"SKY\":"), 6, "\"../SKY\":");
    const auto outside = Open(escaping);
    Require(!outside && outside.error().code == ErrorCode::invalid_metadata,
            "a consolidated key naming a path outside the store was accepted");

    // One node under two spellings. Which document wins does not depend on the order the object is
    // walked in.
    auto doubled = ConsolidatedStore();
    auto& doubled_root = doubled[""];
    doubled_root.insert(doubled_root.find("\"SKY\":"), "\"./SKY\":" + SkyArray() + ",");
    const auto twice = Open(doubled);
    Require(!twice && twice.error().code == ErrorCode::invalid_metadata,
            "a consolidated block listing one node twice was accepted");
}

// The inventory decides once what each node is, so its callers cannot decide it three different
// ways.
void TestTheInventorySaysWhatEachNodeIs() {
    using carta::zarr::internal::NodeKind;
    auto nodes = CompleteStore();
    nodes["SUBDIR"] = RootGroup(false);
    // Says it is an array, and is not one: no shape.
    nodes["BROKEN"] = R"({"zarr_format":3,"node_type":"array","data_type":"float32"})";
    nodes["ODD"] = R"({"zarr_format":3,"node_type":"manifest"})";
    // Present and not a string. Asked for as a string with a default, it threw rather than answered,
    // and the whole discovery went with it.
    nodes["NUMBERED"] = R"({"zarr_format":3,"node_type":3})";

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the inventory store failed to open");
    const auto& inventory = store.value().Inventory();
    Require(static_cast<bool>(inventory), "the inventory could not be taken");

    const auto find = [&](const std::string& name) {
        const auto* entry = store.value().FindNode(name);
        Require(entry != nullptr, "the inventory has no entry for " + name);
        return entry;
    };

    const auto* sky = find("SKY");
    Require(sky->kind == NodeKind::array && sky->array != nullptr && static_cast<bool>(*sky->array) && !sky->reason,
            "an array was not carried with its parsed metadata");
    const auto* group = find("SUBDIR");
    Require(group->kind == NodeKind::group && group->array == nullptr && !group->reason,
            "a group was parsed as an array, or given a reason to be something else");
    const auto* broken = find("BROKEN");
    Require(broken->kind == NodeKind::array && broken->array != nullptr && !*broken->array,
            "an array whose metadata does not parse was not carried with its refusal");
    for (const std::string name : {"ODD", "NUMBERED"}) {
        const auto* entry = find(name);
        Require(entry->kind == NodeKind::unrecognised && entry->array == nullptr && entry->reason.has_value(),
                "a node that does not say it is a group or an array was not unrecognised: " + name);
    }

    // Sorted by name, which is the order discovery's diagnostics come out in.
    Require(std::is_sorted(inventory.value().begin(), inventory.value().end(),
                           [](const auto& left, const auto& right) { return left.name < right.name; }),
            "the inventory was not sorted by name");
    // A name is looked up as it was handed out. Another spelling of it is a name nobody was given.
    Require(store.value().FindNode("./SKY") == nullptr, "a respelled name was looked up as the node");
    Require(store.value().FindNode("NOPE") == nullptr, "a node the store does not hold was found");

    Require(static_cast<bool>(XradioProfile().Discover(store.value())),
            "discovery could not walk a store holding a node whose node_type is not a string");
}

// A node whose document will not parse is diagnosed, not refused. GLOSSARY.md says so of a node
// whose metadata would not parse -- the rest of the dataset is still readable -- so a stray broken
// document beside a perfectly good image does not close the dataset.
void TestANodeThatWillNotParseIsDiagnosedNotRefused() {
    auto nodes = CompleteStore();
    nodes["JUNK"] = "{not valid json";
    nodes["ODD"] = R"({"zarr_format":3,"node_type":"manifest"})";

    Require(Probe(nodes).kind == SchemaMatchKind::match,
            "a document that will not parse refused a store whose image is fine");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the store holding an unparseable node failed to open");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery refused a store holding an unparseable node");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "the image beside an unparseable node was not openable");

    // Neither said what it was, and both say so under the one code for that. Neither is called an
    // array, because neither claimed to be one.
    const auto& said = discovery.value().diagnostics;
    const auto said_of = [&](const std::string& node, carta::zarr::DiagnosticCode code) {
        return std::any_of(said.begin(), said.end(), [&](const auto& diagnostic) {
            return diagnostic.node_path == node && diagnostic.code == code;
        });
    };
    Require(said_of("JUNK", carta::zarr::DiagnosticCode::unrecognised_node),
            "a document that will not parse was not diagnosed as such");
    Require(said_of("ODD", carta::zarr::DiagnosticCode::unrecognised_node),
            "a node_type Zarr does not define was not diagnosed as such");
    Require(!HasDiagnostic(said, carta::zarr::DiagnosticCode::unreadable_array),
            "a node that never said it was an array was called one");

    // The declared size counts the arrays it can see, and a node it cannot read is not one of them.
    const auto size = carta::zarr::internal::DatasetSizeBytes(
        store.value(), std::chrono::steady_clock::now() + std::chrono::seconds(5));
    Require(size && size.value().bytes == 592, "an unparseable node left the dataset without a declared size");
}

// The dataset-level half of a codec that will not parse: the array is diagnosed, and the dataset
// beside it is not. A numeric codec name does not throw past every per-node diagnostic and close
// the store.
void TestAnArrayWithAMalformedCodecIsDiagnosedNotRefused() {
    auto nodes = CompleteStore();
    nodes["JUNK"] = R"({"shape":[2],"data_type":"float32",)"
                    R"("chunk_grid":{"name":"regular","configuration":{"chunk_shape":[2]}},)"
                    R"("codecs":[{"name":"bytes"},{"name":17}],"zarr_format":3,"node_type":"array"})";

    Require(Probe(nodes).kind == SchemaMatchKind::match, "a malformed codec refused a store whose image is fine");
    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the store holding a malformed codec failed to open");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery refused a store holding a malformed codec");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "the image beside a malformed codec was not openable");
    const auto& said = discovery.value().diagnostics;
    Require(std::any_of(said.begin(), said.end(),
                        [](const auto& diagnostic) {
                            return diagnostic.node_path == "JUNK" &&
                                   diagnostic.code == carta::zarr::DiagnosticCode::unreadable_array;
                        }),
            "the array with a malformed codec was not diagnosed as unreadable");
}

// The report latches: once a requirement is unmet, later ones are no-ops. A store with two faults
// is therefore diagnosed once, by the first fault reached -- the behaviour a probe had when every
// check returned early, now stated somewhere rather than emerging from the control flow.
//
// Both faults have to be ones the probe still reaches. A coordinate of the wrong length is not:
// agreeing with it is part of being an image this profile opens, so a store with one has no
// openable image and never gets as far as the requirements below.
void TestFirstFaultIsTheOnlyDiagnostic() {
    auto nodes = CompleteStore();
    // Right length, names another axis.
    nodes["frequency"] = NumericArray("[3]", R"(["x"])", "float64", R"({"units":"Hz","frame":"LSRK"})");
    nodes["polarization"] = PolarizationArray("\"float64\"");  // wrong data type

    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a store with two faults was not reported as invalid");
    Require(probe.diagnostics.size() == 1, "the report diagnosed more than the first unmet requirement");
    Require(probe.diagnostics.front().node_path == "frequency",
            "the diagnostic did not come from the first requirement checked");
    Require(!HasDiagnostic(probe, carta::zarr::DiagnosticCode::unsupported_data_type),
            "a requirement after the first failure was still checked");
}

// Store-level rejections, reported before any profile is consulted.
// Choosing a flag is the half of it that consults the store: what a flag has to be once it is found
// takes no store at all and is checked in tests/flag_test.cc.
//
// The same root, declaring `data_groups` as given.
std::string RootGroupWithDataGroups(const std::string& data_groups) {
    auto root = RootGroup();
    root.insert(root.find("},\"zarr_format\""), ",\"data_groups\":" + data_groups);
    return root;
}

// A declared flag is the image's own statement that its pixels need a mask, so an image naming one
// that cannot be read is closed rather than opened unmasked. Reads apply the mask by default, and an
// unusable mask reported as no mask would show flagged pixels as valid -- the one failure a consumer
// has no way to notice.
void TestADeclaredFlagIsBinding() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    const auto with_declared_flag = [&](bool write_the_flag) {
        auto nodes = CompleteStore();
        nodes["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky","flag":"MASK_0"})");
        if (write_the_flag) {
            nodes["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
        }
        return nodes;
    };

    auto present = Open(with_declared_flag(true));
    Require(static_cast<bool>(present), "the declared-flag store did not open");
    const auto& present_image = present.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(present_image), "SKY was not readable in the declared-flag store");
    auto declared = DetermineFlag(present.value(), present_image.value(), "SKY");
    Require(static_cast<bool>(declared), "a well-formed declared flag was refused");
    Require(declared.value() == "MASK_0", "the declared flag was not the one selected");

    auto absent = Open(with_declared_flag(false));
    Require(static_cast<bool>(absent), "the missing-flag store did not open");
    const auto& absent_image = absent.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(absent_image), "SKY was not readable in the missing-flag store");
    auto missing = DetermineFlag(absent.value(), absent_image.value(), "SKY");
    Require(!missing, "an image declaring a flag variable that does not exist was opened");
}

// The listing's half of the same rule. An image whose declared flag cannot mask it is closed when it
// is opened, so it has to be closed when it is listed too: it was listed openable and chosen as the
// default, and a consumer offered an image it could not open.
//
// A flag is declared on the image's own attribute or by a data group naming the image as its sky, and
// the listing asks both, as DetermineFlag does: it asked only the attribute, so a data group declaring
// a missing flag left the image listed openable and refused when opened.
void TestADeclaredFlagThatCannotMaskClosesTheImageInTheListing() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    const auto listing = [&](bool by_group, const std::string& flag_node, const std::string& data_groups) {
        auto nodes = CompleteStore();
        if (by_group) {
            nodes[""] = RootGroupWithDataGroups(data_groups);
        } else {
            nodes["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky","flag":"MASK_0"})");
        }
        nodes["MODEL"] = SkyArray();
        if (!flag_node.empty()) {
            nodes["MASK_0"] = flag_node;
        }
        auto store = Open(nodes);
        Require(static_cast<bool>(store), "the declared-flag listing store did not open");
        const auto profile = XradioProfile();
        auto discovery = profile.Discover(store.value());
        Require(static_cast<bool>(discovery), "discovery failed on the declared-flag listing store");
        return std::make_pair(discovery.value(), profile.Describe(store.value(), "SKY"));
    };
    const std::string declares_mask = R"({"base":{"sky":"SKY","flag":"MASK_0"}})";
    const auto usable_flag = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");

    for (const bool by_group : {false, true}) {
        const std::string how = by_group ? " declared by a data group" : " declared on the image";
        const auto usable = listing(by_group, usable_flag, declares_mask);
        Require(OpenableImageIds(usable.first.images) == std::vector<std::string>{"SKY", "MODEL"},
                "an image whose flag" + how + " can mask it was not listed openable");

        std::vector<std::tuple<std::string, std::string, std::string>> unusable{
            {"", declares_mask, "a flag that does not exist"},
            {NumericArray("[1,3,2,4,4]", sky_dimensions, "bool", R"({"type":"flag"})"), declares_mask,
             "a flag of another shape"},
            {NumericArray("[1,3,2,4,5]", sky_dimensions, "uint8", R"({"type":"flag"})"), declares_mask,
             "a flag that is not boolean"},
        };
        if (by_group) {
            unusable.emplace_back(usable_flag,
                                  R"({"base":{"sky":"SKY","flag":"MASK_0"},"robust":{"sky":"SKY","flag":"MASK_1"}})",
                                  "two different flags");
        }
        for (const auto& [flag_node, data_groups, what] : unusable) {
            const auto [discovery, sky] = listing(by_group, flag_node, data_groups);
            const auto with = "an image with " + what + how;
            Require(ImageIds(discovery.images) == std::vector<std::string>{"SKY", "MODEL"},
                    with + " was dropped from the listing rather than listed with its reason");
            Require(OpenableImageIds(discovery.images) == std::vector<std::string>{"MODEL"},
                    with + " was listed openable");
            Require(discovery.default_image_id == "MODEL", with + " was chosen as the default");
            Require(HasDiagnostic(discovery.images.front().diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
                    with + " was listed without saying why it will not open");
            Require(!sky && sky.error().code == ErrorCode::invalid_metadata,
                    "describing " + with + " was not refused as invalid metadata");
        }
    }
}

// An image is described by the five axes this profile knows. One with a sixth is not listed
// openable, since every read of it would fail on the rank it had not been told.
void TestAnImageWithAnAxisBeyondTheFiveIsNotOpenable() {
    auto nodes = CompleteStore();
    nodes["MODEL"] = NumericArray("[1,3,2,4,5,1]", R"(["time","frequency","polarization","l","m","extra"])", "float32",
                                  R"({"units":"Jy/beam","type":"sky"})");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the extra-axis store failed to open");
    const auto profile = XradioProfile();
    const auto discovery = profile.Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed on the extra-axis store");
    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY", "MODEL"},
            "the extra-axis image was dropped from the listing rather than listed with its reason");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
            "an image with an axis this profile cannot describe was listed openable");
    Require(HasDiagnostic(discovery.value().images.back().diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
            "the extra-axis image was listed without saying why it will not open");

    const auto model = profile.Describe(store.value(), "MODEL");
    Require(!model && model.error().code == ErrorCode::invalid_metadata,
            "an image with an axis this profile cannot describe was not refused as invalid metadata");
}

carta::zarr::Result<std::string> FlagOfSky(const std::map<std::string, std::string>& nodes) {
    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the data-group store did not open");
    const auto& image = store.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(image), "SKY was not readable in the data-group store");
    return DetermineFlag(store.value(), image.value(), "SKY");
}

// XRADIO's writer ties a flag to its sky image in the root's data_groups and nowhere else: the
// image names no flag. That link is the declaration, and a sky image that has none has no flag,
// however well another fits it.
void TestADataGroupDeclaresItsSkysFlag() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    auto nodes = CompleteStore();
    nodes["FLAG_1"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    nodes["FLAG_2"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");

    nodes[""] =
        RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":"FLAG_2"},"other":{"sky":"OTHER","flag":"FLAG_1"}})");
    auto declared = FlagOfSky(nodes);
    Require(declared && declared.value() == "FLAG_2", "the flag data_groups names for SKY was not the one chosen");

    // A sky image whose group names no flag has none: FLAG_2 fits it and is no one's, and is not taken.
    nodes[""] = RootGroupWithDataGroups(R"({"base":{"sky":"SKY"},"other":{"sky":"OTHER","flag":"FLAG_1"}})");
    auto undeclared = FlagOfSky(nodes);
    Require(undeclared && undeclared.value().empty(), "a sky image that declares no flag was given one");

    // The image's own attribute is its own statement, and outranks a group's.
    auto own = nodes;
    own["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky","flag":"FLAG_1"})");
    own[""] = RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":"FLAG_2"}})");
    auto attribute = FlagOfSky(own);
    Require(attribute && attribute.value() == "FLAG_1", "a data group outranked the image's own flag attribute");

    // Declared flags are binding: one that cannot serve closes the image rather than leaving it unmasked.
    nodes["FLAG_3"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "float32", R"({"type":"flag"})");
    nodes[""] = RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":"FLAG_3"}})");
    auto unusable = FlagOfSky(nodes);
    Require(!unusable && unusable.error().code == ErrorCode::unsupported_data_type,
            "an unusable flag data_groups declares did not close the image");

    // Two groups of the same sky naming different flags contradict each other.
    nodes[""] =
        RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":"FLAG_1"},"robust":{"sky":"SKY","flag":"FLAG_2"}})");
    auto contradictory = FlagOfSky(nodes);
    Require(!contradictory && contradictory.error().code == ErrorCode::invalid_metadata,
            "two different flags declared for one image were not refused");

    // The same flag named by both is one declaration.
    nodes[""] =
        RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":"FLAG_1"},"robust":{"sky":"SKY","flag":"FLAG_1"}})");
    auto agreeing = FlagOfSky(nodes);
    Require(agreeing && agreeing.value() == "FLAG_1", "two groups naming the same flag were refused");

    // "./FLAG_1" names the node FLAG_1, and every rule above holds of it however it is spelt.
    const auto spelt = [&](const std::string& data_groups) {
        nodes[""] = RootGroupWithDataGroups(data_groups);
        return FlagOfSky(nodes);
    };
    auto dotted = spelt(R"({"base":{"sky":"SKY","flag":"./FLAG_2"}})");
    Require(dotted && dotted.value() == "FLAG_2", "a declared ./FLAG_2 was not chosen as the node FLAG_2");
    auto aliases = spelt(R"({"base":{"sky":"SKY","flag":"FLAG_1"},"robust":{"sky":"SKY","flag":"./FLAG_1"}})");
    Require(aliases && aliases.value() == "FLAG_1", "two spellings of one flag were taken for two flags");
    auto dotted_sky = spelt(R"({"base":{"sky":"./SKY","flag":"FLAG_2"}})");
    Require(dotted_sky && dotted_sky.value() == "FLAG_2", "a group whose sky is spelt ./SKY was not SKY's");
    auto own_dotted = nodes;
    own_dotted["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky","flag":"./FLAG_1"})");
    own_dotted[""] = RootGroup();
    auto attribute_dotted = FlagOfSky(own_dotted);
    Require(attribute_dotted && attribute_dotted.value() == "FLAG_1",
            "an image's own ./FLAG_1 was not named as the node it is");
}

// With nothing declared, an image that is not a sky image takes the flag XRADIO names FLAG_<image>,
// and nothing else: not a flag that merely has its shape, however many or few there are. XRADIO
// cannot declare such an image's flag (casangi/xradio#659), and names it so.
void TestAnUndeclaredFlagIsTakenOnlyByItsName() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    auto nodes = CompleteStore();
    nodes["POINT_SPREAD_FUNCTION"] = SkyArray("float32", R"({"type":"point_spread_function"})");
    nodes["FLAG_POINT_SPREAD_FUNCTION"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    nodes["FLAG_OTHER"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    const auto flag_of = [](const std::map<std::string, std::string>& store_nodes, const std::string& image_id) {
        auto store = Open(store_nodes);
        Require(static_cast<bool>(store), "the undeclared-flag store did not open");
        const auto& image = store.value().ReadArrayMetadata(image_id);
        Require(static_cast<bool>(image), image_id + " was not readable in the undeclared-flag store");
        return DetermineFlag(store.value(), image.value(), image_id);
    };

    auto named = flag_of(nodes, "POINT_SPREAD_FUNCTION");
    Require(named && named.value() == "FLAG_POINT_SPREAD_FUNCTION",
            "the flag XRADIO names for a point spread function was not taken, beside another of its shape");

    // The sky declares none, so it has none, though both fit it.
    auto sky = flag_of(nodes, "SKY");
    Require(sky && sky.value().empty(), "a sky image that declares no flag took one that fits it");

    // Without its own, the point spread function has none: a flag of its shape is not its flag.
    auto without = nodes;
    without.erase("FLAG_POINT_SPREAD_FUNCTION");
    auto shaped = flag_of(without, "POINT_SPREAD_FUNCTION");
    Require(shaped && shaped.value().empty(), "a flag that only has the image's shape was taken as its flag");

    // Its own, of another shape, cannot mask it, and is not taken.
    auto other_shape = nodes;
    other_shape["FLAG_POINT_SPREAD_FUNCTION"] =
        NumericArray("[1,3,2,4,4]", sky_dimensions, "bool", R"({"type":"flag"})");
    auto mismatched = flag_of(other_shape, "POINT_SPREAD_FUNCTION");
    Require(mismatched && mismatched.value().empty(), "a FLAG_<image> of another shape was taken as its flag");

    // The image's own flag attribute is its declaration, and outranks the name.
    auto own = nodes;
    own["POINT_SPREAD_FUNCTION"] = SkyArray("float32", R"({"type":"point_spread_function","flag":"FLAG_OTHER"})");
    auto declared = flag_of(own, "POINT_SPREAD_FUNCTION");
    Require(declared && declared.value() == "FLAG_OTHER",
            "a point spread function's own flag attribute was outranked by FLAG_<image>");

    // Nor one that is not typed a flag, nor one that is not boolean.
    for (const auto& [replacement, what] : std::vector<std::pair<std::string, std::string>>{
             {NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", "{}"), "an untyped FLAG_<image>"},
             {NumericArray("[1,3,2,4,5]", sky_dimensions, "uint8", R"({"type":"flag"})"), "a FLAG_<image> of uint8"}}) {
        auto replaced = nodes;
        replaced["FLAG_POINT_SPREAD_FUNCTION"] = replacement;
        auto refused = flag_of(replaced, "POINT_SPREAD_FUNCTION");
        Require(refused && refused.value().empty(), what + " was taken as the image's flag");
    }

    // Nor one a data group gives a sky image as its flag, however the group spells it.
    for (const std::string& spelt : {"FLAG_POINT_SPREAD_FUNCTION", "./FLAG_POINT_SPREAD_FUNCTION"}) {
        auto owned = nodes;
        owned[""] = RootGroupWithDataGroups(R"({"base":{"sky":"SKY","flag":")" + spelt + R"("}})");
        auto taken = flag_of(owned, "POINT_SPREAD_FUNCTION");
        Require(taken && taken.value().empty(),
                "a flag a data group gives the sky, spelt " + spelt + ", was taken as the point spread function's");
    }
}

// A root that says its consolidated metadata is null has none, which is how zarr-python reads it:
// it opens the hierarchy by listing, as a root without the member would, rather than refusing a
// malformed block and closing a store whose every node is fine.
void TestANullConsolidatedBlockIsNoConsolidation() {
    auto nodes = CompleteStore();
    auto root = RootGroup();
    root.insert(root.size() - 1, ",\"consolidated_metadata\":null");
    nodes[""] = root;

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "a root with null consolidated metadata was refused");
    Require(Probe(nodes).kind == SchemaMatchKind::match, "a store with null consolidated metadata was not matched");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(
        static_cast<bool>(discovery) && OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
        "the image of a store with null consolidated metadata was not found by listing");
}

void TestStoreRejections() {
    const auto no_root = Open({{"SKY", SkyArray()}});
    Require(!no_root && no_root.error().code == ErrorCode::not_zarr,
            "a transport with no root node was not rejected as not_zarr");

    const auto version = Open({{"", R"({"zarr_format": 2, "node_type": "group"})"}});
    Require(!version && version.error().code == ErrorCode::unsupported_zarr_version, "Zarr format 2 was not rejected");

    const auto array_root = Open({{"", R"({"zarr_format": 3, "node_type": "array"})"}});
    Require(!array_root && array_root.error().code == ErrorCode::not_zarr, "an array root was not rejected");

    const auto malformed = Open({{"", "{not valid json"}});
    Require(!malformed && malformed.error().code == ErrorCode::invalid_metadata,
            "malformed root metadata was not reported as invalid");

    auto store = Open(CompleteStore());
    Require(static_cast<bool>(store), "the complete store failed to open");
    const auto unknown = carta::zarr::internal::SchemaProfile::For("future.schema");
    Require(!unknown && unknown.error().code == ErrorCode::unsupported_schema, "an unknown schema id was accepted");
}

// How large a dataset is, when the transport cannot say how much room it takes.
//
// An in-memory transport gives up measuring honestly -- it holds no bytes -- so the fallback to the
// declared size is reachable with a map, rather than by pointing at a real store and setting the
// timeout to zero.
void TestSizeFallsBackWhenTheStoreCannotBeMeasured() {
    auto store = Open(CompleteStore());
    Require(static_cast<bool>(store), "OpenStore rejected the in-memory store");

    const auto size = carta::zarr::internal::DatasetSizeBytes(
        store.value(), std::chrono::steady_clock::now() + std::chrono::seconds(5));
    Require(static_cast<bool>(size), "sizing an in-memory store failed");
    Require(size.value().basis == carta::zarr::SizeBasis::declared,
            "a size the transport could not measure must be reported as the declared one");
    // SKY is 120 float32 at 480 bytes; time, frequency, l and m are 1, 3, 4 and 5 float64 at 8, 24,
    // 32 and 40; polarization is two four-byte labels at 8. The arrays, not the store.
    Require(size.value().bytes == 592, "the logical total was " + std::to_string(size.value().bytes) + ", not 592");

    // A deadline that has already passed reaches the same answer by the same route: this transport
    // refuses whatever the clock says, and every way of failing to measure means the declared size.
    const auto expired = carta::zarr::internal::DatasetSizeBytes(
        store.value(), std::chrono::steady_clock::now() - std::chrono::seconds(1));
    Require(expired && expired.value().basis == carta::zarr::SizeBasis::declared && expired.value().bytes == 592,
            "an expired deadline did not reach the same declared size");
}

// And the fallback does not turn every failure into a number: a store with nothing to add up is
// still an error, not a zero.
void TestSizeRefusesAStoreWithNoArrays() {
    auto store = Open({{"", RootGroup()}});
    Require(static_cast<bool>(store), "OpenStore rejected a bare root group");

    const auto size = carta::zarr::internal::DatasetSizeBytes(store.value(), std::chrono::steady_clock::now());
    Require(!size && size.error().code == ErrorCode::invalid_metadata, "a store holding no arrays was given a size");
}

// Every rule that turns metadata and coordinate values into a descriptor, reached with the values
// handed in. This build reads no values, so while describing an image read its own, everything
// below was reachable only from a store written to disk -- and the image_type fallback and the time
// axis's scale and format were not reached at all.
void TestADescriptionIsBuiltFromTheValuesItIsGiven() {
    using carta::zarr::internal::xradio::CoordinateValues;
    using carta::zarr::internal::xradio::DescribeImageFrom;

    auto nodes = CompleteStore();
    nodes["frequency"] = FrequencyArray("[3]", R"({"units":"Hz","frame":"LSRK","reference_frequency":{"data":1.402e9},
                                  "rest_frequency":{"data":1.420405751e9}})");
    nodes["time"] = NumericArray("[1]", R"(["time"])", "float64", R"({"units":"s","scale":"utc","format":"unix"})");
    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the in-memory store failed to open");

    CoordinateValues values;
    values.l = {-0.001, 0.0, 0.001, 0.002};
    values.m = {-0.002, -0.001, 0.0, 0.001, 0.002};
    values.frequency = {1.4e9, 1.401e9, 1.402e9};
    values.time = {1.6e9};
    values.polarization = {"I", "Q"};

    const auto described = DescribeImageFrom(store.value(), "SKY", values);
    Require(static_cast<bool>(described),
            "an image with its values in hand should be described: " + (described ? "" : described.error().message));
    const auto& descriptor = described.value().descriptor;

    Require(descriptor.direction && descriptor.direction->projection == "SIN", "the direction comes from the root");
    Require(descriptor.spectral.has_value(), "three channels make a spectral coordinate");
    const auto& spectral = *descriptor.spectral;
    Require(spectral.values == values.frequency, "the channels are the values handed in");
    Require(spectral.reference_frame == "LSRK", "the frame is frequency's own");
    Require(spectral.reference_value && *spectral.reference_value == 1.402e9,
            "the reference is the frequency named, not the first channel");
    Require(spectral.reference_pixel && *spectral.reference_pixel == 3.0,
            "and its pixel is that channel's, counted from one");
    Require(spectral.rest_frequency && *spectral.rest_frequency == 1.420405751e9, "the rest frequency is read");

    Require(descriptor.polarization && descriptor.polarization->labels == values.polarization,
            "the polarization labels are the ones handed in");
    Require(descriptor.temporal.has_value(), "a time value makes a temporal coordinate");
    Require(descriptor.temporal->values == values.time, "the times are the values handed in");
    Require(descriptor.temporal->scale == "UTC" && descriptor.temporal->format == "UNIX",
            "the time's scale and format come from its attributes, uppercased");
}

// A numeric array spelled out field by field, for the cases where the array's own document has to
// differ from the root's copy in one of them and nothing else.
std::string ArrayDocument(const std::string& shape, const std::string& chunks, const std::string& dimensions,
                          const std::string& data_type = "float64", const std::string& codecs = "",
                          const std::string& fill_value = "") {
    return "{\"shape\":" + shape + ",\"data_type\":\"" + data_type +
           "\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":" + chunks +
           "}},\"attributes\":{},\"dimension_names\":" + dimensions + (codecs.empty() ? "" : ",\"codecs\":" + codecs) +
           (fill_value.empty() ? "" : ",\"fill_value\":" + fill_value) + ",\"zarr_format\":3,\"node_type\":\"array\"}";
}

void RequireRefusedAs(const carta::zarr::Error& error, ErrorCode code, const std::string& node,
                      const std::string& what) {
    Require(error.code == code, what + " was refused as " + carta::zarr::ErrorCodeName(error.code) + " rather than " +
                                    carta::zarr::ErrorCodeName(code) + ": " + error.message);
    Require(error.node_path == node, what + " was refused naming '" + error.node_path + "' rather than '" + node + "'");
}

// The store decides what an image is from the root's copy of each array's document, and every value
// is decoded from the array's own: TensorStore opens an array from its own document, and so does
// the label decoder. Where a rewrite left the two disagreeing about what the library derives from
// the copy -- extent, dimension names and order, data type, how it is chunked, what its attributes
// say its values mean -- the array is not the one the store described, and every read of it says so
// before it asks where the array lives. Here that is before the in-memory transport says it holds no
// array data, which is the only other answer a read gets from it.
void TestAnArrayDisagreeingWithTheRootsCopyIsRefused() {
    const std::vector<std::pair<std::string, std::string>> rewritten{
        {"another extent", ArrayDocument("[6]", "[4]", R"(["l"])")},
        {"another dimension name", ArrayDocument("[4]", "[4]", R"(["x"])")},
        {"another data type", ArrayDocument("[4]", "[4]", R"(["l"])", "float32")},
        {"another chunk", ArrayDocument("[4]", "[2]", R"(["l"])")},
        {"shards", ArrayDocument("[4]", "[4]", R"(["l"])", "float64",
                                 R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[2],
                                      "codecs":[{"name":"bytes","configuration":{"endian":"little"}}],
                                      "index_codecs":[{"name":"bytes","configuration":{"endian":"little"}}]}}])")},
        // The same chunks, each now read out of a shard the copy does not mention: the chunk shape
        // agrees and the shard shape a read is planned around does not.
        {"shards around the same chunks",
         ArrayDocument("[4]", "[4]", R"(["l"])", "float64",
                       R"([{"name":"sharding_indexed","configuration":{"chunk_shape":[4],
                                                             "codecs":[{"name":"bytes","configuration":{"endian":"little"}}],
                                                             "index_codecs":[{"name":"bytes","configuration":{"endian":"little"}}]}}])")},
        // What the values mean: a unit the copy does not give them.
        {"other attributes",
         [] {
             auto document = ArrayDocument("[4]", "[4]", R"(["l"])");
             return document.replace(document.find("\"attributes\":{}"), 15, R"("attributes":{"units":"deg"})");
         }()},
        {"a document that will not parse", "{"},
        {"a group", R"({"zarr_format":3,"node_type":"group"})"},
    };
    for (const auto& [how, own] : rewritten) {
        auto nodes = ConsolidatedStore();
        nodes["l"] = own;
        auto store = Open(nodes);
        Require(static_cast<bool>(store), "the store whose l was rewritten with " + how + " did not open");
        const auto read = store.value().ReadNumericArray("l");
        Require(!read, "l rewritten with " + how + " was read");
        RequireRefusedAs(read.error(), ErrorCode::invalid_metadata, "l", "l rewritten with " + how);
    }

    // The copy names an array whose own document is not there: a store that only a reader honouring
    // the copy alone could open, which ADR 0004 calls malformed.
    auto missing = ConsolidatedStore();
    missing.erase("l");
    auto missing_store = Open(missing);
    Require(static_cast<bool>(missing_store), "the store missing l's own document did not open");
    const auto missing_read = missing_store.value().ReadNumericArray("l");
    Require(!missing_read, "an array whose own document is missing was read");
    RequireRefusedAs(missing_read.error(), ErrorCode::invalid_metadata, "l", "an array whose own document is missing");

    // Labels are decoded without TensorStore, and are held to the same rule.
    auto labels = ConsolidatedStore();
    labels["polarization"] =
        R"({"shape":[2],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":4}},
        "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[1]}},"attributes":{},
        "dimension_names":["polarization"],"zarr_format":3,"node_type":"array"})";
    auto labels_store = Open(labels);
    Require(static_cast<bool>(labels_store), "the store whose labels were rechunked did not open");
    const auto label_read = labels_store.value().ReadStringArray1D("polarization");
    Require(!label_read, "labels chunked differently from the root's copy were read");
    RequireRefusedAs(label_read.error(), ErrorCode::invalid_metadata, "polarization",
                     "labels chunked differently from the root's copy");
}

// What the copy is not held to. A reader decodes with the array's own document, so how its chunks
// are compressed and what an unwritten one reads as are the own document's business, and a copy
// spelling them otherwise -- a later zarr-python writing a default out, a level changed -- is not a
// different array. Passing the check, the read goes on to ask where the array lives.
void TestAnArrayIsNotHeldToHowItsChunksAreEncoded() {
    auto nodes = ConsolidatedStore();
    nodes["l"] = ArrayDocument("[4]", "[4]", R"(["l"])", "float64",
                               R"([{"name":"bytes","configuration":{"endian":"little"}},
                                   {"name":"zstd","configuration":{"level":3,"checksum":false}}])",
                               "\"NaN\"");
    nodes["polarization"] =
        R"({"shape":[2],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":4}},
        "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[2]}},"attributes":{},"fill_value":"V",
        "codecs":[{"name":"bytes","configuration":{"endian":"little"}},{"name":"zstd","configuration":{"level":1}}],
        "dimension_names":["polarization"],"zarr_format":3,"node_type":"array"})";
    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the store whose encodings were respelled did not open");
    const auto values = store.value().ReadNumericArray("l");
    Require(!values && values.error().code == ErrorCode::unsupported_transport,
            "an array encoded differently from the root's copy did not pass the check: " +
                (values ? std::string("it was read") : values.error().message));
    const auto labels = store.value().ReadStringArray1D("polarization");
    Require(!labels && labels.error().code == ErrorCode::unsupported_transport,
            "labels encoded differently from the root's copy did not pass the check: " +
                (labels ? std::string("they were read") : labels.error().message));
}

// Checking an array costs reading its own document, and only where the store had not read it: the
// root's copy saved that read, and an array whose metadata came from its own document has nothing to
// be checked against. Either way it is read once for the life of the store, and probing reads none.
void TestCheckingAnArrayReadsItsOwnDocumentOnce() {
    {
        const auto transport = MakeInMemoryTransport(ConsolidatedStore());
        auto store = carta::zarr::internal::OpenStore(transport);
        Require(static_cast<bool>(store), "the consolidated store did not open");
        Require(static_cast<bool>(XradioProfile().Probe(store.value())), "the consolidated store did not probe");
        Require(transport->reads("l") == 0, "probing a consolidated store read an array's own document");
        (void)store.value().ReadNumericArray("l");
        (void)store.value().ReadNumericArray("l");
        Require(transport->reads("l") == 1, "a consolidated array's own document was read " +
                                                std::to_string(transport->reads("l")) + " times, not once");
    }
    {
        const auto transport = MakeInMemoryTransport(CompleteStore());
        auto store = carta::zarr::internal::OpenStore(transport);
        Require(static_cast<bool>(store), "the plain store did not open");
        (void)store.value().ReadNumericArray("l");
        (void)store.value().ReadNumericArray("l");
        Require(transport->reads("l") == 1, "an array the store read for itself was read " +
                                                std::to_string(transport->reads("l")) + " times, not once");
    }
}

// Describing an image holds the image's own array to the store's copy before any value is read, so
// an image that opens can be read: one whose own document disagrees does not open only to fail
// every read. Its flag is held to it as soon as the flag is known.
void TestDescribingAnImageHoldsItsArraysToTheirOwnDocuments() {
    using carta::zarr::internal::xradio::CoordinateValues;
    using carta::zarr::internal::xradio::DescribeImage;
    using carta::zarr::internal::xradio::DescribeImageFrom;
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";

    auto sky = ConsolidatedStore();
    sky["SKY"] = SkyArray("float64");
    auto sky_store = Open(sky);
    Require(static_cast<bool>(sky_store), "the store whose SKY was rewritten did not open");
    const auto described = DescribeImage(sky_store.value(), "SKY");
    Require(!described, "an image whose own document disagrees with the root's copy was described");
    RequireRefusedAs(described.error(), ErrorCode::invalid_metadata, "SKY",
                     "an image whose own document disagrees with the root's copy");

    auto flagged = CompleteStore();
    flagged["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky","flag":"MASK_0"})");
    flagged["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    auto flagged_nodes = ConsolidatedStore(flagged);
    flagged_nodes["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "uint8", R"({"type":"flag"})");
    auto flag_store = Open(flagged_nodes);
    Require(static_cast<bool>(flag_store), "the store whose flag was rewritten did not open");
    const auto with_flag = DescribeImageFrom(flag_store.value(), "SKY", CoordinateValues{});
    Require(!with_flag, "an image whose flag's own document disagrees with the root's copy was described");
    RequireRefusedAs(with_flag.error(), ErrorCode::invalid_metadata, "MASK_0",
                     "a flag whose own document disagrees with the root's copy");
}

// The pixels and the flag are read through the same check as every other array, whatever was
// described before: a read is the last place an array the store did not describe could get through.
void TestAPixelReadHoldsItsArrayToItsOwnDocument() {
    using carta::zarr::internal::xradio::CoordinateValues;
    using carta::zarr::internal::xradio::DescribeImageFrom;
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";

    auto agreeing = CompleteStore();
    agreeing["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    auto agreeing_store = Open(agreeing);
    Require(static_cast<bool>(agreeing_store), "the agreeing store did not open");
    const auto described = DescribeImageFrom(agreeing_store.value(), "SKY", CoordinateValues{});
    Require(static_cast<bool>(described), "the agreeing store's SKY was not described");
    carta::zarr::ReadRequest request;
    for (std::size_t axis = 0; axis < described.value().descriptor.axes.size(); ++axis) {
        request.axes.push_back(carta::zarr::Range{0, 1, 1});
    }
    const auto selection = carta::zarr::internal::zarr::BuildSelection(
        described.value().descriptor, request, carta::zarr::internal::zarr::DestinationOrder::logical);
    Require(static_cast<bool>(selection), "a one-pixel selection was refused");

    auto rewritten = ConsolidatedStore(agreeing);
    rewritten["SKY"] = SkyArray("float64");
    rewritten["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "uint8", R"({"type":"flag"})");
    auto store = Open(rewritten);
    Require(static_cast<bool>(store), "the rewritten store did not open");
    float pixel = 0.0F;
    const auto pixels = store.value().ReadPixelsInto<float>("SKY", selection.value(), {&pixel, 1}, {});
    Require(!pixels, "pixels of an array disagreeing with the root's copy were read");
    RequireRefusedAs(pixels.error(), ErrorCode::invalid_metadata, "SKY", "pixels of a rewritten array");
    std::uint8_t good = 0;
    const auto mask = store.value().ReadPixelsInto<std::uint8_t>("MASK_0", selection.value(), {&good, 1}, {});
    Require(!mask, "a flag disagreeing with the root's copy was read");
    RequireRefusedAs(mask.error(), ErrorCode::invalid_metadata, "MASK_0", "a rewritten flag");
}

// The coordinates are said once, and every list the profile keeps is read off that table. Each
// plane's are the five an image of it carries, in the order its axes are reported: its own spatial
// pair, then the three both planes share. u and v play the spatial roles as l and m do (ADR 0024).
void TestTheCoordinatesAreOneTable() {
    using carta::zarr::internal::xradio::AxisCount;
    using carta::zarr::internal::xradio::CoordinateKind;
    using carta::zarr::internal::xradio::kCoordinates;
    using carta::zarr::internal::xradio::OnPlane;
    using carta::zarr::internal::xradio::Plane;
    const auto on = [](Plane plane) {
        std::vector<std::string> names;
        for (const auto& coordinate : kCoordinates) {
            if (OnPlane(coordinate, plane)) {
                names.emplace_back(coordinate.name);
            }
        }
        return names;
    };
    Require(on(Plane::sky) == std::vector<std::string>{"l", "m", "frequency", "polarization", "time"},
            "the sky plane's coordinates are not the five, in the logical order");
    Require(on(Plane::aperture) == std::vector<std::string>{"u", "v", "frequency", "polarization", "time"},
            "the aperture plane's coordinates are not u and v and the three shared, in the logical order");
    using carta::zarr::internal::xradio::PlaneCoordinate;
    Require(PlaneCoordinate(carta::zarr::AxisRole::spatial_x, Plane::aperture).name == "u" &&
                PlaneCoordinate(carta::zarr::AxisRole::spatial_y, Plane::aperture).name == "v",
            "u and v do not play the spatial roles on the aperture plane");
    Require(AxisCount(Plane::sky) == 5 && AxisCount(Plane::aperture) == 5, "a plane is not five axes");

    // What a coordinate array has to be, asked of its metadata alone. A two-dimensional one never
    // reaches the probe -- qualification closes every image whose extent cannot agree with it -- so
    // this is the only place the rule is seen to hold.
    using carta::zarr::internal::xradio::CheckCoordinate;
    using carta::zarr::internal::xradio::SkyCoordinate;
    const auto parsed = [](const std::string& document) {
        auto metadata = carta::zarr::internal::zarr::ParseArrayMetadata(nlohmann::json::parse(document), "frequency");
        Require(static_cast<bool>(metadata), "a coordinate document did not parse");
        return metadata.value();
    };
    const auto& frequency = SkyCoordinate(carta::zarr::AxisRole::spectral);
    Require(static_cast<bool>(CheckCoordinate(parsed(NumericArray("[3]", R"(["frequency"])")), frequency)),
            "a well-formed frequency coordinate was refused");
    const auto two = CheckCoordinate(parsed(NumericArray("[3,2]", R"(["frequency","x"])")), frequency);
    Require(!two && two.error().code == ErrorCode::invalid_metadata,
            "a two-dimensional frequency coordinate was accepted");
    for (const auto& coordinate : kCoordinates) {
        Require((coordinate.kind == CoordinateKind::labels) == (coordinate.name == "polarization"),
                "polarization is not the one coordinate holding labels");
    }
}

// A store of schema 0.0.2: it names the version, its base data group names SKY, SKY is typed "sky",
// and its frequency coordinate gives a frame and units. Its reference frequency carries an observer
// as well, as XRADIO writes, which is not read.
std::map<std::string, std::string> Store002() {
    auto nodes = CompleteStore();
    nodes[""] = RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"sky":"SKY"}})");
    nodes["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky"})");
    nodes["frequency"] = NumericArray("[3]", R"(["frequency"])", "float64",
                                      R"({"units":"Hz","frame":"GEO",
                                          "reference_frequency":{"data":1.4e9,"attrs":{"units":"Hz","observer":"lsrk"}}})");
    return nodes;
}

carta::zarr::internal::DescribedImage DescribeWithChannels(std::map<std::string, std::string> nodes,
                                                           const std::string& image_id = "SKY") {
    auto store = Open(std::move(nodes));
    Require(static_cast<bool>(store), "the in-memory store failed to open");
    carta::zarr::internal::xradio::CoordinateValues values;
    values.frequency = {1.4e9, 1.401e9, 1.402e9};
    auto described = carta::zarr::internal::xradio::DescribeImageFrom(store.value(), image_id, values);
    Require(static_cast<bool>(described),
            image_id + " was not described: " + (described ? "" : described.error().message));
    return std::move(described.value());
}

// A store that is what it says, 0.0.2, is read, reports its version and says nothing about it. Its
// frame is frequency's own, not the reference frequency's observer.
void TestAStoreOfSchema002IsRead() {
    const auto probe = Probe(Store002());
    Require(probe.kind == SchemaMatchKind::match, "a 0.0.2 store did not match");
    Require(probe.schema_version == "0.0.2", "a 0.0.2 store reported '" + probe.schema_version + "'");
    Require(probe.diagnostics.empty(), "a 0.0.2 store was matched with something to say");

    const auto described = DescribeWithChannels(Store002());
    Require(described.descriptor.spectral && described.descriptor.spectral->reference_frame == "GEO",
            "a 0.0.2 store's frame is frequency's own frame");
}

// What a store that is not read is refused as, by the probe, by opening the dataset and by describing
// one of its images, which must agree.
// `described` is what describing SKY on its own is refused as, when that is not what the dataset is:
// an image the store contradicts itself about is malformed metadata of the image.
void RequireRefused(const std::map<std::string, std::string>& nodes, carta::zarr::DiagnosticCode diagnostic,
                    ErrorCode error, const std::string& reason, const std::string& what,
                    std::optional<ErrorCode> described_as = std::nullopt) {
    using carta::zarr::internal::ProbeStore;
    using carta::zarr::internal::RequireOpenableDataset;

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the in-memory store failed to open");
    const auto probe = ProbeStore(store.value());
    Require(static_cast<bool>(probe) && probe.value().kind == ProbeKind::invalid_dataset, what + " was not refused");
    const auto& first = probe.value().diagnostics.front();
    Require(first.code == diagnostic, what + " was refused first for something else: " + first.message);
    Require(first.message.find(reason) != std::string::npos,
            what + " was not refused because '" + reason + "': " + first.message);
    const auto openable = RequireOpenableDataset(probe.value(), "/tmp/store");
    Require(!openable && openable.error().code == error,
            what + " was not refused as " + carta::zarr::ErrorCodeName(error));

    const auto described = carta::zarr::internal::xradio::DescribeImageFrom(
        store.value(), "SKY", carta::zarr::internal::xradio::CoordinateValues{});
    Require(!described && described.error().code == described_as.value_or(error),
            "an image of " + what + " was described, or refused as something else");
}

// A store naming no version is what XRADIO 1.2.4 and earlier wrote, and is refused for its version,
// with how to write it again.
void TestAStoreNamingNoVersionIsRefused() {
    auto nodes = CompleteStore();
    nodes[""] = RootWith("");
    RequireRefused(nodes, carta::zarr::DiagnosticCode::unsupported_schema_version,
                   ErrorCode::unsupported_schema_version, "XRADIO 1.2.5 or later", "a store naming no version");

    // Zarr of something else names none either, and is not refused at all.
    auto other = Open({{"", RootWith("")}, {"OTHER", NumericArray("[2]", R"(["x"])")}});
    Require(static_cast<bool>(other), "the in-memory store failed to open");
    const auto probe = carta::zarr::internal::ProbeStore(other.value());
    Require(static_cast<bool>(probe) && probe.value().kind == ProbeKind::zarr_without_supported_schema,
            "Zarr that is not an image and names no version was refused");
}

// A store naming 0.0.2 that is not one -- the version is a claim, which an older store upgraded on
// read can carry without what it requires -- is refused as malformed, and says which requirement it
// failed.
void TestA002StoreThatIsNotOneIsInvalid() {
    const auto invalid = [](const std::map<std::string, std::string>& nodes, const std::string& reason,
                            std::optional<ErrorCode> described_as = std::nullopt) {
        RequireRefused(nodes, carta::zarr::DiagnosticCode::invalid_metadata, ErrorCode::invalid_schema, reason,
                       "a store naming 0.0.2 that is not one", described_as);
    };

    auto no_groups = Store002();
    no_groups[""] = RootWith(R"("schema_version":"0.0.2")");
    invalid(no_groups, "it has no data_groups");

    // The sky its data group names, typed otherwise, is the store contradicting itself, and is said
    // of the image.
    auto intensity = Store002();
    intensity["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"Intensity"})");
    invalid(intensity, "Data group 'base' names it as its sky, but it is typed 'Intensity'",
            ErrorCode::invalid_metadata);

    // A flag is no image, but the group naming it as its sky says it is one: the same contradiction.
    auto flag_sky = Store002();
    flag_sky["SKY"] = SkyArray("bool", R"({"type":"flag"})");
    invalid(flag_sky, "Data group 'base' names it as its sky, but it is typed 'flag'", ErrorCode::invalid_metadata);
    // Beside an image that opens, the store matches, and the sky is listed with why it will not.
    auto beside = flag_sky;
    beside["EXTRA"] = SkyArray();
    const auto beside_probe = Probe(beside);
    Require(beside_probe.kind == SchemaMatchKind::match, "a flag-typed sky beside an openable image closed the store");
    auto beside_store = Open(beside);
    Require(static_cast<bool>(beside_store), "the in-memory store failed to open");
    const auto listed = XradioProfile().Discover(beside_store.value());
    Require(static_cast<bool>(listed) && ImageIds(listed.value().images) == std::vector<std::string>{"SKY", "EXTRA"} &&
                OpenableImageIds(listed.value().images) == std::vector<std::string>{"EXTRA"},
            "a flag-typed sky was not listed unopenable beside the image that opens");
    const auto described = XradioProfile().Describe(beside_store.value(), "SKY");
    Require(!described && described.error().code == ErrorCode::invalid_metadata,
            "describing a flag-typed sky was not refused as invalid metadata");

    auto no_frame = Store002();
    no_frame["frequency"] = FrequencyArray("[3]", R"({"units":"Hz"})");
    invalid(no_frame, "its frequency coordinate has no frame or no units");

    // Nor is a unit taken from the reference frequency in place of frequency's own.
    auto no_units = Store002();
    no_units["frequency"] =
        FrequencyArray("[3]", R"({"frame":"GEO","reference_frequency":{"data":1.4e9,"attrs":{"units":"Hz"}}})");
    invalid(no_units, "its frequency coordinate has no frame or no units");
}

// Schema 0.0.2 does not require a sky image, and XRADIO writes stores without one: an aperture plane
// alone, or a point spread function alone. Such a store is not malformed for it, and is decided by
// what it holds, as any other store is (#40). The sky's own check still applies to a store that names
// one, which TestA002StoreThatIsNotOneIsInvalid holds, and the frequency's to every store.
void TestA002StoreNamingNoSkyIsDecidedByWhatItHolds() {
    using carta::zarr::internal::ProbeStore;
    using carta::zarr::internal::RequireOpenableDataset;

    const auto probe_of = [](const std::map<std::string, std::string>& nodes) {
        auto store = Open(nodes);
        Require(static_cast<bool>(store), "the in-memory store failed to open");
        auto probe = ProbeStore(store.value());
        Require(static_cast<bool>(probe), "ProbeStore reported an error for a well-formed store");
        return std::move(probe.value());
    };
    const auto aperture_plane = R"(["time","frequency","polarization","u","v"])";

    // What XRADIO writes from a CASA image of the aperture plane: its data group names the aperture,
    // its flag is named by nothing, and it has no coordinate system, having no direction. A real one
    // is an image, and opens (ADR 0024).
    std::map<std::string, std::string> aperture{
        {"", RootWith(R"("type":"image_dataset","schema_version":"0.0.2",
                         "data_groups":{"base":{"aperture":"APERTURE"}})",
                      false)},
        {"APERTURE", SkyArray("float32", R"({"units":"Jy/beam","type":"aperture"})", aperture_plane)},
        {"FLAG_APERTURE", SkyArray("bool", R"({"type":"flag"})", aperture_plane)},
        {"time", NumericArray("[1]", R"(["time"])")},
        {"frequency", FrequencyArray()},
        {"polarization", PolarizationArray()},
        {"u", NumericArray("[4]", R"(["u"])", "float64", R"({"units":"lambda"})")},
        {"v", NumericArray("[5]", R"(["v"])", "float64", R"({"units":"lambda"})")},
    };
    const auto aperture_only = probe_of(aperture);
    Require(aperture_only.kind == ProbeKind::supported_dataset && aperture_only.default_image_id == "APERTURE",
            "a store of a real aperture-plane image alone did not open on it");

    // A complex one is not, and a store of nothing else is one this library is not for, not a
    // malformed one: it names no sky, and the schema does not require one.
    auto complex = aperture;
    complex["APERTURE"] = SkyArray("complex64", R"({"units":"Jy/beam","type":"aperture"})", aperture_plane);
    const auto complex_only = probe_of(complex);
    Require(complex_only.kind == ProbeKind::zarr_without_supported_schema,
            "a store of complex aperture-plane variables alone was not reported as one this library does not open");
    Require(!HasDiagnostic(complex_only.diagnostics, carta::zarr::DiagnosticCode::invalid_metadata),
            "a store of complex aperture-plane variables alone was called malformed");
    const auto refused = RequireOpenableDataset(complex_only, "/tmp/store");
    Require(!refused && refused.error().code == ErrorCode::unsupported_schema &&
                refused.error().message == "Complex aperture-plane variables are not openable",
            "a store of complex aperture-plane variables alone was not refused as one this library is not for");

    // Next to a sky, both open, the sky first.
    auto both = Store002();
    both[""] = RootWith(R"("schema_version":"0.0.2",
        "data_groups":{"base":{"sky":"SKY"},"aperture":{"aperture":"APERTURE"}})");
    both["APERTURE"] = SkyArray("float32", R"({"units":"Jy/beam","type":"aperture"})", aperture_plane);
    both["u"] = NumericArray("[4]", R"(["u"])", "float64", R"({"units":"lambda"})");
    both["v"] = NumericArray("[5]", R"(["v"])", "float64", R"({"units":"lambda"})");
    const auto mixed = probe_of(both);
    Require(mixed.kind == ProbeKind::supported_dataset && mixed.default_image_id == "SKY",
            "a store with a sky and an aperture did not open on the sky");
    Require(OpenableImageIds(mixed.images) == std::vector<std::string>{"SKY", "APERTURE"},
            "a store with a sky and an aperture did not offer both");

    // A point spread function alone, whose data group XRADIO's factory writes with no sky, opens.
    auto psf = Store002();
    psf.erase("SKY");
    psf[""] = RootWith(
        R"("schema_version":"0.0.2","data_groups":{"base":{"point_spread_function":"POINT_SPREAD_FUNCTION"}})");
    psf["POINT_SPREAD_FUNCTION"] = SkyArray("float32", R"({"units":"","type":"point_spread_function"})");
    const auto psf_only = probe_of(psf);
    Require(psf_only.kind == ProbeKind::supported_dataset && psf_only.default_image_id == "POINT_SPREAD_FUNCTION",
            "a store of a point spread function alone did not open on it");
    auto psf_store = Open(psf);
    Require(static_cast<bool>(psf_store), "the in-memory store failed to open");
    Require(static_cast<bool>(carta::zarr::internal::xradio::DescribeImageFrom(
                psf_store.value(), "POINT_SPREAD_FUNCTION", carta::zarr::internal::xradio::CoordinateValues{})),
            "a point spread function in a store naming no sky was not described");

    // A store with no sky is held to 0.0.2's frequency all the same: the schema requires its frame and
    // units of every store.
    auto no_frame = psf;
    no_frame["frequency"] = FrequencyArray("[3]", R"({"units":"Hz"})");
    const auto unframed = probe_of(no_frame);
    Require(unframed.kind == ProbeKind::invalid_dataset &&
                unframed.diagnostics.front().code == carta::zarr::DiagnosticCode::invalid_metadata &&
                unframed.diagnostics.front().message.find("its frequency coordinate has no frame or no units") !=
                    std::string::npos,
            "a store naming no sky whose frequency has no frame was not refused for it");
    auto unframed_store = Open(no_frame);
    Require(static_cast<bool>(unframed_store), "the in-memory store failed to open");
    const auto unframed_psf = carta::zarr::internal::xradio::DescribeImageFrom(
        unframed_store.value(), "POINT_SPREAD_FUNCTION", carta::zarr::internal::xradio::CoordinateValues{});
    Require(!unframed_psf && unframed_psf.error().code == ErrorCode::invalid_schema,
            "an image of a store naming no sky whose frequency has no frame was described");
}

// An image of the aperture plane (ADR 0024): u and v are its spatial axes, first, and it has an
// aperture coordinate measured from u's and v's crval rather than a direction. A dataset of such
// images alone needs the aperture plane's coordinates and no coordinate system; an aperture-plane
// variable written without its own coordinates does not open, and the sky images beside it do.
void TestAnApertureImageIsDescribedOnItsPlane() {
    using carta::zarr::AxisRole;
    using carta::zarr::internal::xradio::CoordinateValues;
    using carta::zarr::internal::xradio::DescribeImageFrom;
    const auto aperture_plane = R"(["time","frequency","polarization","u","v"])";
    std::map<std::string, std::string> nodes{
        {"", RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"aperture":"APERTURE"}})", false)},
        {"APERTURE", SkyArray("float32", R"({"units":"Jy","type":"aperture"})", aperture_plane)},
        {"time", NumericArray("[1]", R"(["time"])")},
        {"frequency", FrequencyArray()},
        {"polarization", PolarizationArray()},
        {"u", NumericArray("[4]", R"(["u"])", "float64", R"({"units":"lambda","crval":100.0})")},
        {"v", NumericArray("[5]", R"(["v"])", "float64", R"({"units":"lambda","crval":0.0})")},
    };
    Require(Probe(nodes).kind == SchemaMatchKind::match, "a dataset of an aperture-plane image alone did not match");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the in-memory store failed to open");
    CoordinateValues values;
    values.u = {300.0, 200.0, 100.0, 0.0};
    values.v = {-20.0, -10.0, 0.0, 10.0, 20.0};
    values.frequency = {1.4e9, 1.401e9, 1.402e9};
    const auto described = DescribeImageFrom(store.value(), "APERTURE", values);
    Require(static_cast<bool>(described),
            "an aperture-plane image was not described: " + (described ? "" : described.error().message));
    const auto& descriptor = described.value().descriptor;
    Require(descriptor.axes.at(0).name == "u" && descriptor.axes.at(0).role == AxisRole::spatial_x &&
                descriptor.axes.at(1).name == "v" && descriptor.axes.at(1).role == AxisRole::spatial_y,
            "u and v were not the first two axes, playing the spatial roles");
    Require(!descriptor.direction && descriptor.aperture, "an aperture-plane image was given a direction");
    const auto& aperture = *descriptor.aperture;
    Require(aperture.unit == "lambda" && aperture.reference_value == std::array<double, 2>{100.0, 0.0} &&
                aperture.reference_pixel == std::array<double, 2>{3.0, 3.0} &&
                aperture.increment == std::array<double, 2>{-100.0, 10.0},
            "the aperture coordinate was not measured from u's and v's crval");
    Require(descriptor.diagnostics.empty(), "an aperture coordinate with its crval was described with a diagnostic");

    // XRADIO writes crval on u and v. Without it, the reference value is taken as 0 and that is said.
    auto no_crval = nodes;
    no_crval["u"] = NumericArray("[4]", R"(["u"])", "float64", R"({"units":"lambda"})");
    auto no_crval_store = Open(no_crval);
    Require(static_cast<bool>(no_crval_store), "the in-memory store failed to open");
    const auto uncentred = DescribeImageFrom(no_crval_store.value(), "APERTURE", values);
    Require(static_cast<bool>(uncentred), "an aperture-plane image whose u has no crval was not described");
    Require(uncentred.value().descriptor.aperture->reference_value.at(0) == 0.0,
            "a u with no crval was not given 0 as its reference value");
    const auto& said = uncentred.value().descriptor.diagnostics;
    Require(std::any_of(said.begin(), said.end(),
                        [](const carta::zarr::Diagnostic& diagnostic) {
                            return diagnostic.code == carta::zarr::DiagnosticCode::invalid_metadata &&
                                   diagnostic.node_path == "u";
                        }),
            "a u with no crval was taken as 0 without saying so");

    // Its plane's coordinates are the dataset's to carry.
    auto no_frequency = nodes;
    no_frequency.erase("frequency");
    Require(Probe(no_frequency).kind == SchemaMatchKind::invalid,
            "a dataset of an aperture-plane image with no frequency coordinate was accepted");

    // Beside a sky image, an aperture-plane variable with no u or v does not open, and the sky does.
    auto beside = CompleteStore();
    beside["APERTURE"] = SkyArray("float32", R"({"type":"aperture"})", aperture_plane);
    const auto probe = Probe(beside);
    Require(probe.kind == SchemaMatchKind::match, "a sky image beside a coordinate-less aperture did not match");
    auto beside_store = Open(beside);
    Require(static_cast<bool>(beside_store), "the in-memory store failed to open");
    const auto discovery = XradioProfile().Discover(beside_store.value());
    Require(
        static_cast<bool>(discovery) && OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY"},
        "an aperture-plane variable without its coordinates was offered");
}

// A user keyword that cannot be reported -- a dictionary, as a CASA image's miscinfo can hold -- is
// left out with a diagnostic, and the image still opens with the keywords that can.
void TestAnUnreportableUserKeywordDoesNotCloseTheImage() {
    auto nodes = Store002();
    nodes["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","type":"sky",
        "user":{"instrume":"ALMA","miscinfo":{"nested":1}}})");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::match, "an image with a dictionary among its user keywords did not match");
    const auto described = DescribeWithChannels(nodes);
    const auto& descriptor = described.descriptor;
    Require(descriptor.user_keywords.size() == 1 && descriptor.user_keywords.front().name == "instrume",
            "the reportable user keyword was not the one reported");
    Require(HasDiagnostic(descriptor.diagnostics, carta::zarr::DiagnosticCode::unsupported_data_type),
            "the dictionary was left out without a diagnostic");
}

// Another version is refused as that, ahead of anything else wrong with the store, with an
// error a consumer can tell from a broken file: a newer build may open it.
void TestAnUnreadSchemaVersionIsRefused() {
    using carta::zarr::internal::ProbeStore;
    using carta::zarr::internal::RequireOpenableDataset;

    for (const std::string version : {R"("0.0.1")", R"("0.0.3")", "2"}) {
        auto nodes = Store002();
        nodes[""] = RootWith(R"("schema_version":)" + version + R"(,"data_groups":{"base":{"sky":"SKY"}})");
        auto store = Open(nodes);
        Require(static_cast<bool>(store), "the in-memory store failed to open");
        const auto probe = ProbeStore(store.value());
        Require(static_cast<bool>(probe) && probe.value().kind == ProbeKind::invalid_dataset,
                "schema version " + version + " was not refused");
        Require(probe.value().diagnostics.front().code == carta::zarr::DiagnosticCode::unsupported_schema_version,
                "the refusal of schema version " + version + " did not say it was the version");
        const auto openable = RequireOpenableDataset(probe.value(), "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::unsupported_schema_version,
                "schema version " + version + " was not refused as unsupported_schema_version");

        carta::zarr::internal::xradio::CoordinateValues values;
        const auto described = carta::zarr::internal::xradio::DescribeImageFrom(store.value(), "SKY", values);
        Require(!described && described.error().code == ErrorCode::unsupported_schema_version,
                "an image of schema version " + version + " was described");
    }
}

// XRADIO 1.2.5 and 1.3.0 keep an integer FITS image's type, BITPIX 16 as int16, and write it as it
// is, although schema 0.0.2 allows float32 and float64. What XRADIO's own conversion wrote opens.
void TestAnIntegerSkyImageXradioWroteOpens() {
    auto nodes = Store002();
    nodes["SKY"] = SkyArray("int16", R"({"units":"Jy/beam","type":"sky"})");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::match, "an int16 sky image XRADIO converted from FITS was refused");
    const auto described = DescribeWithChannels(nodes);
    Require(described.descriptor.stored_type == carta::zarr::DataType::int16, "the int16 sky image was not described");
}

// XRADIO's measurement sets name a schema_version too. Zarr of something else is not refused for it.
void TestAnotherDatasetNamingAVersionIsNotAnImage() {
    auto store = Open({{"", RootWith(R"("schema_version":"4.0.0","type":"visibility")")},
                       {"VISIBILITY", NumericArray("[2]", R"(["x"])")}});
    Require(static_cast<bool>(store), "the in-memory store failed to open");
    const auto probe = carta::zarr::internal::ProbeStore(store.value());
    Require(static_cast<bool>(probe) && probe.value().kind == ProbeKind::zarr_without_supported_schema,
            "a dataset that is not an image was refused for its schema version");
    Require(!HasDiagnostic(probe.value().diagnostics, carta::zarr::DiagnosticCode::unsupported_schema_version),
            "a dataset that is not an image was said to have an unsupported schema version");
}

// The data groups' own images are listed first: the base group's, then every other group's by group
// name. The rest follow, XRADIO's own names in its order and any other by name, and still open. A
// variable typed "model", which schema 0.0.2 has no image of, is not one.
void TestAnImagesPartAndWhereItIsListed() {
    auto nodes = Store002();
    nodes[""] = RootWith(R"("schema_version":"0.0.2",
        "data_groups":{"base":{"sky":"SKY"},"model":{"sky":"SKY_MODEL"},"dirty":{"sky":"SKY_DIRTY"}})");
    nodes["SKY_MODEL"] = SkyArray();
    nodes["SKY_DIRTY"] = SkyArray();
    nodes["MASK"] = SkyArray("float32", R"({"units":"","type":"mask"})");
    nodes["EXTRA"] = SkyArray();
    nodes["MODEL"] = SkyArray("float32", R"({"units":"Jy/beam","type":"model"})");
    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the in-memory store failed to open");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "a 0.0.2 store was not listed");
    const auto& images = discovery.value().images;
    Require(ImageIds(images) == std::vector<std::string>{"SKY", "SKY_DIRTY", "SKY_MODEL", "MASK", "EXTRA"},
            "a 0.0.2 store's images were not listed by their data groups, then in XRADIO's order, then by name");
    Require(OpenableImageIds(images) == ImageIds(images), "an image no data group names did not open");
    Require(discovery.value().default_image_id == "SKY", "the base group's sky was not the default");

    // The base group's image does not open, so the default is the next group's, by group name.
    auto base_refused = nodes;
    base_refused["SKY"] = SkyArray("complex64");
    auto refused_store = Open(base_refused);
    Require(static_cast<bool>(refused_store), "the in-memory store failed to open");
    const auto refused = XradioProfile().Discover(refused_store.value());
    Require(static_cast<bool>(refused) && refused.value().default_image_id == "SKY_DIRTY",
            "with the base group's image unopenable, the default was not the next group's");

    // Groups that share one image, and a group with none, rank no image after a group's own: ZED is
    // the z group's, so it comes before SKY_MODEL, which no group has as its own.
    auto shared = Store002();
    shared[""] = RootWith(R"("schema_version":"0.0.2",
        "data_groups":{"base":{"sky":"SKY"},"a":{"sky":"SKY"},"b":{"sky":"SKY"},"c":{"sky":"SKY"},
                       "empty":{},"z":{"sky":"ZED"}})");
    shared["SKY"] = SkyArray("complex64");
    shared["ZED"] = SkyArray();
    shared["SKY_MODEL"] = SkyArray();
    auto shared_store = Open(shared);
    Require(static_cast<bool>(shared_store), "the in-memory store failed to open");
    const auto sharing = XradioProfile().Discover(shared_store.value());
    Require(static_cast<bool>(sharing) &&
                ImageIds(sharing.value().images) == std::vector<std::string>{"SKY", "ZED", "SKY_MODEL"},
            "groups sharing an image put an image no group has as its own before a group's");
    Require(sharing.value().default_image_id == "ZED",
            "with the base group's image unopenable, the default was not the next group's own");

    // The base group's image comes first whatever its name: ZED before group a's SKY, and a base
    // aperture before another group's sky.
    auto base_first = Store002();
    base_first[""] = RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"sky":"ZED"},"a":{"sky":"SKY"}})");
    base_first["ZED"] = SkyArray();
    auto base_first_store = Open(base_first);
    Require(static_cast<bool>(base_first_store), "the in-memory store failed to open");
    const auto zed = XradioProfile().Discover(base_first_store.value());
    Require(static_cast<bool>(zed) && ImageIds(zed.value().images) == std::vector<std::string>{"ZED", "SKY"} &&
                zed.value().default_image_id == "ZED",
            "the base group's ZED was not listed and chosen before group a's SKY");
    auto base_aperture = Store002();
    base_aperture[""] =
        RootWith(R"("schema_version":"0.0.2","data_groups":{"base":{"aperture":"APERTURE"},"a":{"sky":"SKY"}})");
    base_aperture["APERTURE"] =
        SkyArray("float32", R"({"type":"aperture"})", R"(["time","frequency","polarization","u","v"])");
    base_aperture["u"] = NumericArray("[4]", R"(["u"])", "float64", R"({"units":"lambda","crval":0.0})");
    base_aperture["v"] = NumericArray("[5]", R"(["v"])", "float64", R"({"units":"lambda","crval":0.0})");
    auto base_aperture_store = Open(base_aperture);
    Require(static_cast<bool>(base_aperture_store), "the in-memory store failed to open");
    const auto aperture_first = XradioProfile().Discover(base_aperture_store.value());
    Require(static_cast<bool>(aperture_first) && aperture_first.value().default_image_id == "APERTURE",
            "the base group's aperture was not the default before another group's sky");

    // No group names an image that opens, so the default is the first of the rest.
    auto ungrouped = Store002();
    ungrouped["SKY"] = SkyArray("complex64");
    ungrouped["EXTRA"] = SkyArray();
    auto ungrouped_store = Open(ungrouped);
    Require(static_cast<bool>(ungrouped_store), "the in-memory store failed to open");
    const auto rest = XradioProfile().Discover(ungrouped_store.value());
    Require(static_cast<bool>(rest) && rest.value().default_image_id == "EXTRA",
            "an image no data group names was not the default when no group's opens");

    // Without SKY, the dataset opens on a part of the sky rather than on its point spread function.
    auto skyless = Store002();
    skyless.erase("SKY");
    skyless[""] = RootWith(R"("schema_version":"0.0.2","data_groups":{"residual":{"sky":"SKY_RESIDUAL"}})");
    skyless["SKY_RESIDUAL"] = SkyArray();
    skyless["POINT_SPREAD_FUNCTION"] = SkyArray("float32", R"({"units":"","type":"point_spread_function"})");
    auto skyless_store = Open(skyless);
    Require(static_cast<bool>(skyless_store), "the in-memory store failed to open");
    const auto skyless_discovery = XradioProfile().Discover(skyless_store.value());
    Require(static_cast<bool>(skyless_discovery) && skyless_discovery.value().default_image_id == "SKY_RESIDUAL",
            "a dataset with no SKY did not open on its residual");
}

// A continuum image has one channel, and one channel has no spacing to fit an increment from. Schema
// 0.0.2 states the channel's width, which is then the increment, measured from that channel. An axis
// of more channels is described by its channels, whatever width it states.
void TestA002ChannelWidthIsTheIncrementOfOneChannel() {
    using carta::zarr::internal::xradio::CoordinateValues;
    using carta::zarr::internal::xradio::DescribeImageFrom;

    const std::string frequency = R"({"units":"Hz","frame":"LSRK","channel_width":{"data":2.5e5,"dims":[],
        "attrs":{"units":"Hz","type":"quantity"}},"reference_frequency":{"data":1.4e9,"attrs":{"units":"Hz"}}})";
    const auto describe = [](std::map<std::string, std::string> nodes, std::vector<double> channels) {
        auto store = Open(std::move(nodes));
        Require(static_cast<bool>(store), "the in-memory store failed to open");
        CoordinateValues values;
        values.frequency = std::move(channels);
        auto described = DescribeImageFrom(store.value(), "SKY", values);
        Require(static_cast<bool>(described) && described.value().descriptor.spectral.has_value(),
                "a channel should make a spectral coordinate");
        return *described.value().descriptor.spectral;
    };

    const auto continuum_sky = [](const std::string& attributes) {
        return NumericArray("[1,1,2,4,5]", R"(["time","frequency","polarization","l","m"])", "float32", attributes);
    };
    auto continuum = Store002();
    continuum["SKY"] = continuum_sky(R"({"units":"Jy/beam","type":"sky"})");
    continuum["frequency"] = NumericArray("[1]", R"(["frequency"])", "float64", frequency);
    const auto one = describe(continuum, {1.4e9});
    Require(one.increment == 2.5e5, "a 0.0.2 channel's width was not its increment");
    Require(one.reference_pixel == 1.0 && one.reference_value == 1.4e9,
            "the increment of one channel was not measured from it");

    auto cube = Store002();
    cube["frequency"] = NumericArray("[3]", R"(["frequency"])", "float64", frequency);
    const auto three = describe(cube, {1.4e9, 1.401e9, 1.402e9});
    Require(three.increment && std::abs(*three.increment - 1.0e6) < 1.0e-3,
            "an axis of three channels was not described by its channels");
}

// What a sky image measures is read as XRADIO names it.
void TestA002SkyImageSaysWhatItMeasures() {
    auto nodes = Store002();
    nodes["SKY"] = SkyArray("float32", R"({"units":"","type":"sky","sub_type":"SpectralIndex"})");
    Require(DescribeWithChannels(nodes).descriptor.sub_type == "SpectralIndex", "a 0.0.2 sub_type was not read");
    Require(DescribeWithChannels(Store002()).descriptor.sub_type.empty(),
            "a 0.0.2 image naming no sub_type reported one");
}

}  // namespace

int main() {
    try {
        TestCompleteStoreMatches();
        TestDatasetWithoutSkyMatches();
        TestNonMatch();
        TestNothingOpenableSaysWhy();
        TestMissingCoordinateIsInvalid();
        TestCoordinateShapeMismatchIsInvalid();
        TestCoordinateDataTypeIsChecked();
        TestUndecodableLabelsAreInvalid();
        TestMalformedCoordinateSystemIsInvalid();
        TestAMissingCoordinateSystemIsInvalid();
        TestFirstFaultIsTheOnlyDiagnostic();
        TestAnImageWithoutItsPlanesAxesIsInvalid();
        TestDiscoveryClassifiesVariables();
        TestOpenableGate();
        TestAnImageDisagreeingWithACoordinateIsNotOpenable();
        TestNothingOpenableIsInvalidOnlyWhenSomethingIsMalformed();
        TestANestedGroupIsNotABrokenArray();
        TestARefusalCarriesTheVariablesOwnReason();
        TestOneRuleDecidesWhatDatasetIsOpenable();
        TestDefaultImageSkipsUnopenablePreferredImage();
        TestConsolidatedMetadataDiscovery();
        TestAConsolidatedKeyIsFiledUnderItsCanonicalName();
        TestAConsolidatedBlockThatMisnamesItsNodesIsRefused();
        TestTheInventorySaysWhatEachNodeIs();
        TestANodeThatWillNotParseIsDiagnosedNotRefused();
        TestAnArrayWithAMalformedCodecIsDiagnosedNotRefused();
        TestADeclaredFlagIsBinding();
        TestADeclaredFlagThatCannotMaskClosesTheImageInTheListing();
        TestAnImageWithAnAxisBeyondTheFiveIsNotOpenable();
        TestAnUndeclaredFlagIsTakenOnlyByItsName();
        TestADataGroupDeclaresItsSkysFlag();
        TestStoreRejections();
        TestANullConsolidatedBlockIsNoConsolidation();
        TestSizeFallsBackWhenTheStoreCannotBeMeasured();
        TestSizeRefusesAStoreWithNoArrays();
        TestADescriptionIsBuiltFromTheValuesItIsGiven();
        TestAnArrayDisagreeingWithTheRootsCopyIsRefused();
        TestAnArrayIsNotHeldToHowItsChunksAreEncoded();
        TestCheckingAnArrayReadsItsOwnDocumentOnce();
        TestDescribingAnImageHoldsItsArraysToTheirOwnDocuments();
        TestAPixelReadHoldsItsArrayToItsOwnDocument();
        TestTheCoordinatesAreOneTable();
        TestAStoreOfSchema002IsRead();
        TestAStoreNamingNoVersionIsRefused();
        TestA002StoreThatIsNotOneIsInvalid();
        TestA002StoreNamingNoSkyIsDecidedByWhatItHolds();
        TestAnUnreportableUserKeywordDoesNotCloseTheImage();
        TestAnApertureImageIsDescribedOnItsPlane();
        TestAnUnreadSchemaVersionIsRefused();
        TestAnIntegerSkyImageXradioWroteOpens();
        TestAnotherDatasetNamingAVersionIsNotAnImage();
        TestAnImagesPartAndWhereItIsListed();
        TestA002ChannelWidthIsTheIncrementOfOneChannel();
        TestA002SkyImageSaysWhatItMeasures();
        std::cout << "carta-zarr schema profile tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr schema profile tests failed: " << error.what() << '\n';
        return 1;
    }
}
