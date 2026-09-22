/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
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
#include "schema/xradio/flag.h"
#include "store.h"

#include "support/check.h"

#include "support/in_memory_transport.h"

#include <algorithm>
#include <chrono>
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::ProbeKind;
using carta::zarr::SchemaMatchKind;
using carta::zarr::internal::xradio::DetermineFlag;
using carta::zarr::testing::MakeInMemoryTransport;

using carta::zarr::testing::Require;

bool HasDiagnostic(const carta::zarr::SchemaProbeResult& probe, const std::string& code);

bool HasDiagnostic(const std::vector<carta::zarr::Diagnostic>& diagnostics, const std::string& code) {
    for (const auto& diagnostic : diagnostics) {
        if (diagnostic.code == code) {
            return true;
        }
    }
    return false;
}

bool HasDiagnostic(const carta::zarr::SchemaProbeResult& probe, const std::string& code) {
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

std::string RootGroup(bool coordinate_system = true) {
    std::string attributes;
    if (coordinate_system) {
        attributes += R"("coordinate_system_info": {
      "projection": "SIN",
      "reference_direction": {"data": [1.0, 0.5]},
      "native_pole_direction": {"data": [0.0, 1.5707963267948966]},
      "pixel_coordinate_transformation_matrix": [[1.0, 0.0], [0.0, 1.0]]
    })";
    }
    return "{\"attributes\":{" + attributes + "},\"zarr_format\":3,\"node_type\":\"group\"}";
}

std::string NumericArray(const std::string& shape, const std::string& dimensions,
                         const std::string& data_type = "float64", const std::string& attributes = "{}") {
    return "{\"shape\":" + shape + ",\"data_type\":\"" + data_type +
           "\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":" + shape +
           "}},\"attributes\":" + attributes + ",\"dimension_names\":" + dimensions +
           ",\"zarr_format\":3,\"node_type\":\"array\"}";
}

std::string SkyArray(const std::string& data_type = "float32", const std::string& attributes = R"({"units":"Jy/beam"})",
                     const std::string& dimensions = R"(["time","frequency","polarization","l","m"])") {
    return NumericArray("[1,3,2,4,5]", dimensions, data_type, attributes);
}

std::string PolarizationArray(const std::string& data_type = R"({"name": "fixed_length_utf32",
                                    "configuration": {"length_bytes": 4}})") {
    return "{\"shape\":[2],\"data_type\":" + data_type +
           ",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":[2]}},"
           "\"attributes\":{},\"dimension_names\":[\"polarization\"],\"zarr_format\":3,\"node_type\":\"array\"}";
}

// A complete image dataset carries a coordinate array for every axis its image uses, time included.
std::map<std::string, std::string> CompleteStore() {
    return {
        {"", RootGroup()},
        {"SKY", SkyArray()},
        {"time", NumericArray("[1]", R"(["time"])")},
        {"frequency", NumericArray("[3]", R"(["frequency"])")},
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
std::map<std::string, std::string> ConsolidatedStore() {
    auto nodes = CompleteStore();

    std::string metadata;
    for (const auto& child : nodes) {
        if (child.first.empty()) {
            continue;
        }
        if (!metadata.empty()) {
            metadata += ",";
        }
        metadata += "\"" + child.first + "\":" + child.second;
    }
    nodes[""] = "{\"attributes\":{\"coordinate_system_info\":{"
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
    Require(probe.schema_version == "1.2", "unexpected schema version");
}

void TestDatasetWithoutSkyMatches() {
    auto nodes = CompleteStore();
    nodes.erase("SKY");
    nodes["RESIDUAL"] = SkyArray();
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::match,
            "an image dataset with only RESIDUAL was incorrectly rejected");
}

void TestNonMatch() {
    auto store = Open({{"", RootGroup()}, {"OTHER", NumericArray("[2]", R"(["x"])")}});
    Require(static_cast<bool>(store), "a valid non-XRADIO group failed to open");
    const auto probe = carta::zarr::internal::ProbeStore(store.value());
    Require(static_cast<bool>(probe), "ProbeStore reported an error for a valid non-XRADIO group");
    Require(probe.value().kind == ProbeKind::zarr_without_supported_schema,
            "a valid non-XRADIO group was not reported as an unsupported schema");
}

void TestMissingCoordinateIsInvalid() {
    auto nodes = CompleteStore();
    nodes.erase("time");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid,
            "an image dataset missing its time coordinate was not reported as invalid");
    Require(HasDiagnostic(probe.diagnostics, "invalid_metadata"), "the missing coordinate produced no diagnostic");
}

void TestCoordinateShapeMismatchIsInvalid() {
    auto nodes = CompleteStore();
    nodes["frequency"] = NumericArray("[7]", R"(["frequency"])");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a coordinate whose length disagrees with the image was accepted");
}

void TestCoordinateDataTypeIsChecked() {
    auto nodes = CompleteStore();
    // Polarization labels are strings; a numeric polarization coordinate is malformed.
    nodes["polarization"] = PolarizationArray("\"float64\"");
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a numeric polarization coordinate was accepted");
    Require(HasDiagnostic(probe.diagnostics, "unsupported_data_type"),
            "the polarization data type produced no diagnostic");
}

void TestMalformedCoordinateSystemIsInvalid() {
    auto nodes = CompleteStore();
    nodes[""] = R"({"attributes":{"coordinate_system_info":{}},"zarr_format":3,"node_type":"group"})";
    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "malformed coordinate_system_info was accepted");
}

// Discovery only recognizes complete image planes. An incomplete SKY-only store is simply a
// non-match because it offers no complete image plane.
void TestIncompleteImageIsNotMatch() {
    auto missing_axis = CompleteStore();
    missing_axis["SKY"] =
        SkyArray("float32", R"({"units":"Jy/beam"})", R"(["time","frequency","polarization","l","x"])");
    const auto renamed = Probe(missing_axis);
    Require(renamed.kind == SchemaMatchKind::no_match, "an incomplete SKY was treated as an image dataset");

    auto wrong_rank = CompleteStore();
    wrong_rank["SKY"] =
        NumericArray("[1,3,2,4]", R"(["time","frequency","polarization","l"])", "float32", R"({"units":"Jy/beam"})");
    const auto four = Probe(wrong_rank);
    Require(four.kind == SchemaMatchKind::no_match, "a four-dimensional SKY was treated as an image dataset");
}

// Discovery reports every variable carrying a whole plane's axes, and says why the ones it cannot
// open are closed. Flags and the optional (l, m) coordinates are never images at all.
void TestDiscoveryClassifiesVariables() {
    auto nodes = CompleteStore();
    nodes["MODEL"] = SkyArray();
    nodes["COMPLEX"] = SkyArray("complex64");
    nodes["APERTURE"] = SkyArray("float32", "{}", R"(["time","frequency","polarization","u","v"])");
    nodes["MASK_0"] = SkyArray("bool", R"({"type":"flag"})");
    nodes["right_ascension"] = NumericArray("[4,5]", R"(["l","m"])");
    nodes["u"] = NumericArray("[4]", R"(["u"])");
    nodes["v"] = NumericArray("[5]", R"(["v"])");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the discovery store failed to open");
    auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "profile discovery reported an error");

    Require(ImageIds(discovery.value().images) == std::vector<std::string>{"SKY", "MODEL", "APERTURE", "COMPLEX"},
            "discovery did not enumerate the images in display order");
    Require(OpenableImageIds(discovery.value().images) == std::vector<std::string>{"SKY", "MODEL"},
            "discovery did not restrict the openable images to real sky-plane variables");
    Require(discovery.value().default_image_id == "SKY", "discovery did not select the default openable image");
    Require(HasDiagnostic(discovery.value().diagnostics, "unsupported_coordinate_plane"),
            "the aperture-plane variable produced no diagnostic");
    Require(HasDiagnostic(discovery.value().diagnostics, "unsupported_data_type"),
            "the complex variable produced no diagnostic");
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

// One rule decides what is openable, and both sides of the library ask it. Dataset::OpenImage asks
// it of the listing the dataset kept from opening; SchemaProfile::Describe asks it of a store it
// has just inspected. They used to be two copies, and they disagreed about the message -- the
// facade passed the variable's own diagnostic through and the profile answered generically -- so a
// consumer got a worse explanation depending on which way it had arrived.
void TestOneRuleDecidesWhatIsOpenable() {
    auto nodes = CompleteStore();
    nodes["COMPLEX"] = SkyArray("complex64");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the openable-rule store failed to open");
    const auto discovery = XradioProfile().Discover(store.value());
    Require(static_cast<bool>(discovery), "discovery failed on the openable-rule store");
    const auto& images = discovery.value().images;

    Require(static_cast<bool>(carta::zarr::internal::RequireOpenable(images, "SKY")),
            "an openable image was refused");

    const auto complex = carta::zarr::internal::RequireOpenable(images, "COMPLEX");
    Require(!complex && complex.error().code == ErrorCode::unsupported_data_type,
            "a listed but unopenable variable should be refused as an unsupported data type");

    // What the profile said about this variable, rather than a generic refusal. This is the half
    // the two copies used to differ on.
    const auto listed = std::find_if(images.begin(), images.end(),
                                     [](const auto& image) { return image.id == "COMPLEX"; });
    Require(listed != images.end() && !listed->diagnostics.empty(),
            "the complex variable was not listed with a diagnostic to pass on");
    Require(complex.error().message == listed->diagnostics.front().message,
            "the refusal did not carry the variable's own diagnostic");

    const auto absent = carta::zarr::internal::RequireOpenable(images, "NOPE");
    Require(!absent && absent.error().code == ErrorCode::not_found,
            "a variable that was never listed should be reported as missing, not as unopenable");
}

// The dataset-level counterpart, and the three answers it has to keep apart. This used to be
// written out in Dataset::Open, so a probe's refusal meant whatever the facade decided it meant --
// and none of it was reachable without a store on disk shaped to produce each kind.
void TestOneRuleDecidesWhatDatasetIsOpenable() {
    using carta::zarr::internal::RequireOpenableDataset;

    // Nothing matched. Not an error about the file's contents: this library is not for it.
    {
        carta::zarr::ProbeResult probe;
        probe.kind = carta::zarr::ProbeKind::zarr_without_supported_schema;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::unsupported_schema,
                "an unmatched store should be refused as an unsupported schema");
    }

    // Something matched and was malformed. A different answer, and a consumer acts on it.
    {
        carta::zarr::ProbeResult probe;
        probe.kind = carta::zarr::ProbeKind::invalid_dataset;
        probe.diagnostics.push_back(carta::zarr::Diagnostic{"missing_coordinate", "frequency is missing", "SKY"});
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::invalid_metadata,
                "a malformed store should be refused as invalid metadata");
        Require(openable.error().message == "frequency is missing",
                "and it should carry what the probe worked out, not a generic refusal");
    }

    // A probe that had nothing to say falls back rather than reporting an empty message.
    {
        carta::zarr::ProbeResult probe;
        probe.kind = carta::zarr::ProbeKind::not_zarr;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && !openable.error().message.empty(), "a silent probe still needs a message");
    }

    // A profile matched a store with nothing openable in it. The profile is not refusing, so there
    // are no diagnostics; opening it would hand back a dataset a consumer can do nothing with.
    {
        carta::zarr::ProbeResult probe;
        probe.kind = carta::zarr::ProbeKind::supported_dataset;
        const auto openable = RequireOpenableDataset(probe, "/tmp/store");
        Require(!openable && openable.error().code == ErrorCode::invalid_metadata,
                "a matched store with no images should be refused");
    }

    // And the one case that opens.
    {
        carta::zarr::ProbeResult probe;
        probe.kind = carta::zarr::ProbeKind::supported_dataset;
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
                HasDiagnostic(sky->diagnostics, "unsupported_data_type"),
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
    Require(plain->nodes_read().size() > 1,
            "a store without consolidated metadata has to read its children");
}

// The report latches: once a requirement is unmet, later ones are no-ops. A store with two faults
// is therefore diagnosed once, by the first fault reached -- the behaviour a probe had when every
// check returned early, now stated somewhere rather than emerging from the control flow.
void TestFirstFaultIsTheOnlyDiagnostic() {
    auto nodes = CompleteStore();
    nodes["frequency"] = NumericArray("[7]", R"(["frequency"])");  // wrong length
    nodes["polarization"] = PolarizationArray("\"float64\"");      // wrong data type

    const auto probe = Probe(nodes);
    Require(probe.kind == SchemaMatchKind::invalid, "a store with two faults was not reported as invalid");
    Require(probe.diagnostics.size() == 1, "the report diagnosed more than the first unmet requirement");
    Require(probe.diagnostics.front().node_path == "frequency",
            "the diagnostic did not come from the first requirement checked");
    Require(!HasDiagnostic(probe, "unsupported_data_type"), "a requirement after the first failure was still checked");
}

// Store-level rejections, reported before any profile is consulted.
// Choosing a flag is the half of it that consults the store: what a flag has to be once it is found
// takes no store at all and is checked in tests/flag_test.cc.
//
// A declared flag is the image's own statement that its pixels need a mask, so an image naming one
// that cannot be read is closed rather than opened unmasked. Reads apply the mask by default, and an
// unusable mask reported as no mask would show flagged pixels as valid -- the one failure a consumer
// has no way to notice.
void TestADeclaredFlagIsBinding() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    const auto with_declared_flag = [&](bool write_the_flag) {
        auto nodes = CompleteStore();
        nodes["SKY"] = SkyArray("float32", R"({"units":"Jy/beam","flag":"MASK_0"})");
        if (write_the_flag) {
            nodes["MASK_0"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
        }
        return nodes;
    };

    auto present = Open(with_declared_flag(true));
    Require(static_cast<bool>(present), "the declared-flag store did not open");
    const auto& present_image = present.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(present_image), "SKY was not readable in the declared-flag store");
    std::vector<carta::zarr::Diagnostic> diagnostics;
    auto declared = DetermineFlag(present.value(), present_image.value(), "SKY", diagnostics);
    Require(static_cast<bool>(declared), "a well-formed declared flag was refused");
    Require(declared.value() == "MASK_0", "the declared flag was not the one selected");
    Require(diagnostics.empty(), "selecting a declared flag produced a diagnostic");

    auto absent = Open(with_declared_flag(false));
    Require(static_cast<bool>(absent), "the missing-flag store did not open");
    const auto& absent_image = absent.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(absent_image), "SKY was not readable in the missing-flag store");
    auto missing = DetermineFlag(absent.value(), absent_image.value(), "SKY", diagnostics);
    Require(!missing, "an image declaring a flag variable that does not exist was opened");
}

// With nothing declared the store is inspected instead, and a store offering two equally good
// candidates is refused rather than guessed at. The refusal is a diagnostic on the image: the image
// is still readable, just unmasked.
void TestAmbiguousFlagsSelectNone() {
    const std::string sky_dimensions = R"(["time","frequency","polarization","l","m"])";
    auto nodes = CompleteStore();
    nodes["FLAG_1"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");
    nodes["FLAG_2"] = NumericArray("[1,3,2,4,5]", sky_dimensions, "bool", R"({"type":"flag"})");

    auto store = Open(nodes);
    Require(static_cast<bool>(store), "the ambiguous-flag store did not open");
    const auto& image = store.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(image), "SKY was not readable in the ambiguous-flag store");

    std::vector<carta::zarr::Diagnostic> diagnostics;
    auto chosen = DetermineFlag(store.value(), image.value(), "SKY", diagnostics);
    Require(static_cast<bool>(chosen), "ambiguous flags reported an error rather than no mask");
    Require(chosen.value().empty(), "ambiguous flags selected a pixel mask");
    Require(HasDiagnostic(diagnostics, "ambiguous_pixel_mask"), "ambiguous flags did not produce a diagnostic");

    // One candidate is not ambiguous: the same store with FLAG_2 removed selects FLAG_1.
    nodes.erase("FLAG_2");
    auto single = Open(nodes);
    Require(static_cast<bool>(single), "the single-flag store did not open");
    const auto& single_image = single.value().ReadArrayMetadata("SKY");
    Require(static_cast<bool>(single_image), "SKY was not readable in the single-flag store");
    std::vector<carta::zarr::Diagnostic> single_diagnostics;
    auto only = DetermineFlag(single.value(), single_image.value(), "SKY", single_diagnostics);
    Require(static_cast<bool>(only) && only.value() == "FLAG_1", "one matching flag was not selected");
    Require(single_diagnostics.empty(), "one matching flag produced an ambiguity diagnostic");
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
// The measured half of this answer used to be the facade's: it re-parsed the location string and
// walked the directory itself, so the fallback could only be reached from here by pointing at a real
// store and setting the timeout to zero to make the walk give up. An in-memory transport gives up
// honestly instead -- it holds no bytes -- so the branch is reachable with a map.
void TestSizeFallsBackWhenTheStoreCannotBeMeasured() {
    auto store = Open(CompleteStore());
    Require(static_cast<bool>(store), "OpenStore rejected the in-memory store");

    const auto size = carta::zarr::internal::DatasetSizeBytes(store.value(),
                                                              std::chrono::steady_clock::now() +
                                                                  std::chrono::seconds(5));
    Require(static_cast<bool>(size), "sizing an in-memory store failed");
    Require(size.value().basis == carta::zarr::SizeBasis::declared,
            "a size the transport could not measure must be reported as the declared one");
    // SKY is 120 float32 at 480 bytes; time, frequency, l and m are 1, 3, 4 and 5 float64 at 8, 24,
    // 32 and 40; polarization is two four-byte labels at 8. The arrays, not the store.
    Require(size.value().bytes == 592,
            "the logical total was " + std::to_string(size.value().bytes) + ", not 592");

    // A deadline that has already passed reaches the same answer by the same route: this transport
    // refuses whatever the clock says, and every way of failing to measure means the declared size.
    const auto expired = carta::zarr::internal::DatasetSizeBytes(
        store.value(), std::chrono::steady_clock::now() - std::chrono::seconds(1));
    Require(expired && expired.value().basis == carta::zarr::SizeBasis::declared &&
                expired.value().bytes == 592,
            "an expired deadline did not reach the same declared size");
}

// And the fallback does not turn every failure into a number: a store with nothing to add up is
// still an error, not a zero.
void TestSizeRefusesAStoreWithNoArrays() {
    auto store = Open({{"", RootGroup()}});
    Require(static_cast<bool>(store), "OpenStore rejected a bare root group");

    const auto size = carta::zarr::internal::DatasetSizeBytes(store.value(),
                                                              std::chrono::steady_clock::now());
    Require(!size && size.error().code == ErrorCode::invalid_metadata,
            "a store holding no arrays was given a size");
}

}  // namespace

int main() {
    try {
        TestCompleteStoreMatches();
        TestDatasetWithoutSkyMatches();
        TestNonMatch();
        TestMissingCoordinateIsInvalid();
        TestCoordinateShapeMismatchIsInvalid();
        TestCoordinateDataTypeIsChecked();
        TestMalformedCoordinateSystemIsInvalid();
        TestFirstFaultIsTheOnlyDiagnostic();
        TestIncompleteImageIsNotMatch();
        TestDiscoveryClassifiesVariables();
        TestOpenableGate();
        TestOneRuleDecidesWhatIsOpenable();
        TestOneRuleDecidesWhatDatasetIsOpenable();
        TestDefaultImageSkipsUnopenablePreferredImage();
        TestConsolidatedMetadataDiscovery();
        TestADeclaredFlagIsBinding();
        TestAmbiguousFlagsSelectNone();
        TestStoreRejections();
        TestSizeFallsBackWhenTheStoreCannotBeMeasured();
        TestSizeRefusesAStoreWithNoArrays();
        std::cout << "carta-zarr schema profile tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr schema profile tests failed: " << error.what() << '\n';
        return 1;
    }
}
