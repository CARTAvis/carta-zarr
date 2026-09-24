/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <carta-zarr/carta_zarr.h>

#include "support/check.h"

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::SchemaMatchKind;

using carta::zarr::testing::Require;

void Write(const std::filesystem::path& path, const std::string& contents) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    Require(output.is_open(), "Unable to write " + path.string());
    output << contents;
    Require(static_cast<bool>(output), "Unable to finish writing " + path.string());
}

void WriteDoubles(const std::filesystem::path& path, const std::vector<double>& values) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    Require(output.is_open(), "Unable to write " + path.string());
    output.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(double)));
    Require(static_cast<bool>(output), "Unable to finish writing " + path.string());
}

void WriteFloats(const std::filesystem::path& path, const std::vector<float>& values) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    Require(output.is_open(), "Unable to write " + path.string());
    output.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(values.size() * sizeof(float)));
    Require(static_cast<bool>(output), "Unable to finish writing " + path.string());
}

void WriteUtf32(const std::filesystem::path& path, const std::vector<std::string>& values, std::size_t code_points) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary);
    Require(output.is_open(), "Unable to write " + path.string());
    for (const auto& value : values) {
        for (std::size_t index = 0; index < code_points; ++index) {
            const std::uint32_t code_point = index < value.size() ? static_cast<std::uint32_t>(value.at(index)) : 0U;
            output.write(reinterpret_cast<const char*>(&code_point), sizeof(code_point));
        }
    }
    Require(static_cast<bool>(output), "Unable to finish writing " + path.string());
}

bool HasDiagnostic(const std::vector<carta::zarr::Diagnostic>& diagnostics, const std::string& code) {
    return std::any_of(diagnostics.begin(), diagnostics.end(),
                       [&](const auto& diagnostic) { return diagnostic.code == code; });
}

std::vector<std::string> ImageIds(const std::vector<carta::zarr::ImageEntry>& images) {
    std::vector<std::string> ids;
    ids.reserve(images.size());
    for (const auto& image : images) {
        ids.push_back(image.id);
    }
    return ids;
}

std::string RootMetadata(bool coordinate_system = true) {
    return coordinate_system ? R"({
  "attributes": {
    "coordinate_system_info": {
      "projection": "SIN",
      "reference_direction": {"data": [1.0, 0.5]},
      "native_pole_direction": {"data": [0.0, 1.5707963267948966]},
      "pixel_coordinate_transformation_matrix": [[1.0, 0.0], [0.0, 1.0]]
    }
  },
  "zarr_format": 3,
  "node_type": "group"
})"
                             : R"({"zarr_format": 3, "node_type": "group"})";
}

std::string NumericArray(const std::string& shape, const std::string& dimensions,
                         const std::string& data_type = "float64", const std::string& attributes = "{}") {
    return "{\"shape\":" + shape + ",\"data_type\":\"" + data_type +
           "\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":" + shape +
           "}},\"chunk_key_encoding\":{\"name\":\"default\",\"configuration\":{\"separator\":\"/\"}},"
           "\"fill_value\":0,\"codecs\":[{\"name\":\"bytes\",\"configuration\":{\"endian\":\"little\"}}],"
           "\"attributes\":" +
           attributes + ",\"dimension_names\":" + dimensions + ",\"zarr_format\":3,\"node_type\":\"array\"}";
}

std::string SkyArray(std::uint64_t time_length = 1, bool arbitrary_storage_order = false) {
    if (arbitrary_storage_order) {
        return "{\"shape\":[" + std::to_string(time_length) +
               ",5,3,4,2],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_"
               "shape\":[1,5,1,2,1]}},\"attributes\":{\"units\":\"Jy/"
               "beam\"},\"dimension_names\":[\"time\",\"m\",\"frequency\",\"l\",\"polarization\"],\"zarr_format\":3,"
               "\"node_type\":\"array\"}";
    }
    return "{\"shape\":[" + std::to_string(time_length) +
           ",3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_"
           "shape\":[1,1,1,2,5]}},\"attributes\":{\"units\":\"Jy/"
           "beam\"},\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],\"zarr_format\":3,\"node_"
           "type\":\"array\"}";
}

std::string PolarizationArray() {
    return R"({
  "shape": [2],
  "data_type": {"name": "fixed_length_utf32", "configuration": {"length_bytes": 4}},
  "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [2]}},
  "attributes": {"dimension_names": ["polarization"]},
  "dimension_names": ["polarization"],
  "zarr_format": 3,
  "node_type": "array"
})";
}

void CreateValidStore(const std::filesystem::path& root, std::uint64_t time_length = 1, bool coordinate_system = true,
                      bool arbitrary_storage_order = false) {
    Write(root / "zarr.json", RootMetadata(coordinate_system));
    Write(root / "SKY" / "zarr.json", SkyArray(time_length, arbitrary_storage_order));
    Write(root / "time" / "zarr.json", NumericArray("[" + std::to_string(time_length) + "]", "[\"time\"]"));
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", "[\"frequency\"]"));
    Write(root / "l" / "zarr.json", NumericArray("[4]", "[\"l\"]"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", "[\"m\"]"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
}

void TestValidAndTimeAxis(const std::filesystem::path& root) {
    CreateValidStore(root, 1, true, true);
    const auto result = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(result && result.value().kind == SchemaMatchKind::match, "valid XRADIO store did not match");
    Require(result.value().schema_id == carta::zarr::kXradioImageSchema, "unexpected schema id");
    Require(result.value().schema_version == "1.2", "unexpected schema version");

    const auto is_xradio = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(is_xradio && is_xradio.value().kind == SchemaMatchKind::match,
            "a valid store was not matched by the XRADIO profile");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed");
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"SKY"}, "unexpected image ids");
    Require(dataset.value().descriptor().default_image_id == "SKY", "unexpected default image id");

    const auto declared_size = dataset.value().Size(std::chrono::milliseconds(0));
    Require(declared_size && declared_size.value().bytes == 592 &&
                declared_size.value().basis == carta::zarr::SizeBasis::declared,
            "the declared Zarr size calculation was incorrect");

    const auto measured_size = dataset.value().Size(std::chrono::milliseconds(5000));
    Require(measured_size && measured_size.value().bytes > 0 &&
                measured_size.value().basis == carta::zarr::SizeBasis::measured,
            "the measured Zarr size calculation was not used");

    // The two are different questions, and the declared one does not bound the other: this store's
    // arrays declare 592 bytes while the store holds several times that in zarr.json documents
    // alone. Asserted on the store the rest of this case already built, because the point is that
    // nothing about it is unusual -- it is what the relationship between the two is worth.
    Require(declared_size.value().bytes < measured_size.value().bytes,
            "the declared size was not below the measured one, so this store no longer shows that "
            "the declared size is not an upper bound");

    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image),
            "Dataset::OpenImage failed" + (image ? std::string{} : ": " + image.error().message));
    Require(image.value().descriptor().axes.size() == 5, "time axis was not preserved");
    Require(image.value().descriptor().axes.back().length == 1, "unexpected singleton time axis");
    Require(image.value().descriptor().axes.at(0).storage_index == 3 &&
                image.value().descriptor().axes.at(1).storage_index == 1 &&
                image.value().descriptor().axes.at(2).storage_index == 2 &&
                image.value().descriptor().axes.at(3).storage_index == 4 &&
                image.value().descriptor().axes.at(4).storage_index == 0,
            "logical axes were not mapped to the arbitrary storage order");

    // Verify Direction & Coordinates
    const auto& desc = image.value().descriptor();
    Require(desc.direction.has_value(), "DirectionCoordinate missing");
    Require(desc.direction->projection == "SIN", "unexpected projection");
    Require(desc.spectral.has_value(), "SpectralCoordinate missing");
    Require(desc.observation.has_value(), "ObservationInfo missing");
}

void TestTimeGreaterThanOne(const std::filesystem::path& root) {
    CreateValidStore(root, 2);
    const auto result = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(result && result.value().kind == SchemaMatchKind::match,
            "time > 1 must remain valid at library level");
    const auto context = carta::zarr::Context::Create();
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open rejected time > 1");
    const auto image = dataset.value().OpenImage("SKY");
    Require(image && image.value().descriptor().axes.back().length == 2, "time axis length was not preserved");
}

void TestNonMatchAndInvalid(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "OTHER" / "zarr.json", NumericArray("[2]", "[\"x\"]"));
    const auto non_match = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(non_match && non_match.value().kind == SchemaMatchKind::no_match,
            "valid non-XRADIO Zarr was not a non-match");

    CreateValidStore(root);
    std::filesystem::remove(root / "time" / "zarr.json");
    const auto invalid = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(invalid && invalid.value().kind == SchemaMatchKind::invalid,
            "metadata-incomplete XRADIO-like store was not reported as an invalid match");
}

// A store Dataset::Open refuses must say why it refused. The probe has already worked the reason
// out -- often down to the attribute -- and the consumer puts this message in front of whoever
// picked the file, so answering "not a supported dataset" spends a probe and reports nothing.
void TestOpenSaysWhyItRefused(const std::filesystem::path& root) {
    const auto context = carta::zarr::Context::Create();

    // Zarr that no built-in profile claims: nothing is wrong with it, so there is no diagnostic to
    // pass on and the message has to stand on its own.
    const auto unclaimed = root / "unclaimed";
    Write(unclaimed / "zarr.json", RootMetadata());
    Write(unclaimed / "OTHER" / "zarr.json", NumericArray("[2]", "[\"x\"]"));
    const auto not_ours = carta::zarr::Dataset::Open(context.value(), unclaimed.string());
    Require(!not_ours && not_ours.error().code == ErrorCode::unsupported_schema,
            "an unclaimed Zarr store was not reported as an unsupported schema");
    Require(!not_ours.error().message.empty(), "an unclaimed Zarr store was refused without a message");
    const auto unclaimed_probe = carta::zarr::ProbeSchema(unclaimed.string(), carta::zarr::kXradioImageSchema);
    Require(unclaimed_probe && unclaimed_probe.value().diagnostics.empty(),
            "an unclaimed Zarr store was reported with a diagnostic nothing produced");

    // Zarr that no profile claims either, but whose only image the profile recognised and will not
    // open. That is a reason, and the caller is told it rather than that nothing matched.
    const auto complex = root / "complex";
    CreateValidStore(complex);
    auto complex_sky = SkyArray();
    complex_sky.replace(complex_sky.find("float32"), std::string("float32").size(), "complex64");
    Write(complex / "SKY" / "zarr.json", complex_sky);
    const auto refused = carta::zarr::Dataset::Open(context.value(), complex.string());
    Require(!refused && refused.error().code == ErrorCode::unsupported_schema,
            "a store of complex images was not reported as an unsupported schema");
    Require(refused.error().message == "Complex sky-plane variables are not openable",
            "Dataset::Open did not say why a store of complex images was refused: " + refused.error().message);

    // A store this profile claims and then finds malformed. Here the probe does have something to
    // say, and it is what the caller must be told.
    const auto malformed = root / "malformed";
    CreateValidStore(malformed);
    std::filesystem::remove(malformed / "time" / "zarr.json");
    const auto probe = carta::zarr::ProbeSchema(malformed.string(), carta::zarr::kXradioImageSchema);
    Require(probe && probe.value().kind == SchemaMatchKind::invalid && !probe.value().diagnostics.empty(),
            "the malformed store did not probe as invalid with a diagnostic");
    const auto opened = carta::zarr::Dataset::Open(context.value(), malformed.string());
    Require(!opened && opened.error().code == ErrorCode::invalid_metadata,
            "a malformed store was not reported as invalid metadata");
    Require(opened.error().message == probe.value().diagnostics.front().message,
            "Dataset::Open replaced the probe's diagnostic with a message of its own: " +
                opened.error().message);
}

void TestMissingAndUnsupported(const std::filesystem::path& root) {
    const auto missing = carta::zarr::ProbeSchema((root / "missing").string(), carta::zarr::kXradioImageSchema);
    Require(!missing && missing.error().code == ErrorCode::not_found, "missing store error category changed");

    Write(root / "not-zarr" / "zarr.json", "{\"zarr_format\": 2, \"node_type\": \"group\"}");
    const auto unsupported = carta::zarr::ProbeSchema((root / "not-zarr").string(), carta::zarr::kXradioImageSchema);
    Require(!unsupported && unsupported.error().code == ErrorCode::unsupported_zarr_version,
            "unsupported Zarr version error category changed");

    const auto unknown_schema = carta::zarr::ProbeSchema(root.string(), "future.schema");
    Require(!unknown_schema && unknown_schema.error().code == ErrorCode::unsupported_schema,
            "unknown schema error category changed");
}

// What the transport does with a document it can open but whose bytes say nothing. An empty
// zarr.json is not an I/O failure -- the read worked, the file is simply not metadata -- and the
// answer has to come from the parser above the seam rather than from the read.
void TestEmptyAndOversizedMetadata(const std::filesystem::path& root) {
    CreateValidStore(root / "empty-root");
    Write(root / "empty-root" / "zarr.json", "");
    const auto empty_root = carta::zarr::ProbeSchema((root / "empty-root").string(), carta::zarr::kXradioImageSchema);
    Require(!empty_root && empty_root.error().code == ErrorCode::invalid_metadata,
            "an empty root zarr.json should be invalid metadata, not an I/O failure");

    // The same one node down, where it is reached through the node metadata table instead.
    CreateValidStore(root / "empty-child");
    Write(root / "empty-child" / "frequency" / "zarr.json", "");
    const auto empty_child = carta::zarr::ProbeSchema((root / "empty-child").string(), carta::zarr::kXradioImageSchema);
    Require(!empty_child && empty_child.error().code == ErrorCode::invalid_metadata,
            "an empty child zarr.json should be invalid metadata, not an I/O failure");

    // A document larger than any buffer the reader might have sized for one: it is read whole, and
    // the store it describes opens.
    CreateValidStore(root / "padded");
    std::string padded = RootMetadata();
    padded.insert(padded.size() - 1, ",\n  \"note\": \"" + std::string(400000, 'x') + "\"");
    Write(root / "padded" / "zarr.json", padded);
    const auto probe = carta::zarr::ProbeSchema((root / "padded").string(), carta::zarr::kXradioImageSchema);
    Require(probe && probe.value().kind == SchemaMatchKind::match,
            "a root document of several hundred kilobytes was not read whole");
}

void TestReferenceFixture() {
    const std::filesystem::path fixture(CARTA_ZARR_REFERENCE_FIXTURE);
    Require(std::filesystem::exists(fixture), "the XRADIO reference fixture is missing from tests/data");

    const auto result = carta::zarr::ProbeSchema(fixture.string(), carta::zarr::kXradioImageSchema);
    Require(result && result.value().kind == SchemaMatchKind::match,
            "the XRADIO reference fixture did not match");

    // Resource limits must be accepted and applied to every read made through this context.
    carta::zarr::OpenOptions options;
    options.cache_bytes = static_cast<std::size_t>(32U * 1024U * 1024U);
    options.io_threads = 2;
    options.decode_threads = 2;
    const auto context = carta::zarr::Context::Create(options);
    Require(static_cast<bool>(context), "Context::Create rejected valid resource limits");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), fixture.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed on reference fixture");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "OpenImage failed on reference fixture");

    const auto& desc = image.value().descriptor();
    const auto image_ids = ImageIds(dataset.value().descriptor().images);
    Require(image_ids ==
                std::vector<std::string>{"SKY", "MODEL", "RESIDUAL", "MASK_DECONVOLVE", "APERTURE", "COMPLEX"},
            "discovery did not enumerate the multi-image fixture in display order");
    Require(dataset.value().descriptor().default_image_id == "SKY", "default image was not the first openable image");
    // right_ascension and declination are float64 over (l, m) with no type attribute. Matching only
    // "has l and m" would list them as openable images; matching the whole axis set never reads them.
    for (const auto& coordinate : {"right_ascension", "declination", "velocity", "beam_params_label"}) {
        Require(std::find(image_ids.begin(), image_ids.end(), coordinate) == image_ids.end(),
                std::string("optional coordinate ") + coordinate + " was listed as an image");
        Require(!dataset.value().OpenImage(coordinate),
                std::string("optional coordinate ") + coordinate + " was openable as an image");
    }
    Require(HasDiagnostic(dataset.value().descriptor().diagnostics, "unsupported_coordinate_plane"),
            "aperture-plane diagnostic was not reported");
    Require(HasDiagnostic(dataset.value().descriptor().diagnostics, "unsupported_data_type"),
            "complex-dtype diagnostic was not reported");
    Require(desc.direction.has_value(), "DirectionCoordinate missing in reference fixture");
    Require(desc.direction->projection_parameters == std::vector<double>{0.25, -0.5},
            "direction projection parameters were not preserved");
    Require(desc.direction->native_pole_direction.at(0) == 0.0 && desc.direction->native_pole_direction.at(1) == 90.0,
            "native pole direction was not preserved in degrees");
    Require(desc.direction->reference_pixel.at(0) == 2.0 && desc.direction->reference_pixel.at(1) == 3.0,
            "direction reference pixels were not located from the coordinate values");
    Require(desc.image_role == "sky", "image role was not read from the variable's type attribute");
    Require(desc.data_groups == std::vector<std::string>{"base", "deconvolution"},
            "data group references were not reported");
    Require(desc.has_pixel_mask, "the image's own flag attribute did not localize its pixel mask");
    Require(desc.spectral.has_value(), "SpectralCoordinate missing in reference fixture");
    Require(desc.spectral->channel_frequencies == std::vector<double>{1.4e9, 1.401e9, 1.403e9},
            "spectral channel table was not preserved");
    Require(!desc.spectral->reference_pixel.has_value() && !desc.spectral->reference_value.has_value() &&
                !desc.spectral->increment.has_value(),
            "nonuniform spectral coordinates incorrectly exposed a linear description");
    // Withholding the linear description is not enough on its own: the consumer has to know it must
    // build a tabular axis, so the reason is reported rather than left silent.
    Require(HasDiagnostic(desc.diagnostics, "nonuniform_axis"),
            "nonuniform spectral coordinates did not report why they carry no linear description");
    Require(desc.polarization.has_value(), "PolarizationCoordinate missing in reference fixture");
    Require(!desc.polarization->labels.empty(), "Polarization labels empty in reference fixture");
    Require(desc.temporal.has_value(), "TemporalCoordinate missing in reference fixture");
    Require(desc.temporal->values == std::vector<double>{1.6e9} && desc.temporal->unit == "s" &&
                desc.temporal->scale == "UTC" && desc.temporal->format == "UNIX",
            "time coordinate values or attributes were not preserved");

    // The generator writes SKY as unsharded zstd chunks of (1, 1, 1, 2, 5) in stored axis order,
    // which is (2, 5, 1, 1, 1) in the logical order the geometry reports.
    const auto& geometry = image.value().chunk_geometry();
    Require(!geometry.sharded, "reference fixture was reported as sharded");
    Require((geometry.shard_shape == geometry.chunk_shape), "an unsharded image's shard is its chunk");
    Require((geometry.chunk_shape == std::vector<std::uint64_t>{2, 5, 1, 1, 1}),
            "reference fixture chunk shape changed");
    Require(geometry.compressor == "zstd", "reference fixture compressor was not reported as zstd");

    const auto beams = image.value().ReadBeams();
    if (!beams) {
        throw std::runtime_error("ReadBeams failed on reference fixture: " + beams.error().message);
    }

    const auto model = dataset.value().OpenImage("MODEL");
    Require(static_cast<bool>(model), "MODEL image could not be opened");
    Require(model.value().descriptor().image_role == "model", "MODEL descriptor used the wrong image role");
    Require(model.value().descriptor().data_groups == std::vector<std::string>{"base"},
            "MODEL data group metadata was not isolated");
    Require(model.value().descriptor().has_pixel_mask, "MODEL did not use the unique shape-matching flag");

    const auto residual = dataset.value().OpenImage("RESIDUAL");
    Require(static_cast<bool>(residual), "RESIDUAL image could not be opened");
    const auto deconvolution_mask = dataset.value().OpenImage("MASK_DECONVOLVE");
    Require(static_cast<bool>(deconvolution_mask), "MASK_DECONVOLVE was incorrectly treated as a flag");
    const auto aperture = dataset.value().OpenImage("APERTURE");
    Require(!aperture && aperture.error().code == ErrorCode::unsupported_data_type,
            "aperture-plane variable was openable");
    const auto complex = dataset.value().OpenImage("COMPLEX");
    Require(!complex && complex.error().code == ErrorCode::unsupported_data_type, "complex variable was openable");
    const auto flag = dataset.value().OpenImage("MASK_0");
    Require(!flag && flag.error().code == ErrorCode::not_found, "flag variable was exposed as an image");
}

void TestCompatibilityFixture() {
    const std::filesystem::path fixture(CARTA_ZARR_LEGACY_FIXTURE);
    Require(std::filesystem::exists(fixture), "the compatibility fixture is missing");
    const auto result = carta::zarr::ProbeSchema(fixture.string(), carta::zarr::kXradioImageSchema);
    Require(result && result.value().kind == SchemaMatchKind::match,
            "the compatibility fixture did not match");
}

void TestDiscoveryIgnoresNameAllowlist(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "UNLISTED" / "zarr.json", SkyArray());
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for unlisted image");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "unlisted image dataset did not open");
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"SKY", "UNLISTED"},
            "image discovery still used the known-name allowlist");
    const auto image = dataset.value().OpenImage("UNLISTED");
    Require(static_cast<bool>(image), "unlisted sky-plane image was not openable");
}

void TestMetadataCache(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "MODEL" / "zarr.json", SkyArray());

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for metadata cache");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for metadata cache");
    Require(static_cast<bool>(dataset.value().OpenImage("SKY")),
            "initial image open failed while populating metadata cache");

    // A Store is a read-only view. Once metadata has been observed, later schema operations must
    // use the same snapshot instead of rereading a changed metadata file.
    Write(root / "MODEL" / "zarr.json", "{not valid json");
    const auto cached_model = dataset.value().OpenImage("MODEL");
    Require(static_cast<bool>(cached_model), "cached metadata was not reused after the file changed");

    // The node list is cached as well, so nodes created after the first discovery are not visible
    // through an already-open Dataset.
    Write(root / "NEW" / "zarr.json", SkyArray());
    const auto new_image = dataset.value().OpenImage("NEW");
    Require(!new_image && new_image.error().code == ErrorCode::not_found,
            "metadata cache did not preserve the discovered node list");
}

void TestDiscoveryDoesNotDescendIntoArrayChunks(const std::filesystem::path& root) {
    CreateValidStore(root);
    // A chunk directory may contain arbitrary files, including a misleading zarr.json. It must
    // not be treated as a child Zarr node once the parent has been identified as an array.
    Write(root / "SKY" / "c" / "0" / "zarr.json", SkyArray());

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for chunk traversal test");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for chunk traversal test");
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"SKY"},
            "discovery descended into an array's chunk directory");
}

void TestCoordinateCompletion(const std::filesystem::path& root) {
    CreateValidStore(root / "uniform");
    Write(root / "uniform" / "frequency" / "zarr.json",
          NumericArray("[3]", R"(["frequency"])", "float64", R"({"units":"Hz"})"));
    WriteDoubles(root / "uniform" / "frequency" / "c" / "0", {100.0, 102.0, 104.0});
    Write(root / "uniform" / "l" / "zarr.json", NumericArray("[4]", R"(["l"])", "float64"));
    WriteDoubles(root / "uniform" / "l" / "c" / "0", {-0.003, -0.002, 0.0, 0.001});
    Write(root / "uniform" / "m" / "zarr.json", NumericArray("[5]", R"(["m"])", "float64"));
    WriteDoubles(root / "uniform" / "m" / "c" / "0", {-0.004, -0.003, -0.002, -0.001, 0.0});

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for coordinate completion");
    const auto uniform_dataset = carta::zarr::Dataset::Open(context.value(), (root / "uniform").string());
    Require(static_cast<bool>(uniform_dataset), "uniform coordinate dataset did not open");
    const auto uniform_image = uniform_dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(uniform_image), "uniform coordinate image did not open");
    const auto& uniform_desc = uniform_image.value().descriptor();
    Require(uniform_desc.spectral->reference_pixel == 1.0 && uniform_desc.spectral->reference_value == 100.0 &&
                uniform_desc.spectral->increment == 2.0,
            "uniform spectral coordinates did not expose the linear description");
    Require(uniform_desc.direction->reference_pixel.at(0) == 3.0 && uniform_desc.direction->reference_pixel.at(1) == 5.0,
            "exact direction reference pixels were not found by index");

    CreateValidStore(root / "inexact");
    Write(root / "inexact" / "l" / "zarr.json", NumericArray("[4]", R"(["l"])", "float64"));
    WriteDoubles(root / "inexact" / "l" / "c" / "0", {-0.003, -0.002, -0.001, -0.0005});
    Write(root / "inexact" / "m" / "zarr.json", NumericArray("[5]", R"(["m"])", "float64"));
    WriteDoubles(root / "inexact" / "m" / "c" / "0", {-0.004, -0.003, -0.002, -0.001, -0.0005});
    const auto inexact_dataset = carta::zarr::Dataset::Open(context.value(), (root / "inexact").string());
    Require(static_cast<bool>(inexact_dataset), "inexact coordinate dataset did not open");
    const auto inexact_image = inexact_dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(inexact_image), "inexact coordinate image did not open");
    Require(HasDiagnostic(inexact_image.value().descriptor().diagnostics, "inexact_reference_pixel"),
            "inexact direction reference pixel did not produce a diagnostic");
    Require(inexact_image.value().descriptor().direction->reference_pixel.at(0) == 4.0,
            "inexact direction reference pixel did not use linear extrapolation");
}

// The probe accepts a coordinate stored in any real type, so opening one has to accept the same
// set. A float32 coordinate used to probe as a match and then fail to open, which told a consumer
// the dataset was supported and then refused it.
void TestNonDoubleCoordinates(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])", "float32", R"({"units":"Hz"})"));
    WriteFloats(root / "frequency" / "c" / "0", {100.0F, 102.0F, 104.0F});
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])", "float32"));
    WriteFloats(root / "l" / "c" / "0", {-0.003F, -0.002F, -0.001F, 0.0F});
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])", "float32"));
    WriteFloats(root / "m" / "c" / "0", {-0.004F, -0.003F, -0.002F, -0.001F, 0.0F});

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for float32 coordinates");
    const auto supported = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(supported && supported.value().kind == SchemaMatchKind::match,
            "a float32 coordinate was not probed as a supported image dataset");

    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "float32-coordinate dataset did not open");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image),
            "float32-coordinate image did not open" + (image ? std::string{} : ": " + image.error().message));
    const auto& descriptor = image.value().descriptor();
    Require(descriptor.spectral.has_value(), "float32 frequency coordinate produced no spectral description");
    Require(descriptor.spectral->reference_value == 100.0 && descriptor.spectral->increment == 2.0,
            "float32 frequency values were not converted to their double equivalents");
}

// What the listing says about an image reaches a consumer opening it.
//
// The rule itself -- an image whose own frequency axis is a different length from the dataset's
// frequency coordinate is not one this profile opens -- is checked against a store in memory, in
// tests/schema_profile_test.cc. What is checked here is the three steps between that rule and a
// consumer, which are not a restatement of it: the dataset keeps the listing the probe built, the
// entry in it says the variable will not open, and OpenImage refuses it with the reason rather than
// with something of its own.
void TestADisagreeingImageIsRefusedThroughTheDataset(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "MODEL" / "zarr.json",
          "{\"shape\":[1,7,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},\"attributes\":{\"units\":\"Jy/beam\"},"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for a disagreeing image");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "the dataset holding a disagreeing image did not open");

    const auto& images = dataset.value().descriptor().images;
    Require(ImageIds(images) == std::vector<std::string>{"SKY", "MODEL"},
            "the dataset dropped the disagreeing image instead of listing it with its reason");
    const auto listed = std::find_if(images.begin(), images.end(),
                                     [](const auto& image) { return image.id == "MODEL"; });
    Require(listed != images.end() && !listed->openable,
            "the dataset offered an image that opening would refuse");

    Require(static_cast<bool>(dataset.value().OpenImage("SKY")),
            "the image that agrees with the coordinates did not open");

    const auto model = dataset.value().OpenImage("MODEL");
    Require(!model && model.error().code == ErrorCode::invalid_metadata,
            "the reason the listing carried did not reach a consumer opening the image");
}

// The same through the facade, which is where a consumer meets it: a dataset holding one stray
// document that will not parse opens, says which node it could not read, and opens its image. It
// used to be a dataset that did not open at all.
void TestADatasetWithAnUnparseableNodeStillOpens(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "JUNK" / "zarr.json", "{not valid json");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for an unparseable node");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset),
            "a dataset holding one unparseable node did not open" +
                (dataset ? std::string{} : ": " + dataset.error().message));

    const auto& said = dataset.value().descriptor().diagnostics;
    const auto junk = std::find_if(said.begin(), said.end(), [](const carta::zarr::Diagnostic& diagnostic) {
        return diagnostic.code == "unrecognised_node" && diagnostic.node_path.find("JUNK") != std::string::npos;
    });
    Require(junk != said.end(), "the dataset opened without saying which node it could not read");
    Require(static_cast<bool>(dataset.value().OpenImage("SKY")), "the image beside an unparseable node did not open");
}

// A node name is a relative path, and "./MASK_0" names the same node as "MASK_0". Two rules decide
// that and they disagree: Store::NormalizeNodeName accepts a "." component, and
// Transport::ArrayDirectory refuses one. So an image declaring its flag with that spelling describes
// perfectly -- every metadata read goes through the first rule -- and then fails every masked pixel
// read, which goes through the second.
//
// An image that opens, reports a pixel mask, and cannot be read is precisely the outcome ADR 0004
// reopened itself to avoid: a store that can be listed and never opened is a worse answer than one
// that is refused.
void TestANodeNameSpelledWithADotIsTheSameNode(const std::filesystem::path& root) {
    const std::string dimensions = R"(["time","frequency","polarization","l","m"])";
    CreateValidStore(root);
    Write(root / "SKY" / "zarr.json",
          NumericArray("[1,3,2,4,5]", dimensions, "float32", R"({"units":"Jy/beam","flag":"./MASK_0"})"));
    // Written out rather than through NumericArray, which fills with 0: a boolean array's fill
    // value has to be a boolean, and TensorStore refuses the array outright otherwise.
    Write(root / "MASK_0" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"bool\","
          "\"chunk_grid\":{\"name\":\"regular\",\"configuration\":{\"chunk_shape\":[1,3,2,4,5]}},"
          "\"chunk_key_encoding\":{\"name\":\"default\",\"configuration\":{\"separator\":\"/\"}},"
          "\"fill_value\":true,\"codecs\":[{\"name\":\"bytes\"}],"
          "\"attributes\":{\"type\":\"flag\"},\"dimension_names\":" + dimensions +
              ",\"zarr_format\":3,\"node_type\":\"array\"}");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for the dotted node name");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "the dotted-flag dataset did not open");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "SKY did not open in the dotted-flag dataset");
    Require(image.value().descriptor().has_pixel_mask,
            "a flag spelled with a leading ./ was not accepted while describing the image");

    // One pixel on every axis, so this says nothing about order and only asks whether the read
    // reaches the store at all.
    carta::zarr::ReadRequest request;
    request.axes.assign(image.value().descriptor().axes.size(), carta::zarr::Range{0, 1, 1});
    std::vector<float> destination(1, 0.0F);
    const auto read = image.value().Read(request, {destination.data(), destination.size()});
    Require(static_cast<bool>(read),
            std::string("an image whose flag is spelled ./MASK_0 opened and then could not be read: ") +
                (read ? "" : read.error().message));
}


// A sharding codec whose inner chunk shape describes nothing leaves the array unreadable, and this
// pins where that is noticed: at the parse, not after it.
//
// It used to be noticed after. ParseStorageLayout refused such an array, DescribeImage dropped the
// refusal, and the image opened reporting a chunk geometry synthesised from its own shape -- one
// chunk covering the whole image, unsharded, uncompressed, with no diagnostic. That geometry is
// what the reduction plans its blocks from, so a walk would have taken the cube as a single chunk.
//
// Now the array does not parse, so it never becomes an image. What is asserted here is that this
// does not happen in silence: the store reports no images, and the reason names the node.
void TestAShardingCodecThatDescribesNoChunks(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,3,2,4,5]}},\"attributes\":{\"units\":\"Jy/beam\"},"
          "\"codecs\":[{\"name\":\"sharding_indexed\",\"configuration\":{\"chunk_shape\":[1,1,1,0,5],"
          "\"codecs\":[{\"name\":\"bytes\"},{\"name\":\"zstd\"}]}}],"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");
    Write(root / "time" / "zarr.json", NumericArray("[1]", R"(["time"])"));
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])"));

    const auto probe = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(static_cast<bool>(probe), "probing the store failed outright");
    Require(probe.value().kind != SchemaMatchKind::match,
            "a store whose only image cannot be parsed was still matched");

    const auto said = std::find_if(probe.value().diagnostics.begin(), probe.value().diagnostics.end(),
                                   [](const carta::zarr::Diagnostic& diagnostic) {
                                       return diagnostic.code == "unreadable_array";
                                   });
    Require(said != probe.value().diagnostics.end(),
            "the store lost its image without saying which node or why");
    Require(said->node_path.find("SKY") != std::string::npos,
            "the diagnostic did not name the node that was skipped");
}
// A sharded array grids its store by shard; the inner chunk shape lives in the sharding codec.
// An image dataset without SKY is valid: discovery identifies the dataset from its image variables.
// An image that names a beam table is an image with a beam. When the labels that say which
// parameter is which could not be read, the parameters were simply not found and the answer was an
// empty beam list -- indistinguishable from an image that carries no beam at all, which a consumer
// acts on differently: CARTA reports no beam rather than a beam it failed to read.
void TestBeamTableWithUnreadableLabels(const std::filesystem::path& root) {
    const auto sky_with_beam =
        "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
        "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
        "\"attributes\":{\"units\":\"Jy/beam\",\"beam_fit_params\":\"BEAM\"},"
        "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
        "\"zarr_format\":3,\"node_type\":\"array\"}";

    const auto open_sky = [&](const std::filesystem::path& where) {
        const auto context = carta::zarr::Context::Create();
        Require(static_cast<bool>(context), "Context::Create failed for the beam label test");
        const auto dataset = carta::zarr::Dataset::Open(context.value(), where.string());
        Require(static_cast<bool>(dataset), "the beam label dataset did not open");
        const auto image = dataset.value().OpenImage("SKY");
        Require(static_cast<bool>(image), "SKY did not open in the beam label dataset");
        return image.value().ReadBeams();
    };

    // No label array at all, though the image says it has a beam table.
    CreateValidStore(root / "absent");
    Write(root / "absent" / "SKY" / "zarr.json", sky_with_beam);
    Write(root / "absent" / "BEAM" / "zarr.json",
          NumericArray("[1,3,2,3]", R"(["time","frequency","polarization","beam_params_label"])", "float64",
                       R"({"units":"rad"})"));
    const auto absent = open_sky(root / "absent");
    Require(!absent, "a beam table whose labels are missing was reported as an image with no beam");

    // Labels that read, and name something other than the three parameters a beam is made of.
    CreateValidStore(root / "unnamed");
    Write(root / "unnamed" / "SKY" / "zarr.json", sky_with_beam);
    Write(root / "unnamed" / "BEAM" / "zarr.json",
          NumericArray("[1,3,2,3]", R"(["time","frequency","polarization","beam_params_label"])", "float64",
                       R"({"units":"rad"})"));
    Write(root / "unnamed" / "beam_params_label" / "zarr.json",
          R"({"shape":[3],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":24}},
              "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[3]}},"attributes":{},
              "dimension_names":["beam_params_label"],"zarr_format":3,"node_type":"array"})");
    WriteUtf32(root / "unnamed" / "beam_params_label" / "c" / "0", {"one", "two", "three"}, 6);
    const auto unnamed = open_sky(root / "unnamed");
    Require(!unnamed && unnamed.error().code == ErrorCode::invalid_metadata,
            "beam labels naming no beam parameter were reported as an image with no beam");
}

void TestImageDatasetWithoutSky(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "RESIDUAL" / "zarr.json", SkyArray());
    Write(root / "time" / "zarr.json", NumericArray("[1]", R"(["time"])"));
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])"));

    const auto matched = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(matched && matched.value().kind == SchemaMatchKind::match,
            "an image dataset without SKY was not recognized");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for the SKY-less dataset");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for the SKY-less dataset");
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"RESIDUAL"},
            "the SKY-less dataset did not enumerate its only image");
    Require(static_cast<bool>(dataset.value().OpenImage("RESIDUAL")),
            "the only image of a SKY-less dataset could not be opened");
}

// Every image dataset must carry a coordinate array for every axis its image uses, time included:
// both XRADIO readers always write one. A missing coordinate is malformed metadata, and must be
// reported as such rather than as "not this schema".
void TestImageDatasetMissingTimeCoordinate(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "SKY" / "zarr.json", SkyArray());
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])"));

    const auto matched = carta::zarr::ProbeSchema(root.string(), carta::zarr::kXradioImageSchema);
    Require(matched && matched.value().kind == SchemaMatchKind::invalid,
            "an image dataset missing its time coordinate was not reported as an invalid match");
}

void TestShardedStorageLayout(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,3,2,4,5]}},\"attributes\":{\"units\":\"Jy/beam\"},"
          "\"codecs\":[{\"name\":\"sharding_indexed\",\"configuration\":{\"chunk_shape\":[1,1,1,2,5],"
          "\"codecs\":[{\"name\":\"bytes\"},{\"name\":\"blosc\"}]}}],"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");
    Write(root / "time" / "zarr.json", NumericArray("[1]", R"(["time"])"));
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])"));

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for the sharded store");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for the sharded store");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "OpenImage failed for the sharded store");

    // In logical order: stored (time, frequency, polarization, l, m) becomes (l, m, frequency,
    // polarization, time).
    const auto& geometry = image.value().chunk_geometry();
    Require(geometry.sharded, "sharded store was not reported as sharded");
    Require((geometry.shard_shape == std::vector<std::uint64_t>{4, 5, 3, 2, 1}), "shard shape was not reported");
    Require((geometry.chunk_shape == std::vector<std::uint64_t>{2, 5, 1, 1, 1}),
            "inner chunk shape was not taken from the sharding codec");
    Require(geometry.compressor == "blosc", "compressor inside the sharding codec was not reported");
}

}  // namespace

int main() {
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto root = std::filesystem::temp_directory_path() / ("carta-zarr-schema-test-" + std::to_string(suffix));
    try {
        std::filesystem::remove_all(root);
        TestValidAndTimeAxis(root / "valid");
        TestTimeGreaterThanOne(root / "time-two");
        TestNonMatchAndInvalid(root / "classification");
        TestMissingAndUnsupported(root);
        TestOpenSaysWhyItRefused(root / "refusal");
        TestCoordinateCompletion(root / "coordinates");
        TestNonDoubleCoordinates(root / "float32-coordinates");
        TestADisagreeingImageIsRefusedThroughTheDataset(root / "coordinate-disagreement");
        TestADatasetWithAnUnparseableNodeStillOpens(root / "unparseable-node");
        TestDiscoveryIgnoresNameAllowlist(root / "unlisted-image");
        TestMetadataCache(root / "metadata-cache");
        TestDiscoveryDoesNotDescendIntoArrayChunks(root / "array-chunks");
        TestImageDatasetWithoutSky(root / "image-no-sky");
        TestImageDatasetMissingTimeCoordinate(root / "image-no-time");
        TestShardedStorageLayout(root / "sharded");
        TestAShardingCodecThatDescribesNoChunks(root / "unusable-shard");
        TestEmptyAndOversizedMetadata(root / "metadata-bytes");
        TestBeamTableWithUnreadableLabels(root / "beam-labels");
        TestANodeNameSpelledWithADotIsTheSameNode(root / "dotted-node");
        TestReferenceFixture();
        TestCompatibilityFixture();
        std::filesystem::remove_all(root);
        std::cout << "carta-zarr schema probe tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr schema probe tests failed: " << error.what() << '\n';
        std::filesystem::remove_all(root);
        return 1;
    }
}
