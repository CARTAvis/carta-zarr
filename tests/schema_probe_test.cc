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

namespace {

using carta::zarr::ErrorCode;
using carta::zarr::ProbeKind;
using carta::zarr::SchemaMatchKind;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

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

    const auto is_xradio = carta::zarr::IsXradioImage(root.string());
    Require(is_xradio && is_xradio.value(), "IsXradioImage rejected a valid store");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed");
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"SKY"}, "unexpected image ids");
    Require(dataset.value().descriptor().default_image_id == "SKY", "unexpected default image id");

    const auto logical_size = dataset.value().Size(std::chrono::milliseconds(0));
    Require(logical_size && logical_size.value().bytes == 592 && logical_size.value().is_upper_bound,
            "logical Zarr size calculation was incorrect");

    const auto physical_size = dataset.value().Size(std::chrono::milliseconds(5000));
    Require(physical_size && physical_size.value().bytes > 0 && !physical_size.value().is_upper_bound,
            "physical Zarr size calculation was not used");

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
    const auto result = carta::zarr::IsXradioImage(root.string());
    Require(result && result.value(), "time > 1 must remain valid at library level");
    const auto context = carta::zarr::Context::Create();
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open rejected time > 1");
    const auto image = dataset.value().OpenImage("SKY");
    Require(image && image.value().descriptor().axes.back().length == 2, "time axis length was not preserved");
}

void TestNonMatchAndInvalid(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "OTHER" / "zarr.json", NumericArray("[2]", "[\"x\"]"));
    const auto non_match = carta::zarr::IsXradioImage(root.string());
    Require(non_match && !non_match.value(), "valid non-XRADIO Zarr was not a non-match");
    Require(carta::zarr::Probe(root.string()).kind == ProbeKind::zarr_without_supported_schema,
            "Probe did not report a valid unsupported schema");

    CreateValidStore(root);
    std::filesystem::remove(root / "time" / "zarr.json");
    const auto invalid = carta::zarr::IsXradioImage(root.string());
    Require(!invalid && invalid.error().code == ErrorCode::invalid_metadata,
            "metadata-incomplete XRADIO-like store did not report invalid metadata");
    Require(carta::zarr::Probe(root.string()).kind == ProbeKind::invalid_dataset,
            "Probe did not report invalid XRADIO-like metadata");
}

void TestMissingAndUnsupported(const std::filesystem::path& root) {
    const auto missing = carta::zarr::IsXradioImage((root / "missing").string());
    Require(!missing && missing.error().code == ErrorCode::not_found, "missing store error category changed");

    Write(root / "not-zarr" / "zarr.json", "{\"zarr_format\": 2, \"node_type\": \"group\"}");
    const auto unsupported = carta::zarr::IsXradioImage((root / "not-zarr").string());
    Require(!unsupported && unsupported.error().code == ErrorCode::unsupported_zarr_version,
            "unsupported Zarr version error category changed");

    const auto unknown_schema = carta::zarr::ProbeSchema(root.string(), "future.schema");
    Require(!unknown_schema && unknown_schema.error().code == ErrorCode::unsupported_schema,
            "unknown schema error category changed");
}

void TestReferenceFixture() {
    const std::filesystem::path fixture(CARTA_ZARR_REFERENCE_FIXTURE);
    Require(std::filesystem::exists(fixture), "the XRADIO reference fixture is missing from tests/data");

    const auto result = carta::zarr::IsXradioImage(fixture.string());
    Require(result && result.value(), "the XRADIO reference fixture did not match");

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
    Require(dataset.value().descriptor().default_image_id == "SKY", "default image was not the first readable image");
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

    // The generator writes SKY as unsharded zstd chunks of (1, 1, 1, 2, 5) in stored axis order.
    Require(desc.storage.has_value(), "StorageLayout missing in reference fixture");
    Require(!desc.storage->sharded, "reference fixture was reported as sharded");
    Require(desc.storage->shard_shape.empty(), "unsharded reference fixture reported a shard shape");
    Require((desc.storage->chunk_shape == std::vector<std::uint64_t>{1, 1, 1, 2, 5}),
            "reference fixture chunk shape changed");
    Require(desc.storage->compressor == "zstd", "reference fixture compressor was not reported as zstd");

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
    const auto result = carta::zarr::IsXradioImage(fixture.string());
    Require(result && result.value(), "the compatibility fixture did not match");
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
    const auto supported = carta::zarr::IsXradioImage(root.string());
    Require(supported && supported.value(), "a float32 coordinate was not probed as a supported image dataset");

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

// A declared flag is trusted by every read: reads apply the pixel mask by default, and the mask is
// read with the selection built for the image. A flag that is missing, not boolean, or shaped
// differently used to be accepted here and rejected on every read -- or, for a numeric array of the
// right shape, converted to bool and applied as if nonzero meant valid.
void TestDeclaredPixelMaskIsValidated(const std::filesystem::path& root) {
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for declared pixel masks");

    const std::string sky_with_flag =
        "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
        "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
        "\"attributes\":{\"units\":\"Jy/beam\",\"flag\":\"MASK_0\"},"
        "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
        "\"zarr_format\":3,\"node_type\":\"array\"}";
    const std::string mask_dimensions = R"(["time","frequency","polarization","l","m"])";

    // The mask a well-formed image declares: boolean, and the image's own dimensions in order.
    CreateValidStore(root / "good");
    Write(root / "good" / "SKY" / "zarr.json", sky_with_flag);
    Write(root / "good" / "MASK_0" / "zarr.json",
          NumericArray("[1,3,2,4,5]", mask_dimensions, "bool", R"({"type":"flag"})"));
    const auto good = carta::zarr::Dataset::Open(context.value(), (root / "good").string());
    Require(static_cast<bool>(good), "the declared-mask dataset did not open");
    const auto good_image = good.value().OpenImage("SKY");
    Require(static_cast<bool>(good_image), "an image declaring a well-formed mask did not open");
    Require(good_image.value().descriptor().has_pixel_mask, "a well-formed declared mask was not selected");
    Require(good_image.value().descriptor().pixel_mask_id == "MASK_0", "the declared mask was not the one selected");

    const auto open_with_mask = [&](const std::string& name, const std::string& mask_metadata) {
        CreateValidStore(root / name);
        Write(root / name / "SKY" / "zarr.json", sky_with_flag);
        if (!mask_metadata.empty()) {
            Write(root / name / "MASK_0" / "zarr.json", mask_metadata);
        }
        const auto dataset = carta::zarr::Dataset::Open(context.value(), (root / name).string());
        Require(static_cast<bool>(dataset), "the declared-mask dataset " + name + " did not open");
        return dataset.value().OpenImage("SKY");
    };

    const auto absent = open_with_mask("absent", {});
    Require(!absent, "an image declaring a flag variable that does not exist was opened");

    const auto numeric = open_with_mask(
        "numeric", NumericArray("[1,3,2,4,5]", mask_dimensions, "float32", R"({"type":"flag"})"));
    Require(!numeric && numeric.error().code == ErrorCode::unsupported_data_type,
            "a numeric variable was accepted as a pixel mask");

    const auto unmarked = open_with_mask("unmarked", NumericArray("[1,3,2,4,5]", mask_dimensions, "bool"));
    Require(!unmarked && unmarked.error().code == ErrorCode::invalid_metadata,
            "a boolean variable that is not a flag was accepted as a pixel mask");

    const auto reshaped =
        open_with_mask("reshaped", NumericArray("[1,3,2,4,4]", mask_dimensions, "bool", R"({"type":"flag"})"));
    Require(!reshaped && reshaped.error().code == ErrorCode::invalid_metadata,
            "a flag whose shape differs from the image was accepted as a pixel mask");

    const auto transposed = open_with_mask(
        "transposed", NumericArray("[1,3,2,5,4]", R"(["time","frequency","polarization","m","l"])", "bool",
                                   R"({"type":"flag"})"));
    Require(!transposed && transposed.error().code == ErrorCode::invalid_metadata,
            "a flag whose dimension order differs from the image was accepted as a pixel mask");
}

// Optional metadata comes from a file, so its shape is whatever was written, not whatever the
// schema describes. A telescope position holding a string where a number belongs used to throw out
// of nlohmann and past the Result the caller is holding. It is a value the image can do without:
// the image opens, and only the position it could not read is missing.
void TestMalformedObservationMetadata(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
          "\"attributes\":{\"units\":\"Jy/beam\",\"object_name\":\"Zarr test source\","
          "\"telescope\":{\"name\":\"Test scope\",\"direction\":{\"data\":[\"north\",0.0]},"
          "\"distance\":{\"data\":[6371000.0]}}},"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for malformed observation metadata");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "the malformed-observation dataset did not open");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "a malformed telescope position closed an otherwise readable image");

    const auto& observation = image.value().descriptor().observation;
    Require(observation.has_value(), "the image reported no observation metadata at all");
    Require(observation->telescope_name == "Test scope", "the readable telescope metadata was dropped as well");
    Require(!observation->observatory_position.has_value(),
            "a telescope position holding a string was converted rather than skipped");
}

// Coordinates belong to the dataset, not to one image, and discovery lists an image on its
// dimension names alone. A second image whose own frequency axis is a different length was listed
// as readable and then described with the dataset's frequency coordinate as though it were its
// own: an image reporting three channels' worth of coordinates over seven channels of pixels.
void TestImageDisagreeingWithACoordinate(const std::filesystem::path& root) {
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
    Require(ImageIds(dataset.value().descriptor().images) == std::vector<std::string>{"SKY", "MODEL"},
            "the dataset did not list both images");

    const auto sky = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(sky), "the image that agrees with the coordinates did not open");

    const auto model = dataset.value().OpenImage("MODEL");
    Require(!model && model.error().code == ErrorCode::invalid_metadata,
            "an image whose frequency axis disagrees with the frequency coordinate was described anyway");
}

void TestAmbiguousPixelMask(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "MODEL" / "zarr.json", SkyArray());
    Write(root / "FLAG_1" / "zarr.json",
          NumericArray("[1,3,2,4,5]", R"(["time","frequency","polarization","l","m"])", "bool", R"({"type":"flag"})"));
    Write(root / "FLAG_2" / "zarr.json",
          NumericArray("[1,3,2,4,5]", R"(["time","frequency","polarization","l","m"])", "bool", R"({"type":"flag"})"));
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for ambiguous mask");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "ambiguous-mask dataset did not open");
    const auto model = dataset.value().OpenImage("MODEL");
    Require(static_cast<bool>(model), "MODEL image did not open for ambiguous mask");
    Require(!model.value().descriptor().has_pixel_mask, "ambiguous flags selected a pixel mask");
    Require(HasDiagnostic(model.value().descriptor().diagnostics, "ambiguous_pixel_mask"),
            "ambiguous flags did not produce a diagnostic");
}

// The beam table is a four-dimensional array read through a flat buffer, and until this test the
// only thing pinning that indexing was "ReadBeams did not fail". The conformance fixture cannot
// pin it either: every one of its six (frequency, polarization) planes holds the same triple, so
// transposing two strides would go unnoticed. This store gives every plane and parameter a value
// that identifies it: major = 100*frequency + 10*polarization, and minor and the position angle
// one and two above it.
void TestBeamTableIndexing(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
          "\"attributes\":{\"units\":\"Jy/beam\",\"beam_fit_params\":\"BEAM\"},"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");

    // (time, frequency, polarization, beam_params_label) = (1, 3, 2, 3), C order.
    Write(root / "BEAM" / "zarr.json",
          NumericArray("[1,3,2,3]", R"(["time","frequency","polarization","beam_params_label"])", "float64",
                       R"({"units":"rad"})"));
    std::vector<double> beam_values;
    for (std::uint64_t frequency = 0; frequency < 3; ++frequency) {
        for (std::uint64_t polarization = 0; polarization < 2; ++polarization) {
            const double base = (100.0 * static_cast<double>(frequency)) + (10.0 * static_cast<double>(polarization));
            beam_values.push_back(base);
            beam_values.push_back(base + 1.0);
            beam_values.push_back(base + 2.0);
        }
    }
    WriteDoubles(root / "BEAM" / "c" / "0" / "0" / "0" / "0", beam_values);

    // Labels are deliberately not in major/minor/pa order: the reader must locate a parameter by its
    // label, not by its position.
    Write(root / "beam_params_label" / "zarr.json",
          R"({"shape":[3],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":24}},
              "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[3]}},"attributes":{},
              "dimension_names":["beam_params_label"],"zarr_format":3,"node_type":"array"})");
    WriteUtf32(root / "beam_params_label" / "c" / "0", {"minor", "major", "pa"}, 6);

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed for the beam table");
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for the beam table");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "OpenImage failed for the beam table");

    const auto beams = image.value().ReadBeams();
    Require(static_cast<bool>(beams),
            "ReadBeams failed on the beam table" + (beams ? std::string{} : ": " + beams.error().message));
    Require(beams.value().size() == 6, "the beam table did not decode one beam per frequency and polarization");

    for (const auto& beam : beams.value()) {
        const double base =
            (100.0 * static_cast<double>(beam.channel)) + (10.0 * static_cast<double>(beam.polarization));
        // Labels are stored as minor, major, pa, so major sits one past the plane's base value.
        Require(beam.major == base + 1.0, "beam major was read from the wrong element");
        Require(beam.minor == base, "beam minor was read from the wrong element");
        Require(beam.position_angle == base + 2.0, "beam position angle was read from the wrong element");
        Require(beam.unit == "rad", "beam unit was not read from the beam array attributes");
    }
}

// A beam table with more than one time plane used to be read as its first plane only, silently.
// The library reports the whole time axis and leaves any selection to the consumer, and beams now
// follow that too. Time varies slowest, so a single-plane table is unaffected.
void TestBeamTableTimePlanes(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
          "\"attributes\":{\"units\":\"Jy/beam\",\"beam_fit_params\":\"BEAM\"},"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");

    // (time, frequency, polarization, beam_params_label) = (2, 3, 2, 3), C order.
    Write(root / "BEAM" / "zarr.json",
          NumericArray("[2,3,2,3]", R"(["time","frequency","polarization","beam_params_label"])", "float64",
                       R"({"units":"rad"})"));
    std::vector<double> beam_values;
    for (std::uint64_t time = 0; time < 2; ++time) {
        for (std::uint64_t frequency = 0; frequency < 3; ++frequency) {
            for (std::uint64_t polarization = 0; polarization < 2; ++polarization) {
                const double base = (1000.0 * static_cast<double>(time)) + (100.0 * static_cast<double>(frequency)) +
                                    (10.0 * static_cast<double>(polarization));
                beam_values.push_back(base);
                beam_values.push_back(base + 1.0);
                beam_values.push_back(base + 2.0);
            }
        }
    }
    WriteDoubles(root / "BEAM" / "c" / "0" / "0" / "0" / "0", beam_values);
    Write(root / "beam_params_label" / "zarr.json",
          R"({"shape":[3],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":24}},
              "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[3]}},"attributes":{},
              "dimension_names":["beam_params_label"],"zarr_format":3,"node_type":"array"})");
    WriteUtf32(root / "beam_params_label" / "c" / "0", {"minor", "major", "pa"}, 6);

    const auto context = carta::zarr::Context::Create();
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for the multi-plane beam table");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "OpenImage failed for the multi-plane beam table");

    const auto beams = image.value().ReadBeams();
    Require(static_cast<bool>(beams), "ReadBeams failed on the multi-plane beam table");
    Require(beams.value().size() == 12, "the multi-plane beam table did not report every plane");

    for (const auto& beam : beams.value()) {
        const double base = (1000.0 * static_cast<double>(beam.time)) + (100.0 * static_cast<double>(beam.channel)) +
                            (10.0 * static_cast<double>(beam.polarization));
        Require(beam.major == base + 1.0, "beam major was read from the wrong time plane");
        Require(beam.minor == base, "beam minor was read from the wrong time plane");
        Require(beam.position_angle == base + 2.0, "beam position angle was read from the wrong time plane");
    }

    // The first plane still reads back exactly as it did when it was all that was reported.
    Require(beams.value().front().time == 0 && beams.value().front().channel == 0 &&
                beams.value().front().polarization == 0,
            "time did not vary slowest, so single-plane callers would see a different order");
}

// A beam table need not carry a time dimension; ReadBeams treats an absent one as a single
// implicit plane. Addressing the array must not insist on naming a dimension the array lacks.
void TestBeamTableWithoutTimeDimension(const std::filesystem::path& root) {
    CreateValidStore(root);
    Write(root / "SKY" / "zarr.json",
          "{\"shape\":[1,3,2,4,5],\"data_type\":\"float32\",\"chunk_grid\":{\"name\":\"regular\","
          "\"configuration\":{\"chunk_shape\":[1,1,1,2,5]}},"
          "\"attributes\":{\"units\":\"Jy/beam\",\"beam_fit_params\":\"BEAM\"},"
          "\"dimension_names\":[\"time\",\"frequency\",\"polarization\",\"l\",\"m\"],"
          "\"zarr_format\":3,\"node_type\":\"array\"}");

    // (frequency, polarization, beam_params_label) = (3, 2, 3) -- no time dimension at all.
    Write(root / "BEAM" / "zarr.json", NumericArray("[3,2,3]", R"(["frequency","polarization","beam_params_label"])",
                                                    "float64", R"({"units":"rad"})"));
    std::vector<double> beam_values;
    for (std::uint64_t frequency = 0; frequency < 3; ++frequency) {
        for (std::uint64_t polarization = 0; polarization < 2; ++polarization) {
            const double base = (100.0 * static_cast<double>(frequency)) + (10.0 * static_cast<double>(polarization));
            beam_values.push_back(base);
            beam_values.push_back(base + 1.0);
            beam_values.push_back(base + 2.0);
        }
    }
    WriteDoubles(root / "BEAM" / "c" / "0" / "0" / "0", beam_values);
    Write(root / "beam_params_label" / "zarr.json",
          R"({"shape":[3],"data_type":{"name":"fixed_length_utf32","configuration":{"length_bytes":24}},
              "chunk_grid":{"name":"regular","configuration":{"chunk_shape":[3]}},"attributes":{},
              "dimension_names":["beam_params_label"],"zarr_format":3,"node_type":"array"})");
    WriteUtf32(root / "beam_params_label" / "c" / "0", {"minor", "major", "pa"}, 6);

    const auto context = carta::zarr::Context::Create();
    const auto dataset = carta::zarr::Dataset::Open(context.value(), root.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed for the time-less beam table");
    const auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "OpenImage failed for the time-less beam table");

    const auto beams = image.value().ReadBeams();
    Require(static_cast<bool>(beams), "a beam table without a time dimension was not readable" +
                                          (beams ? std::string{} : ": " + beams.error().message));
    Require(beams.value().size() == 6, "the time-less beam table did not report one beam per plane");
    for (const auto& beam : beams.value()) {
        const double base =
            (100.0 * static_cast<double>(beam.channel)) + (10.0 * static_cast<double>(beam.polarization));
        Require(beam.time == 0, "a beam table without a time dimension reported a nonzero time index");
        Require(beam.major == base + 1.0, "beam major was misaddressed without a time dimension");
        Require(beam.minor == base, "beam minor was misaddressed without a time dimension");
        Require(beam.position_angle == base + 2.0, "beam position angle was misaddressed without a time dimension");
    }
}

// A sharded array grids its store by shard; the inner chunk shape lives in the sharding codec.
// An image dataset without SKY is valid: discovery identifies the dataset from its image variables.
void TestImageDatasetWithoutSky(const std::filesystem::path& root) {
    Write(root / "zarr.json", RootMetadata());
    Write(root / "RESIDUAL" / "zarr.json", SkyArray());
    Write(root / "time" / "zarr.json", NumericArray("[1]", R"(["time"])"));
    Write(root / "frequency" / "zarr.json", NumericArray("[3]", R"(["frequency"])"));
    Write(root / "polarization" / "zarr.json", PolarizationArray());
    Write(root / "l" / "zarr.json", NumericArray("[4]", R"(["l"])"));
    Write(root / "m" / "zarr.json", NumericArray("[5]", R"(["m"])"));

    const auto matched = carta::zarr::IsXradioImage(root.string());
    Require(matched && matched.value(), "an image dataset without SKY was not recognized");

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

    const auto matched = carta::zarr::IsXradioImage(root.string());
    Require(!matched && matched.error().code == ErrorCode::invalid_metadata,
            "an image dataset missing its time coordinate was not reported as invalid metadata");
    Require(carta::zarr::Probe(root.string()).kind == ProbeKind::invalid_dataset,
            "Probe did not report the missing time coordinate as an invalid dataset");
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

    const auto& storage = image.value().descriptor().storage;
    Require(storage.has_value(), "sharded store reported no StorageLayout");
    Require(storage->sharded, "sharded store was not reported as sharded");
    Require((storage->shard_shape == std::vector<std::uint64_t>{1, 3, 2, 4, 5}), "shard shape was not reported");
    Require((storage->chunk_shape == std::vector<std::uint64_t>{1, 1, 1, 2, 5}),
            "inner chunk shape was not taken from the sharding codec");
    Require(storage->compressor == "blosc", "compressor inside the sharding codec was not reported");
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
        TestCoordinateCompletion(root / "coordinates");
        TestNonDoubleCoordinates(root / "float32-coordinates");
        TestAmbiguousPixelMask(root / "ambiguous-mask");
        TestImageDisagreeingWithACoordinate(root / "coordinate-disagreement");
        TestMalformedObservationMetadata(root / "observation");
        TestDeclaredPixelMaskIsValidated(root / "declared-mask");
        TestDiscoveryIgnoresNameAllowlist(root / "unlisted-image");
        TestMetadataCache(root / "metadata-cache");
        TestDiscoveryDoesNotDescendIntoArrayChunks(root / "array-chunks");
        TestImageDatasetWithoutSky(root / "image-no-sky");
        TestImageDatasetMissingTimeCoordinate(root / "image-no-time");
        TestShardedStorageLayout(root / "sharded");
        TestBeamTableIndexing(root / "beam-table");
        TestBeamTableTimePlanes(root / "beam-time");
        TestBeamTableWithoutTimeDimension(root / "beam-notime");
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
