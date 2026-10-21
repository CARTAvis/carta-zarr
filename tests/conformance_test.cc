/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// Conformance: assert what this library believes about XRADIO against a store that XRADIO itself
// produced. The fixtures come from tests/data/generate_conformance_fixtures.py, which converts a
// FITS image, and two CASA images with internal masks, with a pinned XRADIO. When that pin is bumped and one of these
// fails, the failure names the assumption the new XRADIO version broke.

#include "support/check.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <carta-zarr/carta_zarr.h>
#include <carta-zarr/descriptor.h>
#include <carta-zarr/read.h>
#include <carta-zarr/reduce.h>

namespace {

using carta::zarr::AxisRole;

const std::filesystem::path kFixture{CARTA_ZARR_CONFORMANCE_FIXTURE};
const std::filesystem::path kFlaggedFixture{CARTA_ZARR_CONFORMANCE_FLAGGED_FIXTURE};
// The schema version the fixtures name.
const std::string kSchemaVersion{CARTA_ZARR_CONFORMANCE_SCHEMA_VERSION};

using carta::zarr::testing::Require;

void RequireClose(double actual, double expected, double tolerance, const std::string& message) {
    Require(std::abs(actual - expected) <= tolerance,
            message + " (expected " + std::to_string(expected) + ", got " + std::to_string(actual) + ")");
}

carta::zarr::Image OpenSky() {
    Require(std::filesystem::exists(kFixture),
            "the conformance fixture is missing; run tests/data/generate_conformance_fixtures.py");

    const auto matched = carta::zarr::ProbeSchema(kFixture.string(), carta::zarr::kXradioImageSchema);
    Require(matched && matched.value().kind == carta::zarr::SchemaMatchKind::match,
            "XRADIO's own output was not recognized as an image dataset");
    Require(matched.value().schema_version == kSchemaVersion,
            "XRADIO's own output reported schema version '" + matched.value().schema_version + "'");

    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    auto dataset = carta::zarr::Dataset::Open(context.value(), kFixture.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed on XRADIO's own output");

    // XRADIO writes right_ascension, declination and velocity as coordinates alongside SKY. Only
    // SKY is typed as an image, so only SKY is one.
    Require(dataset.value().descriptor().images.size() == 1 &&
                dataset.value().descriptor().images.front().id == "SKY" &&
                dataset.value().descriptor().images.front().openable,
            "discovery did not report SKY as the only image of XRADIO's own output");
    auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "SKY could not be opened");
    return image.value();
}

void TestIdentityAndAxes(const carta::zarr::ImageDescriptor& sky) {
    // XRADIO 1.3.0 writes the BTYPE of the FITS image it converted.
    Require(sky.sub_type == "Intensity", "SKY's sub_type was '" + sky.sub_type + "'");
    Require(sky.unit == "Jy/beam", "BUNIT did not survive the conversion as the image unit");
    Require(sky.data_groups == std::vector<std::string>{"base"}, "SKY's data group membership changed");
    Require(sky.pixel_mask_id.empty(), "an image converted from a FITS with no mask reported a pixel mask");

    // The FITS was (RA 5, DEC 4, STOKES 3, FREQ 2); the XRADIO profile reports its own logical order.
    // That order is the profile's rather than a promise of the public API, which reaches an axis by
    // its role, but it is still what this profile does and what every consumer of it has seen.
    constexpr std::array<AxisRole, 5> expected_roles{AxisRole::spatial_x, AxisRole::spatial_y, AxisRole::spectral,
                                                     AxisRole::polarization, AxisRole::time};
    constexpr std::array<std::uint64_t, 5> expected_lengths{5, 4, 2, 3, 1};
    Require(sky.axes.size() == expected_roles.size(), "SKY did not report five axes");
    const auto* expected_role = expected_roles.begin();
    const auto* expected_length = expected_lengths.begin();
    std::size_t index = 0;
    for (const auto& axis : sky.axes) {
        Require(axis.role == *expected_role, "axis " + std::to_string(index) + " has the wrong role");
        Require(axis.length == *expected_length, "axis " + std::to_string(index) + " has the wrong length");
        ++expected_role;
        ++expected_length;
        ++index;
    }
    for (std::size_t role = 0; role < expected_roles.size(); ++role) {
        Require(carta::zarr::AxisIndex(sky.axes, expected_roles.at(role)) == role,
                "AxisIndex did not find the axis playing role " + std::to_string(role) + " where it is");
    }
    Require(!carta::zarr::AxisIndex(sky.axes, AxisRole::other), "AxisIndex found an axis no axis plays");
}

void TestDirection(const carta::zarr::ImageDescriptor& sky) {
    Require(sky.direction.has_value(), "no direction coordinate was reported");
    const auto& direction = *sky.direction;
    Require(direction.projection == "SIN", "projection did not survive the conversion");
    Require(direction.projection_parameters == std::vector<double>{0.0, 0.0}, "projection parameters were dropped");
    Require(direction.reference_frame == "FK5", "RADESYS was not normalized to the canonical reference frame");
    Require(direction.equinox.has_value(), "EQUINOX was dropped");
    RequireClose(*direction.equinox, 2000.0, 1.0e-9, "equinox");
    // LONPOLE 180 and LATPOLE 30 are stored in radians and reported in degrees.
    RequireClose(direction.native_pole_direction.at(0), 180.0, 1.0e-9, "native pole longitude");
    RequireClose(direction.native_pole_direction.at(1), 30.0, 1.0e-9, "native pole latitude");
    RequireClose(direction.reference_value.at(0), 45.0, 1.0e-9, "reference longitude");
    RequireClose(direction.reference_value.at(1), 30.0, 1.0e-9, "reference latitude");
}

// A FITS file with no OBSRA or OBSDEC: XRADIO writes its reference direction as the pointing centre,
// and that is what is reported. Nothing else of its header is left over for `user`.
void TestPointingCentreWithoutObsra(const carta::zarr::ImageDescriptor& sky) {
    Require(sky.observation && sky.observation->pointing_center, "no pointing centre was reported");
    const auto& pointing = *sky.observation->pointing_center;
    RequireClose(pointing.value.at(0), 45.0, 1.0e-9, "pointing longitude, from the reference direction");
    RequireClose(pointing.value.at(1), 30.0, 1.0e-9, "pointing latitude, from the reference direction");
    Require(pointing.reference_frame == "FK5", "the pointing centre's frame was not named as the reference frame is");
    Require(sky.user_keywords.empty(), "a FITS file with no uninterpreted cards reported user keywords");
}

void TestSpectralAndPolarization(const carta::zarr::ImageDescriptor& sky) {
    Require(sky.spectral.has_value(), "no spectral coordinate was reported");
    const auto& spectral = *sky.spectral;
    Require(spectral.values.size() == 2, "the channel table lost a channel");
    RequireClose(spectral.values.at(0), 1.4e9, 1.0, "first channel frequency");
    RequireClose(spectral.values.at(1), 1.401e9, 1.0, "second channel frequency");
    Require(spectral.unit == "Hz", "the frequency unit was not taken from frequency's units");
    Require(spectral.reference_frame == "LSRK", "SPECSYS was not normalized to the canonical spectral system");
    Require(spectral.rest_frequency.has_value(), "RESTFRQ was dropped");
    RequireClose(*spectral.rest_frequency, 1.420405751e9, 1.0, "rest frequency");
    // These channels are evenly spaced, so the optional linear description must be present.
    Require(spectral.increment.has_value(), "evenly spaced channels reported no increment");
    RequireClose(*spectral.increment, 1.0e6, 1.0, "channel increment");

    // The polarization labels are a blosc-compressed fixed_length_utf32 array, which TensorStore
    // cannot read: this is the custom string decoder running against XRADIO's real output.
    Require(sky.polarization.has_value(), "no polarization coordinate was reported");
    Require(sky.polarization->labels == std::vector<std::string>({"I", "Q", "U"}),
            "the compressed polarization label array did not decode to the FITS Stokes ordering");
}

void TestTemporalAndStorage(const carta::zarr::ImageDescriptor& sky, const carta::zarr::ChunkGeometry& geometry) {
    // DATE-OBS 2020-05-31T12:00:00 is MJD 59000.5. The FITS path writes MJD days, not the unix
    // seconds the schema document describes.
    Require(sky.temporal.has_value(), "no time coordinate was reported");
    Require(sky.temporal->values.size() == 1, "the time coordinate lost its only value");
    RequireClose(sky.temporal->values.at(0), 59000.5, 1.0e-6, "time value");
    Require(sky.temporal->unit == "d", "the time unit changed");
    Require(sky.temporal->scale == "UTC", "the time scale was not normalized");
    Require(sky.temporal->format == "MJD", "the time format was not normalized");

    Require(!geometry.sharded, "XRADIO's zarr writer started sharding");
    Require(geometry.compressor == "blosc", "XRADIO's zarr writer changed compressor");
    // Stored as (time, frequency, polarization, l, m) chunks of (1, 2, 3, 5, 4), reported in the
    // descriptor's logical order.
    Require((geometry.chunk_shape == std::vector<std::uint64_t>{5, 4, 2, 3, 1}),
            "the chunk shape changed, or is no longer reported in logical axis order");
}

// The flagged fixture is the same FITS image, made into two CASA images with internal masks and
// converted together, which is what XRADIO writes from CASA: each flag a Zarr bool, true where the
// pixel is bad, and tied to its sky image only through the root's `data_groups`.
//
// Every pixel of the FITS image is finite, and its value is f*60 + s*20 + y*5 + x for frequency f,
// Stokes s and casacore pixel (x, y), so a pixel's value names it. The generator's MASKED_PIXELS
// are the ones these expect gone.
struct MaskedPixel {
    std::uint64_t frequency;
    std::uint64_t polarization;
    float value;
};

constexpr float FitsValue(int x, int y, int stokes, int frequency) {
    return static_cast<float>((frequency * 60) + (stokes * 20) + (y * 5) + x);
}

const std::vector<MaskedPixel> kSkyMasked{{0, 0, FitsValue(1, 2, 0, 0)}, {1, 1, FitsValue(3, 0, 1, 1)}};
const std::vector<MaskedPixel> kOtherMasked{{0, 2, FitsValue(4, 3, 2, 0)}};

constexpr std::uint64_t kPlanePixels = 5 * 4;

std::vector<float> ReadPlane(const carta::zarr::Image& image, std::uint64_t frequency, std::uint64_t polarization,
                             bool masked) {
    carta::zarr::ReadRequest request;
    request.axes = {{0, 5, 1}, {0, 4, 1}, {frequency, 1, 1}, {polarization, 1, 1}, {0, 1, 1}};
    std::vector<float> pixels(kPlanePixels);
    carta::zarr::ReadOptions options;
    options.apply_pixel_mask = masked;
    const auto read = image.Read(request, {pixels.data(), pixels.size()}, options);
    Require(static_cast<bool>(read), "a plane of " + image.descriptor().id + " could not be read");
    return pixels;
}

// How many pixels of each channel the reduction counted over the whole plane, at one polarization.
std::vector<double> CountedPerChannel(const carta::zarr::Image& image, std::uint64_t polarization) {
    const carta::zarr::RegionMask whole{0, 0, 5, 4};
    carta::zarr::SpectralReduceRequest request;
    request.planes.spectral = {0, 2, 1};
    request.planes.polarization = polarization;
    request.regions = {&whole, 1};
    request.statistics = carta::zarr::StatisticSet{carta::zarr::Statistic::num_pixels};
    std::vector<double> counted(2, -1.0);
    const auto reduced = image.ReduceSpectral(request, [&](const carta::zarr::SpectralBlock& block) {
        if (block.complete) {
            for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
                counted.at(block.first_channel + channel) = block.Totals(0, channel).num_pixels;
            }
        }
        return true;
    });
    Require(static_cast<bool>(reduced), "a reduction of " + image.descriptor().id + " failed");
    return counted;
}

void TestFlagsMaskTheirOwnImage(const carta::zarr::Dataset& dataset, const std::string& image_id,
                                const std::string& flag_id, const std::vector<MaskedPixel>& masked) {
    const auto opened = dataset.OpenImage(image_id);
    Require(static_cast<bool>(opened), image_id + " could not be opened");
    const auto& image = opened.value();
    const auto& descriptor = image.descriptor();
    Require(descriptor.pixel_mask_id == flag_id, image_id + " should take " + flag_id +
                                                     ", which data_groups names as its flag, and took '" +
                                                     descriptor.pixel_mask_id + "'");

    for (std::uint64_t polarization = 0; polarization < 3; ++polarization) {
        const auto counted = CountedPerChannel(image, polarization);
        for (std::uint64_t frequency = 0; frequency < 2; ++frequency) {
            std::vector<float> dropped;
            for (const auto& pixel : masked) {
                if (pixel.frequency == frequency && pixel.polarization == polarization) {
                    dropped.push_back(pixel.value);
                }
            }
            const auto where = image_id + " at frequency " + std::to_string(frequency) + ", polarization " +
                               std::to_string(polarization);
            const auto unmasked = ReadPlane(image, frequency, polarization, false);
            const auto pixels = ReadPlane(image, frequency, polarization, true);
            const auto valid = static_cast<std::uint64_t>(
                std::count_if(pixels.begin(), pixels.end(), [](float value) { return std::isfinite(value); }));
            Require(valid == kPlanePixels - dropped.size(), "a masked read of " + where + " kept " +
                                                                std::to_string(valid) + " pixels, expected " +
                                                                std::to_string(kPlanePixels - dropped.size()));
            for (std::size_t i = 0; i < pixels.size(); ++i) {
                const bool flagged = std::find(dropped.begin(), dropped.end(), unmasked.at(i)) != dropped.end();
                Require(std::isfinite(unmasked.at(i)), "an unmasked read of " + where + " lost a pixel");
                Require(flagged ? std::isnan(pixels.at(i)) : pixels.at(i) == unmasked.at(i),
                        "a masked read of " + where + " kept a flagged pixel or dropped a good one");
            }
            Require(counted.at(frequency) == static_cast<double>(kPlanePixels - dropped.size()),
                    "a reduction of " + where + " counted " + std::to_string(counted.at(frequency)) +
                        " pixels, expected " + std::to_string(kPlanePixels - dropped.size()));
        }
    }
}

void TestFlaggedConversion() {
    Require(std::filesystem::exists(kFlaggedFixture),
            "the flagged conformance fixture is missing; run tests/data/generate_conformance_fixtures.py");
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    auto dataset = carta::zarr::Dataset::Open(context.value(), kFlaggedFixture.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed on XRADIO's conversion of masked CASA images");

    std::vector<std::string> images;
    for (const auto& image : dataset.value().descriptor().images) {
        Require(image.openable, image.id + " was discovered but is not openable");
        images.push_back(image.id);
    }
    std::sort(images.begin(), images.end());
    Require((images == std::vector<std::string>{"SKY", "SKY_OTHER"}),
            "discovery should report the two sky images and neither flag");

    TestFlagsMaskTheirOwnImage(dataset.value(), "SKY", "FLAG_SKY", kSkyMasked);
    TestFlagsMaskTheirOwnImage(dataset.value(), "SKY_OTHER", "FLAG_SKY_OTHER", kOtherMasked);
}

// A continuum image of a spectral index, converted from FITS by XRADIO 1.3.0. Its one channel has no
// spacing to give an increment, so the channel width XRADIO writes is the increment; and BTYPE is
// what the image measures.
void TestContinuumConversion() {
    const std::filesystem::path fixture{CARTA_ZARR_CONFORMANCE_CONTINUUM_FIXTURE};
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    auto dataset = carta::zarr::Dataset::Open(context.value(), fixture.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed on XRADIO's continuum image");
    auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "the continuum SKY could not be opened");
    const auto& sky = image.value().descriptor();
    Require(sky.sub_type == "SpectralIndex", "BTYPE 'Spectral index' was read as '" + sky.sub_type + "'");
    Require(sky.spectral.has_value() && sky.spectral->values.size() == 1,
            "the continuum image did not have one channel");
    Require(sky.spectral->increment && *sky.spectral->increment == 1.0e6,
            "the channel width XRADIO wrote was not the increment of the one channel");
    Require(sky.spectral->reference_pixel == 1.0 && sky.spectral->reference_value == 1.4e9,
            "the increment was not measured from the one channel");
}

// The same cube with OBSRA and OBSDEC apart from its reference direction, and cards XRADIO keeps
// uninterpreted: each is reported, by its lower-cased name, with the kind of value it was written as.
void TestKeywordsConversion() {
    const std::filesystem::path fixture{CARTA_ZARR_CONFORMANCE_KEYWORDS_FIXTURE};
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    auto dataset = carta::zarr::Dataset::Open(context.value(), fixture.string());
    Require(static_cast<bool>(dataset), "Dataset::Open failed on XRADIO's image with keywords");
    auto image = dataset.value().OpenImage("SKY");
    Require(static_cast<bool>(image), "the SKY with keywords could not be opened");
    const auto& sky = image.value().descriptor();

    Require(sky.observation && sky.observation->pointing_center, "OBSRA and OBSDEC were not reported");
    const auto& pointing = *sky.observation->pointing_center;
    RequireClose(pointing.value.at(0), 45.5, 1.0e-9, "OBSRA");
    RequireClose(pointing.value.at(1), 30.25, 1.0e-9, "OBSDEC");
    Require(pointing.reference_frame == "FK5", "the pointing centre's frame was not FK5");

    const auto value_of = [&](const std::string& name) -> const carta::zarr::UserKeyword& {
        for (const auto& keyword : sky.user_keywords) {
            if (keyword.name == name) {
                return keyword;
            }
        }
        throw std::runtime_error("user keyword '" + name + "' was not reported");
    };
    Require(std::get<std::string>(value_of("instrume").value) == "Test instrument", "INSTRUME was not a string");
    Require(std::get<double>(value_of("distance").value) == 101.5, "DISTANCE was not a number");
    Require(std::get<std::int64_t>(value_of("ncombine").value) == 3, "NCOMBINE was not an integer");
    Require(std::get<bool>(value_of("pbcor").value), "PBCOR was not a boolean");
    Require(sky.user_keywords.size() == 4, "a card XRADIO interprets was reported as a user keyword");
}

// The amplitude of the cube's transform, an image of the aperture plane as CASA's ia.fft makes one,
// converted by XRADIO 1.3.0. u and v play the spatial roles, so the reads and reductions reach them
// as they reach l and m; it has an aperture coordinate and no direction; its mask is FLAG_APERTURE,
// found by the name XRADIO gives it since no data group can name it; and it has no beams.
void TestApertureConversion() {
    const std::filesystem::path fixture{CARTA_ZARR_CONFORMANCE_APERTURE_FIXTURE};
    const auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    auto dataset = carta::zarr::Dataset::Open(context.value(), fixture.string());
    Require(static_cast<bool>(dataset),
            "Dataset::Open failed on XRADIO's aperture-plane image: " + (dataset ? "" : dataset.error().message));
    Require(dataset.value().descriptor().default_image_id == "APERTURE",
            "the aperture-plane image was not the one the dataset opens on");
    auto opened = dataset.value().OpenImage("APERTURE");
    Require(static_cast<bool>(opened), "APERTURE could not be opened");
    const auto& image = opened.value();
    const auto& aperture = image.descriptor();

    Require(aperture.axes.size() == 5 && aperture.axes.at(0).name == "u" &&
                aperture.axes.at(0).role == AxisRole::spatial_x && aperture.axes.at(0).length == 5 &&
                aperture.axes.at(1).name == "v" && aperture.axes.at(1).role == AxisRole::spatial_y &&
                aperture.axes.at(1).length == 4,
            "u and v were not the spatial axes, first");
    Require(!aperture.direction.has_value(), "an aperture-plane image reported a direction");
    Require(aperture.aperture.has_value(), "an aperture-plane image reported no aperture coordinate");
    const auto& plane = *aperture.aperture;
    Require(plane.unit == "lambda", "the aperture coordinate's unit was '" + plane.unit + "'");
    // u runs from +2 cells to -2 through zero at its third sample, v from -2 to +1 through zero at its
    // third; XRADIO writes zero as each one's crval.
    RequireClose(plane.reference_pixel.at(0), 3.0, 1.0e-9, "u reference pixel");
    RequireClose(plane.reference_pixel.at(1), 3.0, 1.0e-9, "v reference pixel");
    RequireClose(plane.reference_value.at(0), 0.0, 1.0e-9, "u reference value");
    RequireClose(plane.reference_value.at(1), 0.0, 1.0e-9, "v reference value");
    RequireClose(plane.increment.at(0), -11459.155902616463, 1.0e-6, "u increment");
    RequireClose(plane.increment.at(1), 14323.94487827058, 1.0e-6, "v increment");
    Require(aperture.spectral.has_value() && aperture.polarization.has_value(),
            "the aperture-plane image lost its spectral or polarization coordinate");

    Require(aperture.pixel_mask_id == "FLAG_APERTURE", "FLAG_APERTURE was not found as the image's mask");
    // Two pixels are masked: (u 0, v 0) at the first channel and polarization, (u 2, v 1) at the second
    // of each. A plane is read with u fastest.
    const auto unmasked = ReadPlane(image, 0, 0, false);
    const auto masked = ReadPlane(image, 0, 0, true);
    Require(std::isfinite(unmasked.at(0)) && std::isnan(masked.at(0)), "the masked pixel was not masked");
    RequireClose(masked.at(2 + (5 * 2)), 190.0, 1.0e-3, "the amplitude at the centre of the plane");
    const auto first = CountedPerChannel(image, 0);
    const auto second = CountedPerChannel(image, 1);
    Require(first.at(0) == 19.0 && first.at(1) == 20.0 && second.at(0) == 20.0 && second.at(1) == 19.0,
            "a reduction over the aperture plane did not count its masked pixels out");

    const auto beams = image.ReadBeams();
    Require(beams && beams.value().empty(), "an aperture-plane image reported beams");
}

}  // namespace

int main() {
    try {
        const auto image = OpenSky();
        const auto& sky = image.descriptor();
        TestIdentityAndAxes(sky);
        TestDirection(sky);
        TestPointingCentreWithoutObsra(sky);
        TestSpectralAndPolarization(sky);
        TestTemporalAndStorage(sky, image.chunk_geometry());

        const auto beams = image.ReadBeams();
        Require(static_cast<bool>(beams), "BMAJ/BMIN/BPA did not survive as a readable beam");
        Require(!beams.value().empty(), "the beam table decoded to no beams");

        TestFlaggedConversion();
        TestContinuumConversion();
        TestKeywordsConversion();
        TestApertureConversion();

        std::cout << "carta-zarr XRADIO conformance tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr XRADIO conformance tests failed: " << error.what() << '\n';
        return 1;
    }
}
