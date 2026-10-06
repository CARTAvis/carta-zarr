/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_

// The coordinates of an XRADIO image dataset, said once: which there are, what each is called and
// holds, what makes one well formed, whether an image agrees with them, what unit each is in, and
// reading their values.
//
// Said once so that an axis and its coordinate cannot disagree, about the unit above all: one rule
// for it, so that no image reports its spectral axis unitless and its spectral coordinate in Hz.
//
// Who asks what is ADR 0009's division: the probe asks whether each coordinate is
// well formed, which is a fact about the dataset; qualification asks whether an image agrees with
// them, which is a fact about the image; describing an image reads them.
//
// It takes a Store, as qualification does, because that is the narrowest thing its tests can stand
// up for everything but the values (ADR 0006).

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include "../../store.h"
#include "../../zarr/array_metadata.h"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

// What a coordinate's samples are.
enum class CoordinateKind {
    numeric,  // real numbers, read as doubles: l, m, frequency, time
    labels,   // fixed-length UTF-32 strings, decoded by this library: polarization
};

// Which plane's variables carry a coordinate. A variable carrying every coordinate of a plane is on
// that plane: the sky plane's are images, the aperture plane's are listed and not opened.
enum class Plane {
    both,
    sky,
    aperture,
};

struct Coordinate {
    std::string_view name;
    AxisRole role;
    CoordinateKind kind;
    Plane plane;
};

// Every coordinate this profile knows. Those of the sky plane come first, in the order this profile
// reports an image's axes in, whatever order the store holds them in -- the profile's choice rather
// than the library's promise, which is only that an axis can be found by its role (ADR 0012).
//
// Every sky-plane coordinate is required. XRADIO writes optional coordinate arrays over the spatial
// pair alone -- right_ascension and declination are float64 over (l, m) and carry no type attribute
// -- so a rule keyed on "has l and m" offers those to a consumer as openable images (ADR 0001). Both
// XRADIO readers write all five unconditionally, so a dataset missing one is malformed.
inline constexpr std::array<Coordinate, 7> kCoordinates{{
    {"l", AxisRole::spatial_x, CoordinateKind::numeric, Plane::sky},
    {"m", AxisRole::spatial_y, CoordinateKind::numeric, Plane::sky},
    {"frequency", AxisRole::spectral, CoordinateKind::numeric, Plane::both},
    {"polarization", AxisRole::polarization, CoordinateKind::labels, Plane::both},
    {"time", AxisRole::time, CoordinateKind::numeric, Plane::both},
    {"u", AxisRole::other, CoordinateKind::numeric, Plane::aperture},
    {"v", AxisRole::other, CoordinateKind::numeric, Plane::aperture},
}};

constexpr bool OnPlane(const Coordinate& coordinate, Plane plane) noexcept {
    return coordinate.plane == Plane::both || coordinate.plane == plane;
}

// The sky-plane coordinate playing `role`. Every role but `other` has exactly one.
constexpr const Coordinate& SkyCoordinate(AxisRole role) noexcept {
    for (const auto& coordinate : kCoordinates) {
        if (coordinate.role == role && OnPlane(coordinate, Plane::sky)) {
            return coordinate;
        }
    }
    return kCoordinates.front();
}

// How many coordinates a variable on `plane` carries: the whole of what it is described by, so a
// variable with more dimensions than this is not one this profile can read.
std::size_t AxisCount(Plane plane) noexcept;

// Whether a variable names every coordinate of `plane` among its dimensions.
bool CarriesPlane(const zarr::ArrayMetadata& variable, Plane plane);

// Whether a coordinate array is well formed in itself: one-dimensional, over an axis of its own name,
// holding what its kind says -- a real type, or labels in a layout this library can decode. Asked of
// the metadata alone; whether it exists, and whether it could be read, is the asker's to say.
//
// A label array this library cannot decode fails every image that opens, which is why it is a
// question about the coordinate and not left to the read.
Result<void> CheckCoordinate(const zarr::ArrayMetadata& array, const Coordinate& coordinate);

// Why an image cannot be described with the coordinates it names, or nothing when it can. A dataset
// stores each coordinate once and every image references it by dimension name, so an image whose own
// extent disagrees with the coordinate it names would be reported with another image's coordinate
// vector, three channels' worth over seven channels of pixels.
//
// A coordinate the dataset does not carry, or whose metadata will not read, is passed over: every
// image references it, so the probe closes the dataset over it rather than one image.
std::optional<Diagnostic> ExtentDisagreement(const Store& store, const zarr::ArrayMetadata& image,
                                             std::string_view node);

// The unit a coordinate's samples are in, or empty when the dataset does not say. Its units
// attribute; for the spectral coordinate, failing that, its reference frequency's, which is where
// XRADIO writes it. Both the axis and the spectral coordinate report this, so they cannot disagree.
std::string CoordinateUnit(const Store& store, const Coordinate& coordinate);

// The axes of an image, in logical order, each with its length, its unit and where it sits in the
// stored order. The same list for the listing and for a descriptor, from metadata alone.
std::vector<AxisDescriptor> DescribeAxes(const Store& store, const zarr::ArrayMetadata& image);

// The samples of the dataset's coordinates, which is everything a description reads besides
// metadata. Read, every one is there; one left empty is one handed in that way, and describes no
// coordinate of that kind.
struct CoordinateValues {
    std::vector<double> l;
    std::vector<double> m;
    std::vector<double> frequency;
    std::vector<double> time;
    std::vector<std::string> polarization;
};

// Every sky-plane coordinate's samples, each read through the store and so held to its own document
// (Store::VerifyArray).
Result<CoordinateValues> ReadCoordinateValues(const Store& store);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_
