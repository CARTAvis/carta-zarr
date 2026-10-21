/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_

// The coordinates of an XRADIO image dataset, in one place: which there are, what each is called and
// holds, what makes one well formed, whether an image agrees with them, what unit each is in, and
// how their values are read. One place, so that the listing and a descriptor cannot disagree about
// an image's axes, nor two places about a coordinate's unit.
//
// Who asks what follows ADR 0009: the probe asks whether each coordinate is well formed, a fact about
// the dataset; qualification asks whether an image agrees with them, a fact about the image; and
// describing an image reads them.
//
// It takes a Store, as qualification does, because that is the narrowest input its tests can build
// for everything but the values (ADR 0006).

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
    numeric,  // real numbers, read as doubles: l, m, u, v, frequency, time
    labels,   // fixed-length UTF-32 strings, decoded by this library: polarization
};

// Which plane's variables carry a coordinate. An image is on the plane its `type` belongs to, and
// carries every coordinate of it: the sky's, or the aperture plane's, its Fourier conjugate
// (ADR 0024).
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

// Every coordinate this profile knows, in the order it reports an image's axes: the plane's spatial
// pair, then the three every plane shares, whatever order the store holds them in. The order is the
// profile's choice; the library only promises that an axis can be found by its role (ADR 0012). u and
// v play the spatial roles on the aperture plane as l and m do on the sky (ADR 0024).
//
// A dataset with an image of a plane must have every coordinate of that plane. Both XRADIO readers
// write them all, so a dataset missing one is malformed.
inline constexpr std::array<Coordinate, 7> kCoordinates{{
    {"l", AxisRole::spatial_x, CoordinateKind::numeric, Plane::sky},
    {"m", AxisRole::spatial_y, CoordinateKind::numeric, Plane::sky},
    {"u", AxisRole::spatial_x, CoordinateKind::numeric, Plane::aperture},
    {"v", AxisRole::spatial_y, CoordinateKind::numeric, Plane::aperture},
    {"frequency", AxisRole::spectral, CoordinateKind::numeric, Plane::both},
    {"polarization", AxisRole::polarization, CoordinateKind::labels, Plane::both},
    {"time", AxisRole::time, CoordinateKind::numeric, Plane::both},
}};

constexpr bool OnPlane(const Coordinate& coordinate, Plane plane) noexcept {
    return coordinate.plane == Plane::both || coordinate.plane == plane;
}

// The coordinate of `plane` playing `role`. Every role but `other` has exactly one on each plane.
constexpr const Coordinate& PlaneCoordinate(AxisRole role, Plane plane) noexcept {
    for (const auto& coordinate : kCoordinates) {
        if (coordinate.role == role && OnPlane(coordinate, plane)) {
            return coordinate;
        }
    }
    return kCoordinates.front();
}

// The sky-plane coordinate playing `role`.
constexpr const Coordinate& SkyCoordinate(AxisRole role) noexcept {
    return PlaneCoordinate(role, Plane::sky);
}

// How many coordinates a variable on `plane` carries. They are all it is described by, so a variable
// with more dimensions than this is not one this profile can read.
std::size_t AxisCount(Plane plane) noexcept;

// Whether a variable names every coordinate of `plane` among its dimensions.
bool CarriesPlane(const zarr::ArrayMetadata& variable, Plane plane);

// The plane an image of this `type` is on, or nothing when schema 0.0.2 has no image of that type.
std::optional<Plane> ImagePlane(std::string_view type);

// The plane an image is on: the one its `type` belongs to, when it carries every coordinate of that
// plane, and nothing when its type is no image's or it does not.
std::optional<Plane> PlaneOf(const zarr::ArrayMetadata& variable);

// Whether a coordinate array is well formed in itself: one-dimensional, over an axis of its own name,
// and holding what its kind says (a real type, or labels in a layout this library can decode). Asked
// of the metadata alone; whether it exists, and whether it could be read, is for the caller to say.
//
// A label array this library cannot decode fails every image that opens, which is why it is a
// question about the coordinate and not left to the read.
Result<void> CheckCoordinate(const zarr::ArrayMetadata& array, const Coordinate& coordinate);

// Why an image cannot be described with the coordinates it names, or nothing when it can. A dataset
// stores each coordinate once and every image references it by dimension name, so an image whose
// extent disagrees with the coordinate would be reported with the wrong coordinate values, such as
// three channels' worth over seven channels of pixels.
//
// A coordinate the dataset does not carry, or whose metadata will not read, is passed over: every
// image references it, so the probe closes the dataset over it rather than one image.
std::optional<Diagnostic> ExtentDisagreement(const Store& store, const zarr::ArrayMetadata& image,
                                             std::string_view node);

// The unit a coordinate's samples are in (its units attribute), or empty when the dataset does not
// say. A descriptor reports it only in that coordinate.
std::string CoordinateUnit(const Store& store, const Coordinate& coordinate);

// The axes of an image, in logical order, each with its length and where it sits in the stored order.
// The same list for the listing and for a descriptor, from the image's metadata alone.
std::vector<AxisDescriptor> DescribeAxes(const zarr::ArrayMetadata& image);

// The samples of the dataset's coordinates, which is everything a description reads besides metadata.
// ReadCoordinateValues fills every one of the plane's; a test may leave one empty, and then no
// coordinate of that kind is described.
struct CoordinateValues {
    std::vector<double> l;
    std::vector<double> m;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> frequency;
    std::vector<double> time;
    std::vector<std::string> polarization;
};

// The samples of every coordinate of `plane`, each read through the store and so checked against its
// own document (Store::VerifyArray). Those of the other plane are left empty: a dataset need not have
// them.
Result<CoordinateValues> ReadCoordinateValues(const Store& store, Plane plane);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_COORDINATES_H_
