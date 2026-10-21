/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_DESCRIPTOR_H_
#define CARTA_ZARR_DESCRIPTOR_H_

// What the library reports about an image dataset and the images in it: the result of probing a
// location, and the immutable metadata of what it found. Never pixels.

#include "carta-zarr/error.h"

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace carta::zarr {

// What a diagnostic is about. One closed set, so that a code is a name the compiler checks rather
// than a string two places could spell differently -- which is what it was, made from ErrorCode's
// names in some places and written out as literals in others.
//
// Several share a name with an ErrorCode, and mean the same thing: the diagnostic is what the error
// would have been, had the rest of the store not still been readable.
enum class DiagnosticCode {
    invalid_metadata,
    unsupported_data_type,
    ambiguous_schema,
    // A node the store holds that is neither a group nor an array this library recognizes.
    unrecognised_node,
    // An array whose metadata would not parse.
    unreadable_array,
    // A coordinate whose samples are not evenly spaced, so it has no linear description.
    nonuniform_axis,
    // A coordinate with no sample at its reference world value, so its reference pixel was
    // extrapolated rather than read off.
    inexact_reference_pixel,
    // A coordinate with fewer than two distinct samples, which has no linear description at all.
    degenerate_axis,
    // A store naming no schema version, or one this library does not read.
    unsupported_schema_version,
};

// The code's name, as it would be spelled in a message or a log.
inline const char* DiagnosticCodeName(DiagnosticCode code) noexcept {
    switch (code) {
        case DiagnosticCode::invalid_metadata:
            return "invalid_metadata";
        case DiagnosticCode::unsupported_data_type:
            return "unsupported_data_type";
        case DiagnosticCode::ambiguous_schema:
            return "ambiguous_schema";
        case DiagnosticCode::unrecognised_node:
            return "unrecognised_node";
        case DiagnosticCode::unreadable_array:
            return "unreadable_array";
        case DiagnosticCode::nonuniform_axis:
            return "nonuniform_axis";
        case DiagnosticCode::inexact_reference_pixel:
            return "inexact_reference_pixel";
        case DiagnosticCode::degenerate_axis:
            return "degenerate_axis";
        case DiagnosticCode::unsupported_schema_version:
            return "unsupported_schema_version";
    }
    return "unknown";
}

struct Diagnostic {
    DiagnosticCode code = DiagnosticCode::invalid_metadata;
    std::string message;
    std::string node_path;
};

using SchemaId = std::string;

inline constexpr std::string_view kXradioImageSchema = "xradio.image";

enum class SchemaMatchKind {
    no_match,
    match,
    invalid,
};

struct SchemaProbeResult {
    SchemaMatchKind kind = SchemaMatchKind::no_match;
    SchemaId schema_id;
    std::string schema_version;
    std::vector<Diagnostic> diagnostics;
};

// What Context::Create takes: the resources every dataset and image opened through that context
// share. Named for the context rather than for opening, which it was, because it opens nothing --
// Dataset::Open takes a location and the context, and none of these.
struct ContextOptions {
    // How much decoded-chunk cache this context may hold, in bytes.
    //
    // Three answers, and one field because they are one question: no value leaves TensorStore's own
    // default alone, zero asks for a pool that holds nothing, and any other number sizes it. It was
    // two fields -- a size where zero meant "default" and a separate disable_cache -- which made
    // "no cache" sayable twice and, when both were set, resolved silently in favor of the size.
    std::optional<std::size_t> cache_bytes;
    unsigned int io_threads = 0;
    unsigned int decode_threads = 0;
};

enum class AxisRole {
    spatial_x,
    spatial_y,
    spectral,
    polarization,
    time,
    other,
};

struct AxisDescriptor {
    std::string name;
    AxisRole role = AxisRole::other;
    std::uint64_t length = 0;
    std::size_t storage_index = 0;
};

// Where the axis playing `role` sits among `axes`, or none when no axis plays it. Asked of
// ImageDescriptor::axes or ImageEntry::axes, and the answer is also the index of that axis's Range
// in a ReadRequest.
//
// Axes are reached by role because their logical order is a schema profile's choice rather than
// this library's promise. The XRADIO profile happens to report l, m, frequency, polarization, time,
// whatever order the store holds them in; a consumer that indexes by position is assuming that
// profile, which is why no such position is exported. See ADR 0012.
//
// Every role but `other` is played by at most one axis of an image this library describes. For
// `other` this finds the first.
inline std::optional<std::size_t> AxisIndex(const std::vector<AxisDescriptor>& axes, AxisRole role) noexcept {
    for (std::size_t index = 0; index < axes.size(); ++index) {
        if (axes[index].role == role) {
            return index;
        }
    }
    return std::nullopt;
}

// One variable of a dataset's image listing, and what can be said about it without opening it.
struct ImageEntry {
    std::string id;
    bool openable = false;
    std::vector<Diagnostic> diagnostics;
    // The axes Dataset::OpenImage would report for this image -- ImageDescriptor::axes, element for
    // element -- and empty when it is not openable.
    //
    // Here so that a consumer can decide whether it can use an image from the listing alone. One that
    // displays a single time step, say, refuses an image with two, and without these it had to open
    // every image in the dataset to find out, reading every coordinate value to answer a question
    // about shapes. They come from metadata the listing has already parsed.
    std::vector<AxisDescriptor> axes;
};

// What this library knows about an image dataset without opening any image in it: which schema
// profile describes it, which variables in it are images, and whatever that profile had to say about
// a store it nonetheless accepted.
struct DatasetDescriptor {
    SchemaId schema_id;
    std::string schema_version;
    std::vector<ImageEntry> images;
    std::optional<std::string> default_image_id;
    std::vector<Diagnostic> diagnostics;
};

enum class DataType {
    unknown,
    boolean,
    int8,
    uint8,
    int16,
    uint16,
    int32,
    uint32,
    int64,
    uint64,
    float16,
    float32,
    float64,
    complex64,
    complex128,
};

struct DirectionCoordinate {
    std::string projection;
    std::string reference_frame;  // e.g. "FK5", "ICRS", "GALACTIC"
    std::optional<double> equinox;
    // CRPIX, and 1-based as FITS counts it: the first pixel of an axis is 1.0, not 0.0.
    //
    // Stated rather than implied because a consumer that reads it the other way shifts every
    // coordinate it builds by a whole pixel and is told nothing -- which is the failure ADR 0002
    // describes for the projection parameters, arriving by a different route. casacore counts from
    // zero, so a consumer handing this to a DirectionCoordinate subtracts one.
    //
    // SpectralCoordinate::reference_pixel is the same convention, from the same fit.
    std::array<double, 2> reference_pixel{0.0, 0.0};
    std::array<double, 2> reference_value{0.0, 0.0};                                       // CRVAL in degrees
    std::array<double, 2> increment{0.0, 0.0};                                             // CDELT in degrees
    std::array<std::array<double, 2>, 2> transformation_matrix{{{1.0, 0.0}, {0.0, 1.0}}};  // PC matrix
    std::vector<double> projection_parameters;
    std::array<double, 2> native_pole_direction{0.0, 0.0};  // longPole/latPole in degrees
};

// The linear coordinate of an image of the aperture plane, the Fourier conjugate of the sky's: its
// spatial_x and spatial_y axes are u and v rather than l and m, in wavelengths. Both axes share the
// unit, and each array holds u's value first.
//
// Linear by construction, as a direction axis is, so the increment is reported even when the samples
// are not evenly spaced; a diagnostic says so.
struct ApertureCoordinate {
    std::string unit;                                 // "lambda", as XRADIO writes it
    std::array<double, 2> reference_pixel{0.0, 0.0};  // 1-based, as DirectionCoordinate::reference_pixel
    std::array<double, 2> reference_value{0.0, 0.0};
    std::array<double, 2> increment{0.0, 0.0};
};

struct SpectralCoordinate {
    std::string unit;
    std::string
        reference_frame;  // casacore or FITS SPECSYS name, e.g. "LSRK", "BARY", "GEO", "TOPOCENT"; never astropy's
    // CRPIX3, 1-based as DirectionCoordinate::reference_pixel is. Absent, with the two below it,
    // when the channels are not evenly spaced: there is then no linear description to give and a
    // consumer builds a tabular axis from values instead.
    std::optional<double> reference_pixel;
    std::optional<double> reference_value;  // CRVAL3
    std::optional<double> increment;        // CDELT3
    std::optional<double> rest_frequency;
    // One a channel, in `unit`, as TemporalCoordinate::values are in theirs. Named for what is
    // checked of them: the unit is the frequency coordinate's own, which the schema holds only to
    // being given, so nothing here promises Hz.
    std::vector<double> values;
};

struct TemporalCoordinate {
    // As the dataset wrote them: `unit`, `scale` and `format` say how to read them. XRADIO writes
    // MJD days as often as unix seconds, and nothing here converts one to the other.
    std::vector<double> values;
    std::string unit;
    std::string scale;
    std::string format;
};

struct PolarizationCoordinate {
    std::vector<std::string> labels;
};

// A direction on the sky given outside the coordinate system, as a measure: its two angles and the
// frame they are in, each as DirectionCoordinate gives its own.
struct SkyDirection {
    std::string reference_frame;  // as DirectionCoordinate::reference_frame, e.g. "FK5", "ICRS"
    std::optional<double> equinox;
    std::array<double, 2> value{0.0, 0.0};  // longitude and latitude in degrees
};

struct ObservationInfo {
    std::string object_name;
    std::string observer;
    std::string telescope_name;
    std::string time_scale;  // FITS TIMESYS, e.g. "UTC", "TAI"
    // The observation date as an MJD in `time_scale`, from a number in the MJD or unix format. Absent
    // when the dataset wrote none, or wrote it in any other form.
    std::optional<double> observation_mjd;
    std::optional<std::array<double, 3>> observatory_position;  // OBSGEO-X, Y, Z (meters)
    // Where the telescope pointed: FITS OBSRA and OBSDEC. Reported whenever the store has one, even
    // when it equals the reference direction, which is what XRADIO's FITS reader writes for a file
    // with no OBSRA: the store does not say which it was.
    std::optional<SkyDirection> pointing_center;
};

// A keyword the dataset carries that this library does not interpret: a FITS card XRADIO's reader
// left alone, lower-cased, or an entry of a CASA image's miscinfo.
//
// An integer is one JSON wrote as one and a signed 64-bit integer holds. JSON's parser reads an
// integer beyond 64 bits either way as a number, so such a one arrives as the nearest double.
struct UserKeyword {
    std::string name;
    std::variant<bool, std::int64_t, double, std::string> value;
};

/// The read geometry of one image, reported in the logical axis order of ImageDescriptor::axes so
/// that a consumer never has to undo the stored order itself.
///
/// Two granularities, deliberately separate. An inner chunk is what must be decoded to reach any byte
/// inside it, so a consumer sizing a cache reasons about chunk_shape. A shard is one file, holding an
/// index and the inner chunks, so a consumer counting files -- and the opens and stats a shared
/// metadata server answers for them -- reasons about shard_shape. They are equal when the array is not
/// sharded, and can differ by a large factor when it is.
///
/// A shard is not one request. From a local store, the inner chunks one storage read needs from a
/// shard are read as byte ranges of its file, after its index, which is cached: ranges at most a few
/// hundred bytes apart are read as one, and a storage read that needs every inner chunk of a shard
/// reads the whole file at once. Image::Read and the reductions cut what they read into several
/// storage reads (see ReadOptions::read_budget_bytes), and each is its own.
struct ChunkGeometry {
    // The spatial axis the store varies fastest, which is the one a reduction walks along.
    //
    // Reading a plane with the other one fastest means transposing every chunk on the way into the
    // destination, and that is not a rounding error: on the same 2048x2048x16 image stored both ways,
    // a whole-plane spectral profile took a third to a half longer along the wrong axis.
    //
    // So the walk follows the store rather than the other way round, and the runs a reduction makes
    // from a region's raster lie along this axis.
    AxisRole fastest_spatial_axis = AxisRole::spatial_x;
    std::vector<std::uint64_t> chunk_shape;
    std::vector<std::uint64_t> shard_shape;
    bool sharded = false;
    std::string compressor;
};

struct Beam {
    // The plane this beam was fitted on. Every plane is reported; a consumer that handles one time
    // step selects it rather than being handed it.
    std::size_t time = 0;
    std::size_t channel = 0;
    std::size_t polarization = 0;
    double major = 0.0;
    double minor = 0.0;
    double position_angle = 0.0;
    std::string unit;
};

// Which of two questions a reported size answers.
//
// Said rather than resolved into "how close is this", because the library cannot tell: the declared
// answer is what a caller gets when the store could not be measured, so there is no measurement to
// compare it against. Nor can it be inferred -- compression pushes what a store occupies below what
// its arrays declare, while per-node metadata and shard indices push it above, and which wins is a
// property of the store this path could not read. See ADR 0008.
//
// `declared` is first so that a DatasetSize nobody filled in claims the weaker of the two.
enum class SizeBasis {
    // Read from the metadata: every array's shape times its element size, uncompressed. What the
    // dataset says it holds, which is not a bound on what the store occupies.
    declared,
    // Walked: the sum of the sizes of the files under the store.
    measured,
};

struct DatasetSize {
    std::uint64_t bytes = 0;
    SizeBasis basis = SizeBasis::declared;
};

struct ImageDescriptor {
    std::string id;
    // What a sky image measures, as XRADIO image schema 0.0.2 names it from casacore's image type or
    // FITS BTYPE ("Intensity", "SpectralIndex", "OpticalDepth", ...). Empty when the image names none.
    std::string sub_type;
    std::vector<std::string> data_groups;
    DataType stored_type = DataType::unknown;
    // In the image's logical order, which a request's ranges follow. Find one with AxisIndex rather
    // than by position; storage_index says where each one lives on disk.
    std::vector<AxisDescriptor> axes;
    std::string unit;
    // The flag variable supplying this image's pixel mask, empty when it has none. Reported for the
    // same reason `id` is: it names a data variable the consumer may want to see in diagnostics.
    std::string pixel_mask_id;
    // Which plane the image is on: an image of the sky has a direction, and one of the aperture plane
    // an aperture coordinate. Never both. The spatial_x and spatial_y axes are that plane's (ADR 0024).
    std::optional<DirectionCoordinate> direction;
    std::optional<ApertureCoordinate> aperture;
    std::optional<SpectralCoordinate> spectral;
    std::optional<PolarizationCoordinate> polarization;
    std::optional<TemporalCoordinate> temporal;
    std::optional<ObservationInfo> observation;
    // The dataset's `user` keywords and the image's own, the image's winning where both name one, by
    // name. A value that is not a boolean, an integer, a number or a string -- a dictionary or a list,
    // as a CASA image's miscinfo can hold -- is left out, and a diagnostic names it.
    std::vector<UserKeyword> user_keywords;
    std::vector<Diagnostic> diagnostics;
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_DESCRIPTOR_H_
