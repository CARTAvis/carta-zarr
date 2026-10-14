/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_ROW_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_ROW_H_

// Where the fields that differ between versions of XRADIO's image schema live: one row per version
// this library reads, and the rule that picks a store's row (ADR 0022).
//
// A row holds only what differs. A rule every version shares is code, written once where it is
// applied; a row says where to look, and is code only where a version needs a translation rather
// than a place.

#include "carta-zarr/descriptor.h"

#include "../../store.h"

#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

// Where a spectral coordinate's frame is read from.
enum class SpectralFrameSource {
    // frequency.attrs.frame, in casacore's names.
    frame_attribute,
    // frequency.attrs.reference_frequency.attrs.observer, in astropy's names where astropy has the
    // frame and casacore's otherwise. Upper-cased, with GCRS read as casacore's GEO.
    reference_observer,
};

struct SchemaRow {
    // The schema_version a store names to be read by this row; empty for the row of stores that name
    // none.
    std::string_view version;
    SpectralFrameSource spectral_frame;
    // Whether a frequency coordinate with no units of its own takes its reference frequency's.
    bool unit_from_reference_frequency;
    // Whether an observation date written as a number with no format is an MJD, as XRADIO 1.2 has
    // always been read. Schema 0.0.2 requires the format.
    bool unformatted_date_is_mjd;
    // Whether a frequency coordinate's channel_width is the increment of an axis of one channel, which
    // its one channel cannot give.
    bool increment_from_channel_width;
    // Whether a sky image's sub_type, what it measures, is read.
    bool reads_sub_type;
    // XRADIO's names for a dataset's images, in the order they are listed. Any other is listed after
    // them, by name.
    std::vector<std::string_view> known_images;
};

// The rows, newest first. Inline, so that what reads a row's fields needs no store to link.
inline const std::vector<SchemaRow>& SchemaRows() {
    static const std::vector<SchemaRow> rows{
        SchemaRow{"0.0.2",
                  SpectralFrameSource::frame_attribute,
                  false,
                  false,
                  true,
                  true,
                  {"SKY", "SKY_MODEL", "SKY_RESIDUAL", "SKY_DIRTY", "POINT_SPREAD_FUNCTION", "PRIMARY_BEAM", "MASK"}},
        // XRADIO 1.2 named a deconvolution's products both ways: MODEL and RESIDUAL from a list of
        // tclean's images, SKY_MODEL and SKY_RESIDUAL from a dict of them.
        SchemaRow{"",
                  SpectralFrameSource::reference_observer,
                  true,
                  true,
                  false,
                  false,
                  {"SKY", "SKY_MODEL", "SKY_RESIDUAL", "SKY_DIRTY", "MODEL", "RESIDUAL", "DIRTY",
                   "POINT_SPREAD_FUNCTION", "PRIMARY_BEAM", "MASK", "MASK_DECONVOLVE"}},
    };
    return rows;
}

// The row a store is read by, and what choosing it found.
struct ChosenRow {
    // Null when the store names a version no row reads.
    const SchemaRow* row = nullptr;
    // As the store names it, which is what a probe reports; empty when it names none.
    std::string version;
    // Why the row is the unversioned one, or why there is none.
    std::vector<Diagnostic> diagnostics;
};

// A store naming 0.0.2 is read by the 0.0.2 row when it is one: its data groups name a sky image, that
// image is typed "sky", and its frequency coordinate, if it has one, gives a frame and units.
// Otherwise, and for a store naming no version, the unversioned row reads it and a diagnostic says
// why. Any other version is no row's.
ChosenRow ChooseSchemaRow(const Store& store);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_SCHEMA_ROW_H_
