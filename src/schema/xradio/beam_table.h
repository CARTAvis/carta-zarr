/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include "../../zarr/numeric_array.h"

#include <string>
#include <vector>

namespace carta::zarr::internal::xradio {

/**
 * The beams a beam fit parameter table holds, one per (time, channel, polarization).
 *
 * Takes the arrays already read rather than a Store. Reading them takes two value reads, which
 * cannot be tested without an array on disk. Given the values, everything this decides can be tested
 * from a table and its labels: which dimension is which, that the labels name a major axis, a minor
 * axis and a position angle, that a table missing a dimension is an error rather than an empty list,
 * and the order of the output. See ADR 0006.
 *
 * `table` is addressed by dimension name, through the document its values were decoded with. A
 * table with no time dimension reads as a single plane. The unit comes from the table's own `units`
 * attribute, so a caller passes nothing it would have to read twice.
 */
Result<std::vector<Beam>> DescribeBeams(const zarr::NumericArray& table,
                                        const std::vector<std::string>& parameter_labels);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_
