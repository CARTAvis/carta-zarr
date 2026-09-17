/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_

#include "../../zarr/array_metadata.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::xradio {

/**
 * The beams a beam fit parameter table holds, one per (time, channel, polarization).
 *
 * Takes the arrays already read rather than a Store, which is the one thing about this module worth
 * stating: reading them is two value reads, and a caller that has to make them cannot be tested
 * without an array on disk. Handed the values instead, everything this decides -- which dimension is
 * which, that the labels name a major axis, a minor axis and a position angle, that a table missing
 * a dimension is an error rather than an empty list, and the order the three loops emit in -- is
 * checkable from two vectors. See ADR 0006.
 *
 * `values` is the table's flat C-order buffer, addressed through the metadata by dimension name.
 * A table with no time dimension reads as a single plane. The unit comes from the table's own
 * `units` attribute, so a caller passes nothing it would have to read twice.
 *
 * `beam_node` is passed only because ArrayMetadata does not carry the name of the node it came
 * from, and an error about the table should say which variable it is about.
 */
Result<std::vector<Beam>> DescribeBeams(const zarr::ArrayMetadata& beam_metadata, std::string_view beam_node,
                                        const std::vector<std::string>& parameter_labels,
                                        const std::vector<double>& values);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_BEAM_TABLE_H_
