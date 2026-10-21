/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_PLANE_HISTOGRAM_H_
#define CARTA_ZARR_SRC_REDUCE_PLANE_HISTOGRAM_H_

#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include "reducible_image.h"

namespace carta::zarr::internal {

/**
 * Bin every pixel of each plane over a fixed range, one pass over the pixels.
 *
 * Much simpler than ReduceSpectral because the region is always the whole plane: there is nothing
 * to index, nothing to clip, and no mask to consult, so the walk is a band of chunk rows at a time
 * and the binning runs over the read buffer end to end.
 *
 * Like the reduction, it reads in the store's own axis order: binning does not care what order it
 * sees the pixels in, so there is no reason to transpose them.
 *
 * `workers` runs the binning. Each plane's counts are its own row of the block, so the work is split
 * by rows of the plane, with one private histogram per worker, summed into the row at the end.
 * Splitting by plane instead would be limited to the number of planes one read holds, which for a
 * large image is often one.
 */
Result<void> ComputeHistogram(const ReducibleImage& image, const HistogramRequest& request, const HistogramSink& sink,
                              const ReadOptions& options);

/**
 * One histogram for the whole selection in a single pass. See CubeHistogramRequest.
 *
 * `workers` splits the same way, by rows of a plane, but each worker keeps one provisional histogram
 * for the whole walk rather than one per plane: a growing histogram's range is learned from the
 * pixels it has seen, and starting over on each plane would make it learn the range again.
 *
 * The workers' histograms are never merged with each other. The extremes are tracked exactly, so by
 * the end the target bin edges are known, and each worker re-aggregates onto them independently.
 * That makes the serial pass's error, splitting a provisional bin that straddles a target edge, once
 * per worker instead of once. This is why it is the one reduction whose counts are not exactly the
 * serial answer.
 */
Result<CubeHistogramResult> ComputeCubeHistogram(const ReducibleImage& image, const CubeHistogramRequest& request,
                                                 const ReadOptions& options,
                                                 const CubeHistogramProgressCallback& progress);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PLANE_HISTOGRAM_H_
