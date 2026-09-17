/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_TUNING_H_
#define CARTA_ZARR_SRC_REDUCE_TUNING_H_

// Numbers the reductions were tuned to, which a caller neither sets nor needs.
//
// They were public, in the header that describes what an image dataset is, which published them as
// though they were part of the contract. They are not: nothing outside src/reduce/ has ever read
// one. The limits a caller does have to respect -- kMaxHistogramBins and kMaxSpectralRegions --
// stay in carta-zarr/reduce.h, because a request that exceeds either is rejected.

#include "carta-zarr/reduce.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

// Every statistic, in the order a SpectralBlock lays them out.
inline constexpr std::array<Statistic, 6> kStatisticOrder{Statistic::num_pixels, Statistic::nan_count,
                                                          Statistic::sum, Statistic::sum_sq,
                                                          Statistic::min, Statistic::max};

// The provisional resolution a cube histogram bins at when the caller does not choose one.
//
// Sixteen provisional bins for every bin asked for, because what decides the error is how finely
// the walk resolves one target bin, not how many bins it holds in total -- and because the
// provisional histogram is eight bytes a bin and there is one per worker, so the ones nobody needs
// are paid for in cache. A thousand target bins get 16,384 of them, which is 128 kB.
//
// Held between 4,096 and 65,536. Measured on a billion-pixel ASKAP cube against the exact two-pass
// answer, with the walk on twenty-eight threads: 65,536 misplaced 0.003% of pixels in 3.64 s,
// 16,384 misplaced 0.004% in 2.40 s, 8,192 misplaced 0.006% in 2.21 s, 4,096 misplaced 0.011% in
// 2.14 s, and below that both numbers get worse at once. Sixteen is where the time has flattened
// out, with a few times the resolution the error would need.
inline constexpr std::uint32_t kProvisionalBinsPerBin = 16;
inline constexpr std::uint32_t kLeastProvisionalBins = 1u << 12;
inline constexpr std::uint32_t kMostProvisionalBins = 1u << 16;

// The memory one emitted block may occupy. emit_every_channels is reduced to fit it; see
// SpectralReduceRequest::emit_every_channels.
inline constexpr std::size_t kSpectralEmitBudgetBytes = 64u << 20;

// Below this a task is not worth its share of a dispatch, so the work is done on the calling thread
// instead. It is what every reduction hands PlanRowTasks as its floor, and all three of them had
// their own copy of the number -- two under one name and one under another, which is why the two
// looked like different rules rather than the same one written three times.
//
// A statement about pixels, whatever the thing being split is counted in: a plane histogram splits
// rows of a plane, a cube histogram splits rows across the planes of a read, and a spectral
// reduction splits chunk cells. PlanRowTasks multiplies out to pixels before it divides, so all
// three are asking the same question of the same number.
inline constexpr std::uint64_t kLeastPixelsPerTask = 1u << 16;

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_TUNING_H_
