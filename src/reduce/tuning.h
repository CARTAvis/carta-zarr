/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_TUNING_H_
#define CARTA_ZARR_SRC_REDUCE_TUNING_H_

// Numbers the reductions were tuned to, which a caller neither sets nor needs.
//
// They are not part of the public contract: nothing outside src/reduce/ reads them. The limits a
// caller does have to respect, kMaxHistogramBins and kMaxSpectralRegions, are in carta-zarr/reduce.h,
// because a request that exceeds either is rejected.
//
// Four of them were measured on one machine and depend on its caches and thread count rather than on
// what is being computed: how many pixels a task needs, and what the tasks' private accumulators may
// cost together. Each can be replaced when the library is built, through CARTA_ZARR_TUNING_OVERRIDES,
// so that carta-zarr-bench can measure another value on another machine. That is all it is for: an
// override is a measurement, not a setting a deployment chooses. See ADR 0014. The others decide what
// a reduction answers, such as how finely a cube histogram resolves a bin or how much one block may
// hold, and cannot be replaced.

#include "carta-zarr/reduce.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

// Every statistic, in the order a SpectralBlock lays them out.
inline constexpr std::array<Statistic, 7> kStatisticOrder{Statistic::num_pixels, Statistic::nan_count, Statistic::sum,
                                                          Statistic::sum_sq,     Statistic::min,       Statistic::max,
                                                          Statistic::sum_sq_dev};

// The provisional resolution a cube histogram bins at when the caller does not choose one.
//
// Sixteen provisional bins for every bin asked for, because the error depends on how finely the walk
// resolves one target bin, not on how many bins it holds in total. And a provisional histogram is
// eight bytes a bin, with one per worker, so unneeded bins cost cache. A thousand target bins get
// 16,384 provisional ones, which is 128 kB.
//
// Held between 4,096 and 65,536. On a billion-pixel cube against the exact two-pass answer, the
// error rose from 0.003% of pixels at 65,536 to 0.011% at 4,096 while the time fell from 3.64 s to
// 2.14 s, and below 4,096 both get worse at once. At sixteen the time has leveled off, with a few
// times the resolution the error needs.
inline constexpr std::uint32_t kProvisionalBinsPerBin = 16;
inline constexpr std::uint32_t kLeastProvisionalBins = 1u << 12;
inline constexpr std::uint32_t kMostProvisionalBins = 1u << 16;

// The memory one emitted block may occupy. emit_every_channels is reduced to fit it; see
// SpectralReduceRequest::emit_every_channels.
inline constexpr std::size_t kSpectralEmitBudgetBytes = 64u << 20;

// Below this a task is not worth its share of a dispatch, so the work is done on the calling thread
// instead. It is the floor TaskSplit hands PlanRowTasks for every reduction.
//
// It counts pixels, whatever unit the work is split in: a plane histogram splits rows of a plane, a
// cube histogram splits rows across the planes of a read, and a spectral reduction splits chunk
// cells. PlanRowTasks converts to pixels before it divides, so all three use the same number.
//
// At 65,536 a 512 x 512 plane would be four tasks however many cores there were. On a 512x512x7776
// cube, warm, one user, an exact cube histogram binned plane by plane took 11-17% less time at
// 16,384 than at 65,536 on both a desktop and a server, and far longer at 262,144; with eight users
// every value tried was within 6%. 4,096 gained a few percent more, which is less than it risks for
// a small read on a slower dispatch. ADR 0014 has the measurements.
#ifndef CARTA_ZARR_TUNING_LEAST_PIXELS_PER_TASK
#define CARTA_ZARR_TUNING_LEAST_PIXELS_PER_TASK (1u << 14)
#endif
inline constexpr std::uint64_t kLeastPixelsPerTask = CARTA_ZARR_TUNING_LEAST_PIXELS_PER_TASK;
static_assert(kLeastPixelsPerTask > 0, "CARTA_ZARR_TUNING_LEAST_PIXELS_PER_TASK must be positive");

// What all the tasks' private accumulators may cost together, one budget per reduction, each handed
// to TaskSplit with what one accumulator costs. There are three because they bound two different
// things: the first two bound memory, and only the third is ADR 0005's cache cap.
//
// A spectral reduction's partials are the region totals for one slab, allocated and zeroed once per
// slab. The budget stops a reduction over thousands of regions from spending more on the split than
// on the pixels.
#ifndef CARTA_ZARR_TUNING_SPECTRAL_PARTIAL_BUDGET_BYTES
#define CARTA_ZARR_TUNING_SPECTRAL_PARTIAL_BUDGET_BYTES (16u << 20)
#endif
inline constexpr std::size_t kSpectralPartialBudgetBytes = CARTA_ZARR_TUNING_SPECTRAL_PARTIAL_BUDGET_BYTES;
static_assert(kSpectralPartialBudgetBytes > 0, "CARTA_ZARR_TUNING_SPECTRAL_PARTIAL_BUDGET_BYTES must be positive");

// A plane histogram's partials are a copy of the bins per task. A caller may ask for as many as
// kMaxHistogramBins, and a private copy of that for every worker would be hundreds of megabytes for a
// pass that is meant to stream. This keeps a large bin count from turning a split into a large
// allocation.
#ifndef CARTA_ZARR_TUNING_HISTOGRAM_PARTIAL_BUDGET_BYTES
#define CARTA_ZARR_TUNING_HISTOGRAM_PARTIAL_BUDGET_BYTES (64u << 20)
#endif
inline constexpr std::size_t kHistogramPartialBudgetBytes = CARTA_ZARR_TUNING_HISTOGRAM_PARTIAL_BUDGET_BYTES;
static_assert(kHistogramPartialBudgetBytes > 0, "CARTA_ZARR_TUNING_HISTOGRAM_PARTIAL_BUDGET_BYTES must be positive");

// A cube histogram's accumulators are one provisional histogram each, and they are bounded by cache
// rather than memory: this is the cap ADR 0005 is about. A provisional histogram is eight bytes a bin
// (half a megabyte at the default resolution), and every worker writes to its own at random while
// streaming its share of the pixels through the same cache. Past a couple of megabytes of them, the
// pixels evict the histograms and the pass gets slower the more workers it uses.
//
// Measured on a 512x512x7776 ASKAP cube, warm, against 8.8 s for not splitting at all: four workers
// 4.9 s, eight 8.7 s, twenty-eight 16.5 s. At the default resolution this comes out at four.
#ifndef CARTA_ZARR_TUNING_CUBE_ACCUMULATOR_CACHE_BYTES
#define CARTA_ZARR_TUNING_CUBE_ACCUMULATOR_CACHE_BYTES (2u << 20)
#endif
inline constexpr std::size_t kCubeAccumulatorCacheBytes = CARTA_ZARR_TUNING_CUBE_ACCUMULATOR_CACHE_BYTES;
static_assert(kCubeAccumulatorCacheBytes > 0, "CARTA_ZARR_TUNING_CUBE_ACCUMULATOR_CACHE_BYTES must be positive");

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_TUNING_H_
