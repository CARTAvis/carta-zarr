/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_DEVIATIONS_H_
#define CARTA_ZARR_SRC_REDUCE_DEVIATIONS_H_

// The sum of squared deviations from the mean, Statistic::sum_sq_dev: how one set of pixels makes
// it, and how two sets' are combined. See ADR 0018.
//
// Made in two steps because the per-pixel loop cannot afford a division per pixel, which is what
// updating a running mean costs. Within one span, the loop sums the pixels' distances from a shift
// (one of the span's own pixels) and their squares: one subtraction and one multiply-add more than
// it already does, and still vectorized. A span's pixels are close to that shift compared with their
// spread, so turning those sums into a sum of squared deviations loses almost nothing. Everything
// above a span (the spans of a row, the rows of a region, the partials of a task, the reads of a
// block) is combined pairwise as a Spread, which subtracts nothing large.

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace carta::zarr::internal {

// The value a span's distances are taken from: its first finite pixel, or zero when it has none.
// Any of the span's own pixels would do; the first is found soonest.
template <bool kMasked>
float DeviationShift(const float* row, std::uint64_t stride, std::uint64_t count, const std::uint8_t* selected,
                     std::uint64_t mask_stride) {
    for (std::uint64_t i = 0; i < count; ++i) {
        if (kMasked && selected[i * mask_stride] == 0) {
            continue;
        }
        const float value = row[i * stride];
        if (std::isfinite(value)) {
            return value;
        }
    }
    return 0.0F;
}

// How many pixels, their mean, and the sum of their squared deviations from it: all that combining
// two sets of pixels needs.
//
// The mean is kept as a base plus an offset, rather than computed as sum / count or held as one
// double. Either of those resolves a mean near 1e7 only to the 1.9e-9 a double resolves there, and
// the formula below multiplies the error in the difference of two means by the size of the sets. On a
// 256 x 260 plane of 1e7 with a spread of 0.5, sum_sq_dev was off by 1.3e-10 with the mean from the
// sum and by 5e-11 with one double; with a base and an offset, by 7e-16. The base is the first span's
// shift, one of the pixels, so two bases are close and their difference is exact; the offset is
// small and resolved to the last digit.
struct Spread {
    double count = 0.0;
    double base = 0.0;
    double offset = 0.0;
    double sum_sq_dev = 0.0;

    // A span's, from its pixels' distances to `shift`: `count` of them, summing to `distance_sum`,
    // their squares to `distance_sum_sq`.
    //
    // Never below zero. It cannot be in exact arithmetic, and rounding can take it a few units in the
    // last place below only when every pixel has the same value, where the answer is zero.
    static Spread OfSpan(double count, double shift, double distance_sum, double distance_sum_sq) {
        if (count <= 0.0) {
            return {};
        }
        return {count, shift, distance_sum / count,
                std::max(distance_sum_sq - (distance_sum * distance_sum / count), 0.0)};
    }

    // Takes in another set of pixels, by the formula of Chan, Golub and LeVeque. The term that
    // accounts for the two means differing is the square of their difference, which is small exactly
    // when it matters little. The base stays where the first set put it.
    void Merge(const Spread& other) {
        if (other.count <= 0.0) {
            return;
        }
        if (count <= 0.0) {
            *this = other;
            return;
        }
        const double total = count + other.count;
        const double between = (other.base - base) + (other.offset - offset);
        sum_sq_dev += other.sum_sq_dev + (between * between * (count * other.count / total));
        offset += between * (other.count / total);
        count = total;
    }
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_DEVIATIONS_H_
