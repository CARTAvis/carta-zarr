/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_DEVIATIONS_H_
#define CARTA_ZARR_SRC_REDUCE_DEVIATIONS_H_

// The sum of squared deviations from the mean, Statistic::sum_sq_dev: how one set of pixels makes
// it, and how two sets' put together. See ADR 0018.
//
// Made in two steps because the per-pixel loop cannot afford a division per pixel, which is what
// updating a running mean costs. Within one span the loop sums the pixels' distances from a shift,
// one of the span's own pixels, and their squares -- one subtraction and one multiply-add more than
// it already does, and still vectorised. A span's pixels are close to that shift against their
// spread, so the subtraction that turns those sums into a sum of squared deviations loses next to
// nothing. Everything above a span -- the spans of a row, the rows of a region, the partials of a
// task, the reads of a block -- is put together pairwise as a Spread, which subtracts nothing large.

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace carta::zarr::internal {

// The value a span's distances are taken from: its first finite pixel, or zero when it has none,
// which leaves nothing to take distances of. Any of the span's own pixels would do; the first is the
// one found soonest.
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

// How many pixels, their mean, and the sum of their squared deviations from it: as much as putting
// two sets of pixels together needs.
//
// The mean is kept, and kept as a base and an offset from it, rather than had from a sum over a count
// or held as one double. Both of those carry a mean of pixels near 1e7 only to the 1.9e-9 a double
// resolves there, and the formula below multiplies the error in the difference of two means by the
// size of the sets: on a 256 x 260 plane of 1e7 with a spread of 0.5, a mean had from the sum left
// sum_sq_dev 1.3e-10 out and one held as a double 5e-11; as a base and an offset, 7e-16. The base
// is the first span's shift, one of the pixels, so two bases are near each other and their difference
// is exact; the offset is small, and resolved to the last digit.
struct Spread {
    double count = 0.0;
    double base = 0.0;
    double offset = 0.0;
    double sum_sq_dev = 0.0;

    // A span's, from its pixels' distances to `shift`: `count` of them, summing to `distance_sum`,
    // their squares to `distance_sum_sq`.
    //
    // Never below zero. It cannot be in exact arithmetic, and rounding can take it only a few units
    // of the last place below, for pixels that are all one value -- whose answer is zero.
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
