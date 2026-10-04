/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_GROWING_HISTOGRAM_H_
#define CARTA_ZARR_SRC_REDUCE_GROWING_HISTOGRAM_H_

// The one-pass histogram's arithmetic, on its own.
//
// It was already a deep module -- two methods over a range that seeds itself, doubles, merges and
// re-aggregates -- but the seam was file-local inside plane_histogram.cc, so the only way to reach
// it was ComputeCubeHistogram over a fixture on disk, compared against a two-pass oracle to a
// tolerance. It is the most numerically subtle code in this library and it has had three fixes --
// 9c21397 on how a straddling bin is split, 6ba0ce2 and then a narrower seed on where the range
// seeds itself -- and each was a five-line exact-count test from here.
//
// Nothing about the arithmetic changed in the move, and ComputeCubeHistogram uses the class exactly
// as before. One thing about its codegen did, and it is recorded here rather than left to be
// rediscovered: `Grow` has to stay out of line.
//
// The class was file-local in an anonymous namespace, which told the compiler that nothing outside
// that translation unit could call `Add`, and it inlined `Add` into the per-pixel loop. Given
// external linkage it does the opposite -- it folds the cold `Grow`, which allocates a vector, into
// `Add`, and the result is too large to inline into the loop that calls it a billion times. Marking
// `Grow` as not-inline restores the choice the anonymous namespace used to make for free.
//
// Measured on build-release against tests/data/.../pixels_wide (512 x 520 x 4, four threads, nine
// repeats), best median ComputeCubeHistogram of three runs: 1.49 ms file-local, 1.82 ms moved here
// as-is, 1.38 ms moved here with `Grow` out of line. As ADR 0005 warns about the visitor, none of
// this is visible in a Debug build or in any test that asserts.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace carta::zarr::internal {

// A histogram whose range grows to fit whatever arrives.
//
// Doubling to one side and merging bins in pairs keeps every count: the old range becomes one half
// of the new one, and old bins 2m and 2m+1 become bin m of that half. Resolution halves each time,
// which is the only thing given up -- and only for as many doublings as the data's dynamic range
// actually needs.
class GrowingHistogram {
public:
    explicit GrowingHistogram(std::size_t bins) : _counts(bins, 0) {}

    void Add(float value) {
        const double v = value;
        if (!_seeded) {
            // A first range around the first pixel seen, as narrow as a range can usefully be: one
            // bin per float, at the spacing of the floats there. The range only ever grows, so a
            // guess that is too wide is permanent while one that is too narrow costs a few merges and
            // then fits -- each merge doubles it, so reaching any spread the pixels have is a few
            // dozen, once, against a cube of billions of pixels.
            //
            // It was seeded twice too wide before this. Anchored at one, it spent almost all of the
            // resolution on the empty space between a Jansky and the hundredths of one the pixels
            // actually are. Anchored at the pixel's own magnitude, it was 2000 wide around a first
            // pixel of 1000, and pixels a sixteenth apart there shared one provisional bin and came
            // out spread evenly over every target bin between them.
            //
            // A first pixel of exactly zero seeds at the least denormal, which is as narrow as it
            // gets: about a hundred and fifty merges to reach a spread of one. A range of zero width
            // would never grow at all.
            const auto at = static_cast<float>(std::abs(v));
            const double spacing =
                static_cast<double>(std::nextafter(at, std::numeric_limits<float>::infinity())) - static_cast<double>(at);
            _width = spacing;
            _lower = v - (spacing * static_cast<double>(_counts.size() / 2));
            _seeded = true;
        }
        while (v < _lower) {
            Grow(false);
        }
        while (v >= _lower + (_width * static_cast<double>(_counts.size()))) {
            Grow(true);
        }
        auto bin = static_cast<std::size_t>((v - _lower) / _width);
        if (bin >= _counts.size()) {
            bin = _counts.size() - 1;
        }
        ++_counts[bin];
    }

    // Re-aggregate over the range the caller wants, splitting each provisional bin between the
    // target bins it overlaps in proportion to how much of it each one covers.
    //
    // Giving the whole of it to the bin its centre falls in is the obvious thing and it is what
    // this used to do, but the error that leaves is a bias rather than a wobble: every provisional
    // bin straddling a target edge leans the same way, and nothing averages it out. Splitting
    // assumes the pixels are spread evenly inside one provisional bin, which is the assumption the
    // caller's own percentile already makes between target bins.
    //
    // Largest remainder, so the split still adds up: every pixel the walk binned comes out in some
    // target bin, which is what lets the caller compare the total against its own pixel count.
    std::vector<std::uint64_t> Aggregate(std::size_t bins, double lower, double upper) const {
        std::vector<std::uint64_t> out(bins, 0);
        std::uint64_t total = 0;
        for (const auto count : _counts) {
            total += count;
        }
        if (!(lower < upper)) {
            // Every pixel had the same value, so there is one bin it can be in.
            out.front() = total;
            return out;
        }
        const double target_width = (upper - lower) / static_cast<double>(bins);
        const auto last_bin = static_cast<std::ptrdiff_t>(bins) - 1;
        // Held outside the loop: a provisional bin usually overlaps two target bins, and this would
        // otherwise allocate for every one of tens of thousands of them.
        std::vector<double> shares;
        for (std::size_t i = 0; i < _counts.size(); ++i) {
            const std::uint64_t count = _counts[i];
            if (count == 0) {
                continue;
            }
            const double from = _lower + (_width * static_cast<double>(i));
            const double to = from + _width;
            auto first = static_cast<std::ptrdiff_t>(std::floor((from - lower) / target_width));
            auto last = static_cast<std::ptrdiff_t>(std::floor((to - lower) / target_width));
            first = std::clamp<std::ptrdiff_t>(first, 0, last_bin);
            last = std::clamp<std::ptrdiff_t>(last, 0, last_bin);
            if (first == last) {
                out[static_cast<std::size_t>(first)] += count;
                continue;
            }

            // Normalised by what the target range actually covers, not by the provisional width: a
            // bin hanging over either end would otherwise leave a remainder the size of the part
            // outside, and those pixels are inside the range by construction.
            shares.assign(static_cast<std::size_t>(last - first + 1), 0.0);
            double covered = 0.0;
            for (std::ptrdiff_t bin = first; bin <= last; ++bin) {
                const double low = std::max(from, lower + (target_width * static_cast<double>(bin)));
                const double high = std::min(to, lower + (target_width * static_cast<double>(bin + 1)));
                const double piece = high > low ? high - low : 0.0;
                shares[static_cast<std::size_t>(bin - first)] = piece;
                covered += piece;
            }
            if (!(covered > 0.0)) {
                out[static_cast<std::size_t>(first)] += count;
                continue;
            }

            std::uint64_t placed = 0;
            for (std::size_t offset = 0; offset < shares.size(); ++offset) {
                const double share = static_cast<double>(count) * (shares[offset] / covered);
                const auto whole = static_cast<std::uint64_t>(share);
                out[static_cast<std::size_t>(first) + offset] += whole;
                placed += whole;
                // Reused as the fractional part, which is what decides who gets the leftovers.
                shares[offset] = share - static_cast<double>(whole);
            }
            while (placed < count) {
                std::size_t best = 0;
                for (std::size_t offset = 1; offset < shares.size(); ++offset) {
                    if (shares[offset] > shares[best]) {
                        best = offset;
                    }
                }
                ++out[static_cast<std::size_t>(first) + best];
                shares[best] = -1.0;
                ++placed;
            }
        }
        return out;
    }

private:
    // Not inlined, and that is load-bearing rather than a hint: see the head of this file. It runs
    // once per doubling of the range against a billion calls to Add, so nothing is lost by the
    // call.
    [[gnu::noinline]] void Grow(bool upward) {
        const std::size_t half = _counts.size() / 2;
        std::vector<std::uint64_t> merged(_counts.size(), 0);
        for (std::size_t m = 0; m < half; ++m) {
            merged[upward ? m : half + m] = _counts[2 * m] + _counts[(2 * m) + 1];
        }
        if (!upward) {
            _lower -= _width * static_cast<double>(_counts.size());
        }
        _width *= 2.0;
        _counts.swap(merged);
    }

    std::vector<std::uint64_t> _counts;
    double _lower = 0.0;
    double _width = 1.0;
    bool _seeded = false;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_GROWING_HISTOGRAM_H_
