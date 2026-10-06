/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_TESTS_SUPPORT_CHECK_H_
#define CARTA_ZARR_TESTS_SUPPORT_CHECK_H_

// The assertion every test in this repository is written with.
//
// It was written out in twenty-four files, identically, because it is five lines and copying it
// was cheaper than reaching for it. Five lines times twenty-four is still the only thing every
// test here agrees about, so it is worth one place.
//
// A throw rather than an abort: a test's main catches it, names the suite, and returns 1, which is
// what lets one binary hold a suite rather than a case. The suites do that differently from each
// other -- some loop over fixtures and name the one that failed -- so this file stops at the
// assertion and does not try to be a runner as well.

#include <stdexcept>
#include <string>

namespace carta::zarr::testing {

inline void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace carta::zarr::testing

#endif  // CARTA_ZARR_TESTS_SUPPORT_CHECK_H_
