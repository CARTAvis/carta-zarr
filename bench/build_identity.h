/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_BENCH_BUILD_IDENTITY_H_
#define CARTA_ZARR_BENCH_BUILD_IDENTITY_H_

// Which build of the bench is measuring: what decides how it reads, as a run key has to say it.
//
// A trial is skipped on --resume when its key is in the CSV, and a sweep resumes when what it
// measures is unchanged; neither knew which build had measured. The commit a bench reports is the
// one at configure time, which says nothing of an edit since or of the library it loads -- and a
// reader changed since a trial ran is exactly what a rerun is for.
//
// So a build is the bytes it runs: the executable's and the carta-zarr library's it loaded, found
// where the loader found it. The same bytes linked again are the same build, and anything that
// changes how a read is done changes them. Read once per process, a few megabytes.

#include <string>

namespace carta::zarr::bench {

// Sixteen hex digits. A file that cannot be read is named in what is hashed instead of its bytes,
// so the answer is still a build's, if a coarser one.
std::string BuildIdentity();

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_BUILD_IDENTITY_H_
