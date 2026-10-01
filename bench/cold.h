/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_BENCH_COLD_H_
#define CARTA_ZARR_BENCH_COLD_H_

// Emptying the caches between a trial and the storage, so that it reads from the disks.
//
// Only the administrator's command can reach a parallel filesystem's servers. drop_caches and
// fadvise empty this host's page cache and nothing beyond it, which on Lustre or BeeGFS leaves the
// servers' caches warm; the CSV records which was used so that the report can say so.

#include "options.h"

#include <optional>
#include <string>

namespace carta::zarr::bench {

// The method a run uses: the one asked for, or for auto the best this host allows. An error when the
// one asked for cannot work here.
struct ColdChoice {
    ColdMethod method = ColdMethod::off;
    std::string error;
};

ColdChoice ChooseColdMethod(std::optional<ColdMethod> asked, const std::string& command);

// Empties the caches by `method`. An empty string when it worked, otherwise what went wrong.
std::string DropCaches(ColdMethod method, const std::string& command, const std::string& dataset);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_COLD_H_
