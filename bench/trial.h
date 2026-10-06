/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_BENCH_TRIAL_H_
#define CARTA_ZARR_BENCH_TRIAL_H_

// One trial: the caches emptied, then one process per user, each with a context of its own, released
// together and timed against one clock.
//
// Every trial forks, even with one process. A process is what carta-controller gives each user, so a
// fresh one per trial is a fresh backend: no cache, no thread pool and no heap left over from the
// trial before. It is also why the parent never touches the library -- forking a process that has
// started TensorStore's threads leaves the child holding locks no thread will release.

#include "options.h"
#include "record.h"

#include <string>
#include <vector>

namespace carta::zarr::bench {

struct TrialOutcome {
    // Formatted CSV rows, process by process and operation by operation.
    std::vector<std::string> rows;
    unsigned ok = 0;
    unsigned timeouts = 0;
    unsigned errors = 0;
    // From the release to the last operation's end, over every process: what the slowest user waited.
    double makespan_s = 0.0;
};

TrialOutcome RunTrial(const RunOptions& options, Mode mode, unsigned trial, ColdMethod cold, const Row& base);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_TRIAL_H_
