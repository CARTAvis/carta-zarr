/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_TASK_SPLIT_H_
#define CARTA_ZARR_SRC_REDUCE_TASK_SPLIT_H_

#include "reduce/tuning.h"
#include "work_pool.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

/**
 * How a reduction divides one read's arithmetic among the pool, given what each task's private
 * accumulator costs.
 *
 * Every reduction gives each task an accumulator of its own and adds them up afterwards, and every
 * one of them worked out the same three things for itself: how many accumulators it could afford,
 * which is a budget divided by what one costs and never more than the pool can run; how many tasks
 * a read was worth, which PlanRowTasks answers against kLeastPixelsPerTask; and how to hand each
 * task its contiguous run of rows. Three copies of the first, with the budget a constexpr inside the
 * function, is how the numbers ADR 0005 turns on came to live in three places.
 *
 * The accumulators themselves stay with the reductions, because they live for different spans: a
 * spectral reduction's are reset every slab and merged in order, a plane histogram's are zeroed per
 * plane and summed, and a cube histogram's last the whole walk and meet only at the end.
 *
 * What this does keep is ADR 0005's rule that there are never more tasks than accumulators. `Tasks`
 * never answers more than `most`, and a reduction reaches the pool only through `Run`, so a body
 * indexing its accumulator by task never shares one with another body.
 *
 * `Run` takes its body as a template parameter, and when a read is one task it calls it directly on
 * the calling thread. The pool's own std::function is paid once per task, outside the pixel loop,
 * which is the arrangement ADR 0005 records and the reason this must not grow a type-erased entry
 * point.
 *
 * Cheap to make -- a reference and a count -- so a reduction whose accumulator's size depends on the
 * read makes one per read.
 */
class TaskSplit {
public:
    // `budget_bytes` is what every task's accumulator may cost together, and `accumulator_bytes` is
    // what one costs; zero is taken as one byte. See tuning.h for the budgets and why they differ.
    TaskSplit(WorkPool& workers, std::size_t budget_bytes, std::size_t accumulator_bytes)
        : _workers(&workers),
          _most(std::min(workers.size(),
                         std::max<std::size_t>(1, budget_bytes / std::max<std::size_t>(1, accumulator_bytes)))) {}

    // The most tasks any read will be split into, which is how many accumulators a reduction that
    // allocates them once has to hold. At least one.
    std::size_t most() const noexcept { return _most; }

    // How many tasks `units` units of `unit_pixels` pixels each are worth: one when they are too few
    // to be worth waking anyone for, and never more than `most`.
    std::size_t Tasks(std::uint64_t unit_pixels, std::uint64_t units) const {
        return PlanRowTasks(unit_pixels, units, _most, kLeastPixelsPerTask);
    }

    // Call `body(task, first, last)` for each of `tasks` contiguous runs of [0, units), and return
    // once all of them are done. `tasks` comes from Tasks. One task runs on the calling thread
    // without going near the pool; a task whose run is empty is not called at all.
    template <typename Body>
    void Run(std::size_t tasks, std::uint64_t units, const Body& body) const {
        if (tasks <= 1) {
            body(std::size_t{0}, std::uint64_t{0}, units);
            return;
        }
        _workers->Run(tasks, [&](std::size_t task, std::size_t) {
            const auto share = TaskRows(task, tasks, units);
            if (share.first == share.last) {
                return;
            }
            body(task, share.first, share.last);
        });
    }

private:
    WorkPool* _workers;
    std::size_t _most;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_TASK_SPLIT_H_
