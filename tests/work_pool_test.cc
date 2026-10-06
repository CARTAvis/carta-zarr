/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The library's own worker threads, and the rule that decides whether to wake them.
//
// This is the only concurrency carta-zarr owns, and the committed pixel fixtures are far too small
// to reach it -- PlanRowTasks answers 1 for a plane of twenty pixels, which is the correct answer.
// reduce_synthetic_test reaches the reductions' parallel paths with planes made up for the purpose;
// the pool itself, and the rules that decide how it is used, are tested here directly.

#include "work_pool.h"

#include "support/check.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using carta::zarr::internal::PlanRowTasks;
using carta::zarr::internal::TaskRows;
using carta::zarr::internal::WorkPool;

using carta::zarr::testing::Require;

void EveryTaskRunsExactlyOnce() {
    for (const std::size_t threads : {std::size_t{1}, std::size_t{2}, std::size_t{8}}) {
        WorkPool pool(threads);
        for (const std::size_t tasks : {std::size_t{0}, std::size_t{1}, std::size_t{3}, std::size_t{1000}}) {
            std::vector<int> seen(tasks, 0);
            pool.Run(tasks, [&](std::size_t task, std::size_t) { ++seen[task]; });
            for (std::size_t task = 0; task < tasks; ++task) {
                Require(seen[task] == 1, "task " + std::to_string(task) + " ran " + std::to_string(seen[task]) +
                                             " times with " + std::to_string(threads) + " threads");
            }
        }
    }
}

// The contract a private accumulator relies on: no two bodies running at the same moment are given
// the same worker index, and no index is ever outside [0, size()).
void WorkerIndicesAreExclusiveWhileRunning() {
    WorkPool pool(8);
    const std::size_t size = pool.size();
    std::vector<std::atomic<int>> occupied(size);
    for (auto& slot : occupied) {
        slot.store(0);
    }
    std::atomic<bool> collided{false};
    std::atomic<bool> out_of_range{false};

    pool.Run(4096, [&](std::size_t, std::size_t worker) {
        if (worker >= size) {
            out_of_range.store(true);
            return;
        }
        if (occupied[worker].fetch_add(1) != 0) {
            collided.store(true);
        }
        // Long enough that two bodies sharing a slot would overlap here rather than pass in turn.
        for (volatile int spin = 0; spin < 2000; ++spin) {
        }
        occupied[worker].fetch_sub(1);
    });

    Require(!out_of_range.load(), "a body was given a worker index outside [0, size())");
    Require(!collided.load(), "two concurrent bodies shared a worker index");
}

// A pool of one starts no threads, so the body has to run on the caller's thread or not at all.
void SingleThreadedPoolRunsInline() {
    WorkPool pool(1);
    Require(pool.size() == 1, "a pool of one reports a size other than one");
    const auto caller = std::this_thread::get_id();
    bool elsewhere = false;
    pool.Run(16, [&](std::size_t, std::size_t worker) {
        if (std::this_thread::get_id() != caller || worker != 0) {
            elsewhere = true;
        }
    });
    Require(!elsewhere, "a pool of one ran a body off the calling thread");
}

// Runs are sequential and reusable: the second one must not inherit the first one's counters.
void RepeatedRunsDoNotLeak() {
    WorkPool pool(4);
    std::atomic<std::uint64_t> total{0};
    for (int round = 0; round < 200; ++round) {
        total.store(0);
        pool.Run(97, [&](std::size_t task, std::size_t) { total.fetch_add(task); });
        Require(total.load() == (96U * 97U) / 2, "round " + std::to_string(round) + " summed wrong");
    }
}

// Two callers at once.
//
// The pool belongs to a Context, and a Context is shared for the whole process, so a cube histogram
// on one thread and a region profile on another really do arrive here together. The pool holds the
// body, the task counter and the running count in one set of fields; before those calls were
// serialised, the second caller overwrote them while the first caller's workers were still reading
// them, and a worker would run one body with the other's task count.
//
// Each job checks its own arithmetic: every task of its own run happens exactly once, and no task
// index arrives that its own run did not ask for.
void ConcurrentRunsDoNotShareState() {
    WorkPool pool(8);
    std::atomic<bool> wrong_task_index{false};
    std::atomic<bool> wrong_count{false};
    std::atomic<bool> uneven{false};

    const auto job = [&](std::size_t tasks, int rounds) {
        for (int round = 0; round < rounds; ++round) {
            std::vector<std::atomic<int>> seen(tasks);
            for (auto& count : seen) {
                count.store(0);
            }
            std::atomic<std::size_t> ran{0};
            pool.Run(tasks, [&](std::size_t task, std::size_t) {
                if (task >= tasks) {
                    // A task index this run never handed out: it came from the other call's count.
                    wrong_task_index.store(true);
                    return;
                }
                seen[task].fetch_add(1);
                ran.fetch_add(1);
            });
            if (ran.load() != tasks) {
                wrong_count.store(true);
            }
            for (const auto& count : seen) {
                if (count.load() != 1) {
                    uneven.store(true);
                }
            }
        }
    };

    // Unserialised, the first thing that goes wrong is not a wrong number: a call clobbers the
    // count of workers still out there and the other call waits for it forever. A test that hangs
    // reports nothing, so the wait is given a deadline it will never reach when this works -- the
    // whole file runs in hundredths of a second.
    std::thread watchdog([] {
        std::this_thread::sleep_for(std::chrono::seconds(30));
        std::cerr << "work pool test failed: concurrent runs deadlocked\n";
        std::_Exit(1);
    });
    watchdog.detach();

    // Different task counts, so a body running under the other call's count is out of range one way
    // and short the other.
    std::thread first(job, std::size_t{37}, 400);
    std::thread second(job, std::size_t{53}, 400);
    first.join();
    second.join();

    Require(!wrong_task_index.load(), "a body was given a task index its own run never asked for");
    Require(!wrong_count.load(), "a run executed a different number of bodies than it had tasks");
    Require(!uneven.load(), "a task ran other than exactly once");
}

void TheSplitRuleIsConservative() {
    // Nothing to split.
    Require(PlanRowTasks(1000, 0, 8, 1U << 16U) == 0, "no rows should ask for no tasks");
    // A plane smaller than one task's worth stays in place. This is the pixel fixture's case.
    Require(PlanRowTasks(5, 4, 28, 1U << 16U) == 1, "a twenty-pixel plane should not be split");
    // A pool of one never splits, however large the plane.
    Require(PlanRowTasks(4742, 7763, 1, 1U << 16U) == 1, "a pool of one should not split");
    // Never more pieces than rows.
    Require(PlanRowTasks(1U << 20U, 3, 28, 1U << 16U) == 3, "a three-row plane should give three tasks");
    // Never more pieces than the pool.
    Require(PlanRowTasks(4742, 7763, 28, 1U << 16U) == 28, "a large plane should fill the pool");
    // The affordability cap bites between those two.
    Require(PlanRowTasks(64, 64, 28, 1U << 16U) == 1, "4096 pixels should not be split 28 ways");
    Require(PlanRowTasks(256, 512, 28, 1U << 16U) == 2, "131072 pixels should give two tasks");
}

// What the histograms do with the plan: contiguous pieces, in task order, that cover every row
// exactly once. Asked of TaskRows, which is what they call, rather than of a copy of its arithmetic.
void RowRangesCoverEveryRowOnce() {
    for (std::uint64_t rows : {std::uint64_t{1}, std::uint64_t{7}, std::uint64_t{28}, std::uint64_t{7763}}) {
        for (const std::size_t tasks : {std::size_t{1}, std::size_t{3}, PlanRowTasks(1U << 20U, rows, 28, 1U << 16U)}) {
            Require(tasks >= 1, "a plan of zero tasks for a non-empty plane");
            std::vector<int> seen(rows, 0);
            std::uint64_t next = 0;
            for (std::size_t task = 0; task < tasks; ++task) {
                const auto range = TaskRows(task, tasks, rows);
                Require(range.first <= range.last && range.last <= rows, "a range stays inside the plane");
                if (range.first == range.last) {
                    continue;
                }
                Require(range.first == next, "ranges follow one another in task order");
                next = range.last;
                for (std::uint64_t row = range.first; row < range.last; ++row) {
                    ++seen[row];
                }
            }
            Require(next == rows, "the last range ends at the last row");
            for (std::uint64_t row = 0; row < rows; ++row) {
                Require(seen[row] == 1, "row " + std::to_string(row) + " of " + std::to_string(rows) + " was covered " +
                                            std::to_string(seen[row]) + " times");
            }
        }
    }
    // More tasks than rows leaves the ones past the end with nothing, not with a row past the end.
    Require(TaskRows(4, 5, 3).first == TaskRows(4, 5, 3).last, "a task past the last row is empty");
}

// A body that throws -- a histogram growing its bins is an allocation, and an allocation can fail
// -- reaches the thread that called Run, after every worker has stopped touching the caller's
// buffers. Left to propagate on a pool thread it would end the process: nothing above a
// std::thread's function catches.
void AThrowingBodyReachesTheCaller() {
    WorkPool pool(4);
    const auto caller = std::this_thread::get_id();
    for (const bool on_caller : {false, true}) {
        std::atomic<int> running{0};
        std::atomic<int> after_failure{0};
        std::atomic<bool> failed{false};
        bool caught = false;
        try {
            pool.Run(1000, [&](std::size_t, std::size_t) {
                ++running;
                if (failed.load()) {
                    ++after_failure;
                }
                const bool here = (std::this_thread::get_id() == caller) == on_caller;
                if (here && !failed.exchange(true)) {
                    --running;
                    throw std::bad_alloc();
                }
                std::this_thread::sleep_for(std::chrono::microseconds(20));
                --running;
            });
        } catch (const std::bad_alloc&) {
            caught = true;
        }
        const std::string where = on_caller ? "the calling thread" : "a pool thread";
        Require(caught, "an exception thrown on " + where + " did not reach the caller of Run");
        Require(running.load() == 0, "Run returned while a body thrown alongside was still running");
        // Bodies already claimed may finish, but the failure stops the rest being handed out.
        Require(after_failure.load() < 100, "tasks kept being started after one failed on " + where);

        // The pool is as it was: the next Run runs every task.
        std::atomic<std::size_t> ran{0};
        pool.Run(1000, [&](std::size_t, std::size_t) { ++ran; });
        Require(ran.load() == 1000, "a pool whose body threw did not run every task of the next Run");
    }
}

}  // namespace

int main() {
    try {
        EveryTaskRunsExactlyOnce();
        WorkerIndicesAreExclusiveWhileRunning();
        SingleThreadedPoolRunsInline();
        RepeatedRunsDoNotLeak();
        ConcurrentRunsDoNotShareState();
        TheSplitRuleIsConservative();
        RowRangesCoverEveryRowOnce();
        AThrowingBodyReachesTheCaller();
    } catch (const std::exception& error) {
        std::cerr << "work pool test failed: " << error.what() << "\n";
        return 1;
    }
    std::cout << "work pool tests passed\n";
    return 0;
}
