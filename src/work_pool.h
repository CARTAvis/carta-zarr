/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_WORK_POOL_H_
#define CARTA_ZARR_SRC_WORK_POOL_H_

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace carta::zarr::internal {

/**
 * The library's own worker threads, for the per-pixel work a reduction does after a read returns.
 *
 * This is the only concurrency carta-zarr owns. Everything else is TensorStore's: its read pool
 * decompresses chunks, and until this existed the pixels those chunks produced were then visited on
 * one thread, which measured as 78% of a cube scan on a 28-core machine while the decode pool sat
 * idle (implementation-plan.md 6.2.3 and 6.2.5).
 *
 * Deliberately small, and deliberately not a general executor:
 *
 * - Threads are started once and parked, because the alternative -- starting them per slab -- costs
 *   more than the work when a slab is a few hundred thousand pixels.
 * - `Run` blocks until every worker has finished, so a caller's buffers stay alive for the whole
 *   span and there is no future to own. The walk is a loop over slabs, and it wants each slab done
 *   before it reuses the buffer for the next one.
 * - The body is handed a worker index as well as a task index, so a caller that needs private
 *   accumulation can address one slot per worker without a map or a lock.
 * - An exception a body throws reaches the caller of `Run`, rethrown there once every worker has
 *   stopped. A body is mostly arithmetic over a buffer the caller owns, but not only: a histogram
 *   growing its bins allocates, and an allocation can fail. Thrown on a pool thread and left there,
 *   it ended the process, since nothing above a std::thread's function catches; on the caller's
 *   thread it reaches the public entry point, which reports it as an error like any other.
 */
class WorkPool {
public:
    // `threads` of zero asks for one worker per hardware thread. One means everything runs inline
    // on the calling thread and no threads are started at all.
    //
    // Throws std::system_error if the system refuses a thread. Nothing is left running when it
    // does, and the caller gets no pool rather than a smaller one: a pool that quietly started
    // three of thirty-two workers would report the difference only as everything being slow.
    explicit WorkPool(std::size_t threads);
    ~WorkPool();

    WorkPool(const WorkPool&) = delete;
    WorkPool& operator=(const WorkPool&) = delete;

    // How many bodies may run at once, counting the calling thread. Always at least one.
    std::size_t size() const noexcept {
        return _workers.empty() ? 1 : _workers.size() + 1;
    }

    /**
     * Run `body(task, worker)` for every task in [0, tasks), and return once all of them are done.
     *
     * `worker` is in [0, size()) and no two bodies running at the same time share one, so a private
     * accumulator indexed by it needs no lock. Tasks are claimed from a shared counter rather than
     * divided up front, so an uneven one does not leave workers waiting on the slowest slice.
     *
     * Safe to call from more than one thread: calls are serialised, and one of them waits. That is
     * not a detail a caller can ignore on a pool it does not own -- see `_run` below.
     *
     * If a body throws, no task not yet claimed is started, the ones already running finish, and
     * the first exception is rethrown here. Which tasks ran is then unspecified, so a caller
     * treats whatever they were writing as lost -- which a caller unwinding does anyway.
     */
    void Run(std::size_t tasks, const std::function<void(std::size_t task, std::size_t worker)>& body);

private:
    void Worker(std::size_t index);
    // Take tasks from the shared counter until there are none left, running each one as `worker`.
    // Both a pool thread and the thread that called Run do exactly this, which is what makes the
    // caller worker 0 rather than a special case.
    void DrainTasks(const std::function<void(std::size_t task, std::size_t worker)>& body, std::size_t worker);
    // Wake every worker, tell it to return, and join it. Idempotent, and the only way a worker ever
    // ends -- one that is parked in _wake has nothing else to wake it.
    void StopWorkers() noexcept;

    std::vector<std::thread> _workers;

    // One Run at a time. The body, the task counter and the count of workers still out there are
    // one set of fields, so a second caller arriving while the first is still running would hand
    // that first caller's workers its own body and its own task count -- a body indexing buffers
    // sized for the other call.
    //
    // It is not a hypothetical: the pool belongs to a Context, a Context is shared for the whole
    // process, and the caller runs a cube histogram on one thread while a region profile runs on
    // another. Held for one call only, and a call is the arithmetic over a single read, so the
    // waiting is bounded by that and not by the walk.
    //
    // A body must not call Run. None does -- they are arithmetic over a buffer the caller already
    // owns -- and this is the reason to keep it that way.
    std::mutex _run;

    std::mutex _mutex;
    std::condition_variable _wake;
    std::condition_variable _done;

    const std::function<void(std::size_t, std::size_t)>* _body = nullptr;
    std::size_t _tasks = 0;
    std::size_t _next = 0;
    std::size_t _running = 0;
    std::size_t _generation = 0;
    bool _stopping = false;
    // The first exception a body of the current Run threw, for Run to rethrow.
    std::exception_ptr _failure;
};

/**
 * How many contiguous pieces a plane of `rows` rows of `row_pixels` each should be split into.
 *
 * Separate from the pool, and pure, because it is the half of the split that has an answer worth
 * checking: a piece smaller than `least_pixels` costs more to hand out than to run, and a plane
 * cannot be split into more pieces than it has rows. Returning 1 means the caller should bin in
 * place and not go near the pool at all -- which is the right answer for a tiny image, and is why
 * a fixture of a few dozen pixels never reaches the parallel path.
 */
std::size_t PlanRowTasks(std::uint64_t row_pixels, std::uint64_t rows, std::size_t max_tasks,
                         std::uint64_t least_pixels);

// The rows [first, last) that task `task` of `tasks` takes, when `rows` rows are cut into contiguous
// pieces in task order. Every row is in exactly one piece; a task past the last row gets an empty
// one, first == last, rather than a row that is not there. `tasks` is at least one.
//
// One rule rather than a copy per caller, because where the pieces fall decides the order a task
// adds up its rows in -- so a caller that cut them differently would agree with the others on
// everything but the last bits of a floating-point sum.
struct RowRange {
    std::uint64_t first = 0;
    std::uint64_t last = 0;
};
RowRange TaskRows(std::size_t task, std::size_t tasks, std::uint64_t rows);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_WORK_POOL_H_
