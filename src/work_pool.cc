/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "work_pool.h"

#include <algorithm>
#include <utility>

namespace carta::zarr::internal {

WorkPool::WorkPool(std::size_t threads) {
    if (threads == 0) {
        const auto hardware = std::thread::hardware_concurrency();
        threads = hardware > 0 ? static_cast<std::size_t>(hardware) : 1;
    }
    // The calling thread is one of the workers, so it is the one that is not started here. A pool
    // of one starts nothing and runs everything inline, which is what a consumer asking for a
    // single-threaded library gets.
    _workers.reserve(threads - 1);
    try {
        for (std::size_t index = 1; index < threads; ++index) {
            _workers.emplace_back([this, index] { Worker(index); });
        }
    } catch (...) {
        // A constructor that throws never runs its own destructor, but its members are destroyed
        // anyway -- and ~thread on a joinable thread is std::terminate, not an error. So the
        // workers already started have to be stopped here or the process dies during unwinding,
        // taking the consumer with it. They would also never end on their own: each one is parked
        // in _wake, and nothing else is going to set _stopping.
        StopWorkers();
        throw;
    }
}

WorkPool::~WorkPool() {
    StopWorkers();
}

void WorkPool::StopWorkers() noexcept {
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _stopping = true;
    }
    _wake.notify_all();
    for (auto& worker : _workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    // Leaves nothing joinable behind, which is the postcondition the constructor's failure path
    // needs and the destructor gets for free.
    _workers.clear();
}

void WorkPool::DrainTasks(const std::function<void(std::size_t, std::size_t)>& body, std::size_t worker) {
    while (true) {
        std::size_t task = 0;
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (_next >= _tasks) {
                return;
            }
            task = _next++;
        }
        // Outside the lock: a body is the arithmetic a reduction came for, and holding the counter
        // across it would serialise the very thing the pool exists to spread.
        try {
            body(task, worker);
        } catch (...) {
            // Kept for Run to rethrow on the caller's thread, and every task not yet claimed given
            // up: the call has failed, and the sooner the workers stop, the sooner it says so.
            std::lock_guard<std::mutex> lock(_mutex);
            if (!_failure) {
                _failure = std::current_exception();
            }
            _next = _tasks;
            return;
        }
    }
}

void WorkPool::Worker(std::size_t index) {
    std::size_t seen = 0;
    while (true) {
        const std::function<void(std::size_t, std::size_t)>* body = nullptr;
        {
            std::unique_lock<std::mutex> lock(_mutex);
            _wake.wait(lock, [&] { return _stopping || _generation != seen; });
            if (_stopping) {
                return;
            }
            seen = _generation;
            body = _body;
        }

        DrainTasks(*body, index);

        {
            std::lock_guard<std::mutex> lock(_mutex);
            --_running;
            if (_running == 0) {
                _done.notify_one();
            }
        }
    }
}

void WorkPool::Run(std::size_t tasks, const std::function<void(std::size_t, std::size_t)>& body) {
    if (tasks == 0) {
        return;
    }
    if (_workers.empty() || tasks == 1) {
        // Nothing to gain from waking anyone: the caller's thread is a worker too, and it is
        // already here.
        for (std::size_t task = 0; task < tasks; ++task) {
            body(task, 0);
        }
        return;
    }

    // Everything below writes the pool's one set of fields and then waits for workers reading
    // them, so only one call may be inside it at a time. See _run.
    std::lock_guard<std::mutex> serialised(_run);
    std::exception_ptr failure;

    {
        std::lock_guard<std::mutex> lock(_mutex);
        _body = &body;
        _tasks = tasks;
        _next = 0;
        // Counted here rather than by each worker as it wakes, so that a Run which finds every
        // task already taken by the calling thread still waits for exactly the workers it woke.
        _running = _workers.size();
        ++_generation;
    }
    _wake.notify_all();

    // The calling thread takes tasks as worker 0 rather than waiting, which is what keeps a pool of
    // N threads worth N and not N-1.
    DrainTasks(body, 0);

    {
        std::unique_lock<std::mutex> lock(_mutex);
        _done.wait(lock, [&] { return _running == 0; });
        _body = nullptr;
        _tasks = 0;
        failure = std::exchange(_failure, nullptr);
    }
    if (failure) {
        std::rethrow_exception(failure);
    }
}

std::size_t PlanRowTasks(std::uint64_t row_pixels, std::uint64_t rows, std::size_t max_tasks,
                         std::uint64_t least_pixels) {
    if (rows == 0 || row_pixels == 0 || max_tasks <= 1) {
        return rows == 0 ? 0 : 1;
    }
    const std::uint64_t pixels = row_pixels * rows;
    const std::uint64_t affordable = std::max<std::uint64_t>(1, pixels / std::max<std::uint64_t>(1, least_pixels));
    return static_cast<std::size_t>(std::min<std::uint64_t>({static_cast<std::uint64_t>(max_tasks), affordable, rows}));
}

RowRange TaskRows(std::size_t task, std::size_t tasks, std::uint64_t rows) {
    const std::uint64_t per_task = (rows + tasks - 1) / tasks;
    const std::uint64_t first = static_cast<std::uint64_t>(task) * per_task;
    return RowRange{std::min(first, rows), std::min(first + per_task, rows)};
}

}  // namespace carta::zarr::internal
