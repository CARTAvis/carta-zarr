# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# Run with cmake -P in TensorStore's source tree, as its PATCH_COMMAND, when CARTA_ZARR_THREAD_SANITIZER
# is on.
#
# TensorStore hands a task to another worker through a work-stealing queue whose push publishes it
# with a release fence and then a relaxed store of the queue's bottom, which a thief loads with
# acquire. ThreadSanitizer does not model a fence on its own, so it sees no ordering between the
# thread that queued a task and the one that stole it, and reports everything the task touches as
# raced: dozens of reports a run, inside TensorStore, and none of them real.
#
# Making that store a release store says the same thing in a form it follows: the fence already orders
# everything before it ahead of the store, and the store now carries that order itself. With it, the
# threading tests run with no report and no suppression. Only a ThreadSanitizer build is patched; every
# other build compiles TensorStore as released.

set(queue "tensorstore/internal/container/single_producer_queue.h")
set(fenced "    std::atomic_thread_fence(std::memory_order_release);\n    bottom_.store(b + 1, std::memory_order_relaxed);")
set(released "    std::atomic_thread_fence(std::memory_order_release);\n    bottom_.store(b + 1, std::memory_order_release);")

file(READ "${queue}" source)
string(FIND "${source}" "${released}" already)
if(already GREATER_EQUAL 0)
    return()
endif()
string(FIND "${source}" "${fenced}" at)
if(at LESS 0)
    message(FATAL_ERROR "${queue} no longer publishes a task the way this patch expects; read it again "
                        "before building TensorStore for ThreadSanitizer")
endif()
string(REPLACE "${fenced}" "${released}" source "${source}")
file(WRITE "${queue}" "${source}")
