/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// An allocation that fails inside a public entry point is a Result, as every other failure there is,
// and not a std::bad_alloc a consumer that checks Results never thought to catch.
//
// The global operator new is replaced so that one allocation, counted on the calling thread only,
// can be made to fail once. Only the calling thread's are counted because the library's own threads
// allocate whenever they like; and only once, because a report is a string and making it allocates.

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <new>
#include <string>
#include <typeinfo>

#include <carta-zarr/carta_zarr.h>

#include "support/check.h"

namespace {

// Allocations made on this thread since counting began, and which one of them is to fail; zero for
// none.
thread_local bool counting = false;
thread_local std::size_t allocations = 0;
thread_local std::size_t failing = 0;

void* Allocate(std::size_t size) {
    if (counting) {
        ++allocations;
        if (allocations == failing) {
            failing = 0;
            throw std::bad_alloc();
        }
    }
    if (void* memory = std::malloc(size == 0 ? 1 : size)) {
        return memory;
    }
    throw std::bad_alloc();
}

}  // namespace

void* operator new(std::size_t size) {
    return Allocate(size);
}
void* operator new[](std::size_t size) {
    return Allocate(size);
}
void operator delete(void* memory) noexcept {
    std::free(memory);
}
void operator delete[](void* memory) noexcept {
    std::free(memory);
}
void operator delete(void* memory, std::size_t) noexcept {
    std::free(memory);
}
void operator delete[](void* memory, std::size_t) noexcept {
    std::free(memory);
}

namespace {

using carta::zarr::testing::Require;

// How many allocations `call` makes on this thread when none fails.
template <typename Call>
std::size_t AllocationsOf(Call&& call) {
    allocations = 0;
    counting = true;
    call();
    counting = false;
    return allocations;
}

// Runs `call` with its `which`th allocation failing, and says what reached the caller.
template <typename Call>
std::string WithFailure(std::size_t which, Call&& call) {
    allocations = 0;
    failing = which;
    counting = true;
    std::string outcome;
    try {
        outcome = call() ? "a value" : "an error";
    } catch (const std::exception& error) {
        outcome = std::string("a thrown ") + typeid(error).name();
    }
    counting = false;
    failing = 0;
    return outcome;
}

void TestNewCachePool(const carta::zarr::Context& context) {
    const auto call = [&]() { return static_cast<bool>(context.NewCachePool(std::size_t{1} << 20)); };
    const auto total = AllocationsOf(call);
    Require(total > 0, "NewCachePool allocated nothing to fail");
    // The first, which is the library asking for the pool, and the last, which is the handle it
    // returns.
    for (const auto which : {std::size_t{1}, total}) {
        const auto outcome = WithFailure(which, call);
        Require(outcome == "an error", "NewCachePool with allocation " + std::to_string(which) + " of " +
                                           std::to_string(total) + " failing gave " + outcome);
    }
}

}  // namespace

int main() {
    auto context = carta::zarr::Context::Create({});
    Require(static_cast<bool>(context), "Context::Create failed");
    TestNewCachePool(*context);
    std::cout << "allocation failure tests passed\n";
    return 0;
}
