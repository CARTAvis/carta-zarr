/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// What a cache pool decoded is freed when the last holder of the pool lets go -- not when the image
// read through it does. A pool is made for one task, a moment or a cube histogram, and an image may
// stay open through any number of them; whatever outlived each would add up.
//
// Measured as the bytes this process has live, counted by replacing the global operator new. Every
// allocation is counted, on whatever thread made it, at the size the allocator actually gave it.

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <vector>

#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include "support/check.h"

#include <carta-zarr/carta_zarr.h>

namespace {

std::atomic<std::int64_t> live{0};

std::size_t SizeOf(void* memory) {
#if defined(__APPLE__)
    return malloc_size(memory);
#else
    return malloc_usable_size(memory);
#endif
}

void* Allocate(std::size_t size) {
    void* memory = std::malloc(size == 0 ? 1 : size);
    if (memory == nullptr) {
        throw std::bad_alloc();
    }
    live += static_cast<std::int64_t>(SizeOf(memory));
    return memory;
}

void* AllocateAligned(std::size_t size, std::align_val_t alignment) {
    void* memory = nullptr;
    const auto align = std::max(static_cast<std::size_t>(alignment), sizeof(void*));
    if (posix_memalign(&memory, align, size == 0 ? 1 : size) != 0) {
        throw std::bad_alloc();
    }
    live += static_cast<std::int64_t>(SizeOf(memory));
    return memory;
}

void Free(void* memory) noexcept {
    if (memory != nullptr) {
        live -= static_cast<std::int64_t>(SizeOf(memory));
        std::free(memory);
    }
}

}  // namespace

void* operator new(std::size_t size) {
    return Allocate(size);
}
void* operator new[](std::size_t size) {
    return Allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
    return AllocateAligned(size, alignment);
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
    return AllocateAligned(size, alignment);
}
void operator delete(void* memory) noexcept {
    Free(memory);
}
void operator delete[](void* memory) noexcept {
    Free(memory);
}
void operator delete(void* memory, std::size_t) noexcept {
    Free(memory);
}
void operator delete[](void* memory, std::size_t) noexcept {
    Free(memory);
}
void operator delete(void* memory, std::align_val_t) noexcept {
    Free(memory);
}
void operator delete[](void* memory, std::align_val_t) noexcept {
    Free(memory);
}
void operator delete(void* memory, std::size_t, std::align_val_t) noexcept {
    Free(memory);
}
void operator delete[](void* memory, std::size_t, std::align_val_t) noexcept {
    Free(memory);
}

namespace {

using carta::zarr::testing::Require;

// The wide fixture is 512 x 520 x 4 x 2 float32: 8.5 MB decoded, read whole through each pool.
void TestWhatAPoolDecodedGoesWithThePool() {
    auto context = carta::zarr::Context::Create();
    Require(static_cast<bool>(context), "Context::Create failed");
    const auto dataset = carta::zarr::Dataset::Open(*context, CARTA_ZARR_PIXEL_FIXTURE_WIDE);
    Require(static_cast<bool>(dataset), "Dataset::Open failed on the wide fixture");
    const auto image = dataset->OpenImage("SKY");
    Require(static_cast<bool>(image), "SKY could not be opened");

    carta::zarr::ReadRequest whole;
    std::size_t elements = 1;
    for (const auto& axis : image->descriptor().axes) {
        whole.axes.push_back({0, axis.length, 1});
        elements *= axis.length;
    }
    std::vector<float> pixels(elements);
    const auto decoded = static_cast<std::int64_t>(elements * sizeof(float));

    // One read first, so that what the store keeps for as long as it is open -- metadata, the handle
    // read through the session's pool -- is in the baseline rather than counted against the pools.
    Require(static_cast<bool>(image->Read(whole, {pixels.data(), pixels.size()})), "the first read failed");
    const std::int64_t before = live;
    constexpr int kPools = 5;
    {
        std::vector<carta::zarr::CachePool> pools;
        for (int i = 0; i < kPools; ++i) {
            auto pool = context->NewCachePool(std::size_t{64} << 20);
            Require(static_cast<bool>(pool), "NewCachePool failed");
            carta::zarr::ReadOptions options;
            options.control.cache_pool = *pool;
            const auto read = image->Read(whole, {pixels.data(), pixels.size()}, options);
            Require(static_cast<bool>(read), "a read through a pool failed");
            pools.push_back(*pool);
        }
        Require(live - before > decoded * (kPools - 1),
                "the pools did not keep what they decoded, so this test shows nothing about letting go of it: " +
                    std::to_string(live - before) + " bytes");
    }
    const std::int64_t after = live;
    Require(after - before < decoded, "letting go of every pool kept " + std::to_string(after - before) +
                                          " bytes, against " + std::to_string(decoded) + " decoded by each");
}

}  // namespace

int main() {
    try {
        TestWhatAPoolDecodedGoesWithThePool();
    } catch (const std::exception& error) {
        std::cerr << "pool release test failed: " << error.what() << "\n";
        return 1;
    }
    std::cout << "pool release tests passed\n";
    return 0;
}
