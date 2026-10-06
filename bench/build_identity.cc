/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#include "build_identity.h"

#include <carta-zarr/carta_zarr.h>
#include <dlfcn.h>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#include <array>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <system_error>

namespace carta::zarr::bench {
namespace {

constexpr std::uint64_t kFnvOffset = 0xCBF29CE484222325ULL;
constexpr std::uint64_t kFnvPrime = 0x100000001B3ULL;

void Mix(std::uint64_t& hash, const char* bytes, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        hash ^= static_cast<unsigned char>(bytes[i]);
        hash *= kFnvPrime;
    }
}

std::string ExecutablePath() {
#if defined(__APPLE__)
    std::uint32_t size = 0;
    (void)_NSGetExecutablePath(nullptr, &size);
    std::string path(size, '\0');
    if (_NSGetExecutablePath(path.data(), &size) != 0) {
        return {};
    }
    path.resize(path.find('\0') == std::string::npos ? path.size() : path.find('\0'));
    return path;
#else
    std::error_code ignored;
    return std::filesystem::read_symlink("/proc/self/exe", ignored).string();
#endif
}

// The carta-zarr this process loaded: where the loader found the library ProbeSchema lives in --
// exported and defined there, unlike the header's inline functions, whose copies are the caller's.
// The executable itself when it was linked in statically.
std::string LibraryPath() {
    Dl_info info{};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): dladdr takes an address.
    if (dladdr(reinterpret_cast<const void*>(&carta::zarr::ProbeSchema), &info) == 0 || info.dli_fname == nullptr) {
        return {};
    }
    return info.dli_fname;
}

void MixFile(std::uint64_t& hash, const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (path.empty() || !file) {
        const std::string missing = "unreadable:" + path;
        Mix(hash, missing.data(), missing.size() + 1);
        return;
    }
    std::array<char, 1 << 16> buffer{};
    while (file) {
        file.read(buffer.data(), buffer.size());
        Mix(hash, buffer.data(), static_cast<std::size_t>(file.gcount()));
    }
    // Where one file ends, so that two whose bytes run together differently are told apart.
    Mix(hash, "\0", 1);
}

}  // namespace

std::string BuildIdentity() {
    static const std::string identity = [] {
        std::uint64_t hash = kFnvOffset;
        MixFile(hash, ExecutablePath());
        MixFile(hash, LibraryPath());
        std::array<char, 17> text{};
        std::snprintf(text.data(), text.size(), "%016" PRIx64, hash);
        return std::string(text.data());
    }();
    return identity;
}

}  // namespace carta::zarr::bench
