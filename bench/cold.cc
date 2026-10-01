/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "cold.h"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <system_error>

namespace carta::zarr::bench {

namespace {

constexpr const char* kDropCaches = "/proc/sys/vm/drop_caches";

bool HaveFadvise() {
#if defined(__linux__)
    return true;
#else
    return false;
#endif
}

bool CanDropCaches() {
    return geteuid() == 0 && access(kDropCaches, W_OK) == 0;
}

std::string RunCommand(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status == -1) {
        return "could not run the drop-cache command: " + std::string(std::strerror(errno));
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        return "the drop-cache command failed with status " + std::to_string(status);
    }
    return "";
}

std::string WriteDropCaches() {
    sync();
    std::ofstream file(kDropCaches);
    file << "3\n";
    file.flush();
    if (!file) {
        return std::string("could not write ") + kDropCaches;
    }
    return "";
}

// Every regular file under the dataset, told it will not be needed. Dirty pages are not dropped by
// this, which is why the dataset is synced first: the generator may have only just written it.
std::string Fadvise(const std::string& dataset) {
#if defined(__linux__)
    std::error_code walk;
    std::size_t failed = 0;
    std::string first_failure;
    for (auto entry = std::filesystem::recursive_directory_iterator(dataset, walk);
         !walk && entry != std::filesystem::recursive_directory_iterator(); entry.increment(walk)) {
        if (!entry->is_regular_file()) {
            continue;
        }
        const int fd = open(entry->path().c_str(), O_RDONLY | O_CLOEXEC);
        int error = fd < 0 ? errno : 0;
        if (fd >= 0) {
            fdatasync(fd);
            error = posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
            close(fd);
        }
        if (error != 0 && failed++ == 0) {
            first_failure = entry->path().string() + ": " + std::strerror(error);
        }
    }
    if (walk) {
        return "could not walk " + dataset + ": " + walk.message();
    }
    if (failed > 0) {
        return "fadvise failed on " + std::to_string(failed) + " files, the first " + first_failure;
    }
    return "";
#else
    (void)dataset;
    return "posix_fadvise is not available on this system";
#endif
}

}  // namespace

ColdChoice ChooseColdMethod(std::optional<ColdMethod> asked, const std::string& command) {
    if (!asked) {
        if (!command.empty()) {
            return {ColdMethod::command, ""};
        }
        if (CanDropCaches()) {
            return {ColdMethod::drop_caches, ""};
        }
        if (HaveFadvise()) {
            return {ColdMethod::fadvise, ""};
        }
        return {ColdMethod::off, ""};
    }
    switch (*asked) {
        case ColdMethod::drop_caches:
            if (!CanDropCaches()) {
                return {*asked, std::string("--cold drop-caches needs root and a writable ") + kDropCaches};
            }
            break;
        case ColdMethod::fadvise:
            if (!HaveFadvise()) {
                return {*asked, "--cold fadvise needs posix_fadvise, which this system does not have"};
            }
            break;
        case ColdMethod::command:
        case ColdMethod::off:
            break;
    }
    return {*asked, ""};
}

std::string DropCaches(ColdMethod method, const std::string& command, const std::string& dataset) {
    switch (method) {
        case ColdMethod::command:
            return RunCommand(command);
        case ColdMethod::drop_caches:
            return WriteDropCaches();
        case ColdMethod::fadvise:
            return Fadvise(dataset);
        case ColdMethod::off:
            return "";
    }
    return "";
}

}  // namespace carta::zarr::bench
