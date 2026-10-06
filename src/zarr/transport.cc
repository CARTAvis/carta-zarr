/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "transport.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace carta::zarr::internal {
namespace {

// A walk of the store that follows directory links, and where it must not. What a link points at is
// in the store as much as what sits there -- a dataset's arrays are as often as not links to where
// their bytes are, and its consolidated copy names what is under a linked group -- so a walk that
// stopped at links lost those nodes and that size. A link to a directory the walk is already inside
// is a cycle: it is neither listed, which would name that directory again, nor followed.
//
// Holds the real path of each directory the walk is inside, the root first. An ordinary directory's
// is its parent's and its name; only a link's is asked of the filesystem.
//
// A directory that refuses to be read is either passed over or the walk's error, as the walk says. A
// listing passes over it: no node there could be read, and a later read says why. A size cannot: one
// that leaves a directory out is not the store's size, however it is labelled.
enum class Unreadable { skipped, failed };

class LinkedWalk {
public:
    // `root` is resolved already: the transport's root always is.
    LinkedWalk(const std::filesystem::path& root, Unreadable unreadable)
        : iterator(root,
                   unreadable == Unreadable::skipped ? std::filesystem::directory_options::skip_permission_denied |
                                                           std::filesystem::directory_options::follow_directory_symlink
                                                     : std::filesystem::directory_options::follow_directory_symlink,
                   error),
          _ancestors{root} {}

    // Whether the directory the walk is at is one to take: false, with recursion into it turned off,
    // for a link back to a directory the walk is inside or one that cannot be resolved.
    bool Enter() {
        const auto& path = iterator->path();
        _ancestors.resize(static_cast<std::size_t>(iterator.depth()) + 1);
        std::filesystem::path real = _ancestors.back() / path.filename();
        std::error_code entry_error;
        if (iterator->is_symlink(entry_error)) {
            real = std::filesystem::canonical(path, entry_error);
            if (entry_error || std::find(_ancestors.begin(), _ancestors.end(), real) != _ancestors.end()) {
                iterator.disable_recursion_pending();
                return false;
            }
        }
        _ancestors.push_back(std::move(real));
        return true;
    }

    std::error_code error;
    std::filesystem::recursive_directory_iterator iterator;

private:
    std::vector<std::filesystem::path> _ancestors;
};

Result<std::filesystem::path> NormalizeLocation(std::string_view location) {
    if (location.empty()) {
        return Error{ErrorCode::invalid_argument, "Zarr location must not be empty"};
    }

    std::string const location_string(location);
    std::filesystem::path path;
    if (location_string.rfind("file://", 0) == 0) {
        path = std::filesystem::path(location_string.substr(7));
    } else if (location_string.find("://") != std::string::npos) {
        return Error{ErrorCode::unsupported_transport, "Only local filesystem and file:// Zarr stores are supported"};
    } else {
        path = std::filesystem::path(location_string);
    }

    if (path.empty()) {
        return Error{ErrorCode::invalid_argument, "Zarr location must not be empty"};
    }
    // Resolved against the working directory here and nowhere else. A store keeps its root for as
    // long as the images opened from it live, so a root left relative makes every later read depend
    // on the process still being where it was when the consumer opened the file -- and it is the
    // pixel reads that would fail, long after the open that looked fine. Doing it once at the root
    // also keeps it off the per-read path, where it would be a handful of lstat calls per cursor step.
    std::error_code error;
    auto resolved = std::filesystem::weakly_canonical(std::filesystem::absolute(path, error), error);
    if (error) {
        return Error{ErrorCode::io_error, "Unable to resolve Zarr location: " + error.message(), path.string()};
    }
    return resolved;
}

class FilesystemTransport final : public Transport {
public:
    explicit FilesystemTransport(std::filesystem::path root) : _root(std::move(root)) {}

    Result<std::string> ReadNodeBytes(std::string_view node) const override {
        const std::filesystem::path metadata_path = node.empty() ? _root / "zarr.json" : _root / node / "zarr.json";
        std::error_code error;
        if (!std::filesystem::exists(metadata_path, error)) {
            if (error) {
                return Error{ErrorCode::io_error, "Unable to inspect Zarr node metadata: " + error.message(),
                             metadata_path.string()};
            }
            return Error{ErrorCode::not_found, "Zarr node is missing zarr.json", metadata_path.string()};
        }

        // Opened at the end so that the document can be sized, then read in one go. Streaming
        // rdbuf() into an ostringstream instead moves it a character at a time through two stream
        // buffers and a growing string, which measured as a quarter of the cost of opening a
        // consolidated store -- more than the read it was there to perform.
        std::ifstream input(metadata_path, std::ios::binary | std::ios::ate);
        if (!input.is_open()) {
            return Error{ErrorCode::io_error, "Unable to read Zarr metadata", metadata_path.string()};
        }
        const auto size = input.tellg();
        if (size < 0) {
            return Error{ErrorCode::io_error, "Unable to read Zarr metadata", metadata_path.string()};
        }
        input.seekg(0);

        // An empty file is not an I/O failure: it reaches the JSON parser above the seam and is
        // reported as invalid metadata, which is what it is. Asking for no bytes would set failbit
        // and say nothing, so it is not asked for.
        std::string bytes(static_cast<std::size_t>(size), '\0');
        if (!bytes.empty()) {
            input.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
            if (input.bad()) {
                return Error{ErrorCode::io_error, "Unable to read Zarr metadata", metadata_path.string()};
            }
            // A file that shrank between being sized and being read hands back what was there.
            bytes.resize(static_cast<std::size_t>(input.gcount()));
        }
        return bytes;
    }

    Result<std::vector<std::string>> ListNodes() const override {
        std::vector<std::string> nodes;
        LinkedWalk walk(_root, Unreadable::skipped);
        auto& iterator = walk.iterator;
        std::error_code& error = walk.error;
        const std::filesystem::recursive_directory_iterator end;
        for (; iterator != end; iterator.increment(error)) {
            if (error) {
                return Error{ErrorCode::io_error, "Unable to enumerate Zarr metadata: " + error.message(),
                             _root.string()};
            }

            // Only a directory is ever a node, and only the directory is listed. A group's own
            // zarr.json is a file inside a directory that has already been named, so counting
            // files as well named every group twice -- once as itself and once as its metadata's
            // parent. Store sorts and de-duplicates the listing, which is why that cost a second
            // walk of the tree rather than a wrong answer, and why nothing noticed.
            std::error_code entry_error;
            if (!iterator->is_directory(entry_error) || entry_error) {
                continue;
            }

            if (!walk.Enter()) {
                continue;
            }
            const auto path = iterator->path();

            const auto metadata_path = path / "zarr.json";
            // A directory carrying no readable zarr.json is not a node. Whether the answer was no
            // or the question could not be asked makes no difference here: neither is a node, and
            // a directory that could not be inspected is left for a later read to diagnose.
            std::error_code metadata_error;
            if (!std::filesystem::is_regular_file(metadata_path, metadata_error)) {
                continue;
            }

            // Named by where it sits in the store, which the walk spelt from the root: relative()
            // resolves links first, so an array linked in from elsewhere came out as
            // ../elsewhere/SKY, a name the store refuses.
            const auto relative_path = path.lexically_relative(_root);
            if (relative_path.empty()) {
                continue;
            }
            nodes.push_back(relative_path.generic_string());

            // Array chunks are descendants of the array node, but are not Zarr nodes. Inspect only
            // the small node header to avoid walking millions of chunk files. Malformed metadata is
            // left for Store::ReadNodeMetadata to diagnose; in that case traversal remains conservative.
            std::ifstream input(metadata_path);
            nlohmann::json metadata;
            if (input.is_open()) {
                try {
                    input >> metadata;
                    if (metadata.is_object() && metadata.value("node_type", "") == "array") {
                        iterator.disable_recursion_pending();
                    }
                } catch (const std::exception&) {
                    // The metadata parser above the transport seam reports the definitive error.
                }
            }
        }
        if (error) {
            return Error{ErrorCode::io_error, "Unable to enumerate Zarr metadata: " + error.message(), _root.string()};
        }
        return nodes;
    }

    Result<std::uint64_t> StoredSizeBytes(std::chrono::steady_clock::time_point deadline) const override {
        std::uint64_t total = 0;
        LinkedWalk walk(_root, Unreadable::failed);
        auto& iterator = walk.iterator;
        std::error_code& error = walk.error;
        if (error) {
            return Error{ErrorCode::io_error, "Unable to enumerate the Zarr store: " + error.message(), _root.string()};
        }

        const std::filesystem::recursive_directory_iterator end;
        while (iterator != end) {
            // Checked per entry rather than per directory: a store is mostly chunk files, so the
            // entries are where the time goes and a deadline checked a level up would overrun a
            // wide array by however long that array takes.
            if (std::chrono::steady_clock::now() >= deadline) {
                return Error{ErrorCode::cancelled, "The Zarr store size deadline expired", _root.string()};
            }

            std::error_code entry_error;
            if (iterator->is_directory(entry_error) && !entry_error) {
                walk.Enter();
            } else if (iterator->is_regular_file(entry_error)) {
                const auto file_size = iterator->file_size(entry_error);
                if (entry_error) {
                    return Error{ErrorCode::io_error, "Unable to size a Zarr store file: " + entry_error.message(),
                                 iterator->path().string()};
                }
                if (file_size > std::numeric_limits<std::uint64_t>::max() - total) {
                    return Error{ErrorCode::io_error, "The Zarr store's size overflows uint64_t", _root.string()};
                }
                total += file_size;
            } else if (entry_error) {
                return Error{ErrorCode::io_error, "Unable to inspect a Zarr store entry: " + entry_error.message(),
                             iterator->path().string()};
            }

            iterator.increment(error);
            if (error) {
                return Error{ErrorCode::io_error, "Unable to enumerate the Zarr store: " + error.message(),
                             _root.string()};
            }
        }
        return total;
    }

    // A backslash is refused here although NormalizeNodeName, the store's rule, accepts one: that
    // is TensorStore's opinion rather than this transport's. Its file kvstore reads a backslash as
    // a separator, so the array of a node named "SKY\2" would be looked for in "SKY/2". On POSIX the
    // name is a legal file name and the node's metadata reads, so such a variable is listed and
    // opens, and every read of it is refused here, saying why. No XRADIO writer names a variable
    // that way; a listing that knew would need the store to know what TensorStore can locate.
    Result<std::filesystem::path> ArrayDirectory(std::string_view node) const override {
        const std::filesystem::path relative(node);
        if (relative.empty() || relative.is_absolute() || relative.has_root_name() ||
            node.find('\\') != std::string_view::npos) {
            return Error{
                ErrorCode::invalid_argument,
                "Invalid Zarr array path " + std::string(node) +
                    (node.find('\\') != std::string_view::npos ? ": TensorStore reads a backslash as a path separator"
                                                               : ""),
                std::string(node)};
        }
        for (const auto& component : relative) {
            if (component == "." || component == "..") {
                return Error{ErrorCode::invalid_argument, "Invalid Zarr array path " + std::string(node),
                             std::string(node)};
            }
        }
        return _root / relative;
    }

private:
    std::filesystem::path _root;
};

}  // namespace

Result<TransportPtr> OpenFilesystemTransport(std::string_view location) {
    auto path_result = NormalizeLocation(location);
    if (!path_result) {
        return path_result.error();
    }
    const std::filesystem::path& path = path_result.value();

    std::error_code error;
    if (!std::filesystem::exists(path, error)) {
        if (error) {
            return Error{ErrorCode::io_error, "Unable to inspect Zarr location: " + error.message()};
        }
        return Error{ErrorCode::not_found, "Zarr location does not exist", path.string()};
    }
    if (error || !std::filesystem::is_directory(path, error)) {
        return Error{ErrorCode::not_zarr, "Zarr location is not a directory", path.string()};
    }

    return TransportPtr(std::make_shared<const FilesystemTransport>(path));
}

}  // namespace carta::zarr::internal
