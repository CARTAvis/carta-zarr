/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "carta-zarr/carta_zarr.h"

#include "chunk_blocks.h"
#include "read/pieces.h"
#include "reduce/plane_histogram.h"
#include "reduce/spectral_reduce.h"
#include "reduce/store_slab_source.h"
#include "schema/profile.h"
#include "store.h"
#include "work_pool.h"
#include "zarr/array_metadata.h"
#include "zarr/pixel_reader.h"
#include "zarr/store_context.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <cstdint>
#include <limits>
#include <mutex>
#include <optional>
#include <vector>
#include <unordered_map>
#include <utility>

namespace carta::zarr {
namespace {

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

// Every public entry point reports its failures as a Result, and a consumer that checks one should
// never have to catch as well. Underneath, though, metadata comes from a file and buffers are sized
// from what it says: nlohmann throws on a value that is not the type it is read as, and the standard
// library throws on a length it cannot allocate. This is where that becomes an Error.
template <typename Function>
auto Guarded(ErrorCode code, std::string node, Function&& function) -> decltype(function()) {
    try {
        return function();
    } catch (const std::exception& error) {
        return MakeError(code, error.what(), std::move(node));
    }
}

// A probe that rejects a store has already worked out why, often down to the attribute, and a
// consumer shows whatever comes back here to whoever picked the file. Reporting "not a supported
// dataset" instead throws that away and names a schema profile the store may have nothing to do
// with. The fallback is for the case the probe genuinely had nothing to say.
std::string RejectionMessage(const std::vector<Diagnostic>& diagnostics, std::string fallback) {
    return diagnostics.empty() ? std::move(fallback) : diagnostics.front().message;
}

bool TryComputeDirectorySize(std::string_view location, std::chrono::milliseconds timeout, std::uint64_t& size) {
    const std::string location_string(location);
    std::filesystem::path root_path;
    if (location_string.rfind("file://", 0) == 0) {
        root_path = std::filesystem::path(location_string.substr(7));
    } else if (location_string.find("://") != std::string::npos) {
        return false;
    } else {
        root_path = std::filesystem::path(location_string);
    }

    const auto deadline = std::chrono::steady_clock::now() + timeout;
    std::uint64_t directory_size = 0;
    std::error_code error;
    std::filesystem::recursive_directory_iterator iterator(
        root_path, std::filesystem::directory_options::skip_permission_denied, error);
    if (error) {
        return false;
    }
    const std::filesystem::recursive_directory_iterator end;
    while (iterator != end) {
        if (std::chrono::steady_clock::now() >= deadline) {
            return false;
        }

        std::error_code entry_error;
        if (iterator->is_regular_file(entry_error)) {
            if (entry_error) {
                return false;
            }
            const auto file_size = iterator->file_size(entry_error);
            if (entry_error || file_size > std::numeric_limits<std::uint64_t>::max() - directory_size) {
                return false;
            }
            directory_size += file_size;
        } else if (entry_error) {
            return false;
        }

        iterator.increment(error);
        if (error) {
            return false;
        }
    }

    size = directory_size;
    return true;
}

}  // namespace

class Context::Impl {
public:
    Impl(OpenOptions options, internal::StoreContextPtr store_context)
        : options(options),
          store_context(std::move(store_context)),
          // decode_threads is the consumer's statement of how much of this machine the library may
          // use, so it sizes both pools rather than only TensorStore's. The two are busy at
          // different moments -- a slab is read and then visited -- so sizing each at the whole
          // budget does not double the demand. Zero means one worker per hardware thread, which is
          // what TensorStore's own default does with the same number.
          workers(std::make_shared<internal::WorkPool>(options.decode_threads)) {}

    OpenOptions options;
    // Shared by every dataset and image opened through this context, so that its cache and
    // concurrency limits apply to all reads rather than being rebuilt per read.
    internal::StoreContextPtr store_context;
    // The per-pixel work of a reduction. See WorkPool.
    std::shared_ptr<internal::WorkPool> workers;
};

Context::Context(std::shared_ptr<Impl> impl) : _impl(std::move(impl)) {}
Context::~Context() = default;

Result<Context> Context::Create(const OpenOptions& options) {
    // Guarded like every other public entry point, and for a reason the others do not have: building
    // an Impl starts the worker threads, and a system that refuses one throws std::system_error.
    // Without this, the one call a consumer makes before it can do anything else is also the only
    // one that could throw at it. io_error is the nearest existing code for "the machine would not
    // give us what we asked for"; the message says which resource it was.
    return Guarded(ErrorCode::io_error, {}, [&]() -> Result<Context> {
        auto store_context = internal::MakeStoreContext(options);
        if (!store_context) {
            return store_context.error();
        }
        return Context{std::make_shared<Impl>(options, std::move(store_context.value()))};
    });
}

class Image::Impl {
public:
    Impl(std::shared_ptr<Context::Impl> context, std::string location, internal::SchemaProfile profile,
         std::shared_ptr<internal::Store> store, ImageDescriptor descriptor, ChunkGeometry geometry)
        : context(std::move(context)),
          location(std::move(location)),
          profile(profile),
          store(std::move(store)),
          descriptor(std::move(descriptor)),
          geometry(std::move(geometry)) {}

    std::shared_ptr<Context::Impl> context;
    std::string location;
    // The profile that described this image, rather than the name of one to look up again. It is a
    // pointer into a table that outlives every store, so holding it costs nothing and removes an
    // error path that could only fire if a descriptor named a profile the library does not have.
    internal::SchemaProfile profile;
    std::shared_ptr<internal::Store> store;
    ImageDescriptor descriptor;
    ChunkGeometry geometry;
};

namespace {

ChunkGeometry BuildChunkGeometry(const ImageDescriptor& descriptor, const StorageLayout& layout) {
    ChunkGeometry geometry;
    geometry.sharded = layout.sharded;
    geometry.compressor = layout.compressor;

    const auto rank = descriptor.axes.size();
    geometry.chunk_shape.resize(rank);
    geometry.shard_shape.resize(rank);
    geometry.grid_shape.resize(rank);
    for (std::size_t logical = 0; logical < rank; ++logical) {
        const auto& axis = descriptor.axes.at(logical);
        const auto stored = axis.storage_index;
        const auto chunk =
            stored < layout.chunk_shape.size() ? layout.chunk_shape.at(stored) : axis.length;
        const auto shard =
            stored < layout.shard_shape.size() ? layout.shard_shape.at(stored) : chunk;
        geometry.chunk_shape.at(logical) = chunk;
        geometry.shard_shape.at(logical) = shard == 0 ? chunk : shard;
        geometry.grid_shape.at(logical) = chunk == 0 ? 0 : (axis.length + chunk - 1) / chunk;
        if (stored != logical) {
            geometry.transpose_required = true;
        }
    }

    // The last stored dimension varies fastest, so of the two spatial axes the one with the larger
    // storage index is the one a plane is contiguous along.
    std::size_t x_stored = 0;
    std::size_t y_stored = 0;
    for (const auto& axis : descriptor.axes) {
        if (axis.role == AxisRole::spatial_x) {
            x_stored = axis.storage_index;
        } else if (axis.role == AxisRole::spatial_y) {
            y_stored = axis.storage_index;
        }
    }
    geometry.fastest_spatial_axis = y_stored > x_stored ? AxisRole::spatial_y : AxisRole::spatial_x;
    return geometry;
}

}  // namespace

Image::Image(std::shared_ptr<Impl> impl) : _impl(std::move(impl)) {}
Image::~Image() = default;

const ImageDescriptor& Image::descriptor() const noexcept {
    static const ImageDescriptor empty_descriptor;
    return _impl ? _impl->descriptor : empty_descriptor;
}

const ChunkGeometry& Image::chunk_geometry() const noexcept {
    static const ChunkGeometry empty_geometry;
    return _impl ? _impl->geometry : empty_geometry;
}

Result<std::size_t> Image::Read(const ReadRequest& request, BufferView<float> destination) const {
    return Read(request, destination, ReadOptions{});
}

Result<std::size_t> Image::Read(const ReadRequest& request, BufferView<float> destination,
                                const ReadOptions& options) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<std::size_t> {
        if (!_impl) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        return internal::ReadInPieces(*_impl->store, _impl->descriptor, _impl->geometry, request, destination,
                                      options);
    });
}

Result<std::size_t> Image::ReadPixelMask(const ReadRequest& request, BufferView<std::uint8_t> destination) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<std::size_t> {
        if (!_impl) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        if (!_impl->descriptor.has_pixel_mask) {
            return MakeError(ErrorCode::not_found, "This image has no pixel mask", _impl->descriptor.id);
        }

        auto selection = internal::zarr::BuildSelection(_impl->descriptor, request);
        if (!selection) {
            return selection.error();
        }
        const auto elements = internal::zarr::SelectionElementCount(selection.value());
        if (elements == 0 || elements > destination.size) {
            return MakeError(ErrorCode::invalid_argument, "Destination buffer is too small for the request",
                             _impl->descriptor.id);
        }

        auto read = _impl->store->ReadPixelMaskBytes(_impl->descriptor.pixel_mask_id, selection.value(),
                                                     destination.data, static_cast<std::size_t>(elements),
                                                     ReadOptions{});
        if (!read) {
            return read.error();
        }
        return static_cast<std::size_t>(elements);
    });
}

Result<void> Image::ReduceSpectral(const SpectralReduceRequest& request, const SpectralSink& sink) const {
    return ReduceSpectral(request, sink, ReadOptions{});
}

Result<void> Image::ReduceSpectral(const SpectralReduceRequest& request, const SpectralSink& sink,
                                   const ReadOptions& options) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<void> {
        if (!_impl || !_impl->store) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        const internal::StoreSlabSource source(*_impl->store, _impl->descriptor);
        return internal::ReduceSpectral(source, _impl->descriptor, _impl->geometry, request, sink, options,
                                        *_impl->context->workers);
    });
}

Result<void> Image::ComputeHistogram(const HistogramRequest& request, const HistogramSink& sink) const {
    return ComputeHistogram(request, sink, ReadOptions{});
}

Result<void> Image::ComputeHistogram(const HistogramRequest& request, const HistogramSink& sink,
                                     const ReadOptions& options) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<void> {
        if (!_impl || !_impl->store) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        const internal::StoreSlabSource source(*_impl->store, _impl->descriptor);
        return internal::ComputeHistogram(source, _impl->descriptor, _impl->geometry, request, sink, options,
                                          *_impl->context->workers);
    });
}

Result<CubeHistogramResult> Image::ComputeCubeHistogram(const CubeHistogramRequest& request) const {
    return ComputeCubeHistogram(request, ReadOptions{});
}

Result<CubeHistogramResult> Image::ComputeCubeHistogram(const CubeHistogramRequest& request,
                                                        const ReadOptions& options) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<CubeHistogramResult> {
        if (!_impl || !_impl->store) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        const internal::StoreSlabSource source(*_impl->store, _impl->descriptor);
        return internal::ComputeCubeHistogram(source, _impl->descriptor, _impl->geometry, request, options,
                                              *_impl->context->workers);
    });
}

Result<std::vector<Beam>> Image::ReadBeams() const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::invalid_metadata, node, [&]() -> Result<std::vector<Beam>> {
        if (!_impl) {
            return MakeError(ErrorCode::invalid_argument, "Image handle is empty");
        }
        if (!_impl->store) {
            return MakeError(ErrorCode::invalid_argument, "Image store is unavailable");
        }
        return _impl->profile.ReadBeams(*_impl->store, _impl->descriptor.id);
    });
}

class Dataset::Impl {
public:
    Impl(std::shared_ptr<Context::Impl> context, std::string location, DatasetDescriptor descriptor,
         internal::SchemaProfile profile, internal::Store store)
        : context(std::move(context)),
          location(std::move(location)),
          descriptor(std::move(descriptor)),
          profile(profile),
          store(std::make_shared<internal::Store>(std::move(store))) {}

    std::shared_ptr<Context::Impl> context;
    std::string location;
    DatasetDescriptor descriptor;
    // The profile the probe matched. descriptor.schema_id names it for the consumer; this is the
    // one the library asks, resolved where the match happened rather than at each use.
    internal::SchemaProfile profile;
    std::shared_ptr<internal::Store> store;
    mutable std::mutex mutex;
    mutable std::unordered_map<std::string, ImageDescriptor> image_descriptors;
};

Dataset::Dataset(std::shared_ptr<Impl> impl) : _impl(std::move(impl)) {}
Dataset::~Dataset() = default;

Result<Dataset> Dataset::Open(const Context& context, std::string_view location) {
    try {
        if (!context._impl) {
            return MakeError(ErrorCode::invalid_argument, "Context handle is empty");
        }

        // TensorStore resources are shared by Context, while array handles are scoped to this
        // Dataset and the Images that retain its Store.
        auto store_context = context._impl->store_context->CloneForStore();
        auto store_result = internal::OpenStore(location, std::move(store_context));
        if (!store_result) {
            return store_result.error();
        }
        auto probe_result = internal::ProbeStore(store_result.value());
        if (!probe_result) {
            return probe_result.error();
        }
        const auto& probe = probe_result.value();
        if (probe.kind != ProbeKind::supported_dataset) {
            const ErrorCode code =
                probe.kind == ProbeKind::invalid_dataset ? ErrorCode::invalid_metadata : ErrorCode::unsupported_schema;
            return MakeError(
                code, RejectionMessage(probe.diagnostics, "No built-in schema profile matched the Zarr store"),
                std::string(location));
        }

        if (probe.images.empty()) {
            return MakeError(ErrorCode::invalid_metadata, "Supported schema has no image variables",
                             std::string(location));
        }
        auto profile = internal::SchemaProfile::For(probe.schema_id);
        if (!profile) {
            return profile.error();
        }
        // The descriptor is the probe's answer without the question it was answering, so it is taken
        // rather than rebuilt field by field -- and taken by move, since the probe result dies here.
        DatasetDescriptor descriptor = std::move(static_cast<DatasetDescriptor&>(probe_result.value()));
        return Dataset{std::make_shared<Impl>(context._impl, std::string(location), std::move(descriptor),
                                              profile.value(), std::move(store_result.value()))};
    } catch (const std::exception& error) {
        return MakeError(ErrorCode::invalid_metadata, error.what(), std::string(location));
    }
}

const DatasetDescriptor& Dataset::descriptor() const noexcept {
    static const DatasetDescriptor empty_descriptor;
    return _impl ? _impl->descriptor : empty_descriptor;
}

Result<DatasetSize> Dataset::Size(std::chrono::milliseconds directory_size_timeout) const {
    const std::string node = _impl ? _impl->location : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<DatasetSize> {
        if (!_impl) {
            return MakeError(ErrorCode::invalid_argument, "Dataset handle is empty");
        }

        std::uint64_t physical_size = 0;
        if (TryComputeDirectorySize(_impl->location, directory_size_timeout, physical_size)) {
            return DatasetSize{physical_size, false};
        }

        auto logical_size = _impl->store->ComputeTotalArraySizeBytes();
        if (!logical_size) {
            return logical_size.error();
        }
        return DatasetSize{logical_size.value(), true};
    });
}

Result<Image> Dataset::OpenImage(std::string_view image_id) const {
    const std::string node(image_id);
    return Guarded(ErrorCode::invalid_metadata, node, [&]() -> Result<Image> {
        if (!_impl) {
            return MakeError(ErrorCode::invalid_argument, "Dataset handle is empty");
        }
        std::scoped_lock const lock(_impl->mutex);
        const std::string image_name(image_id);
        const auto entry = std::find_if(_impl->descriptor.images.begin(), _impl->descriptor.images.end(),
                                        [&](const ImageEntry& image) { return image.id == image_name; });
        if (entry == _impl->descriptor.images.end()) {
            return MakeError(ErrorCode::not_found, "Image variable was not found", image_name);
        }
        if (!entry->readable) {
            const auto message = entry->diagnostics.empty() ? "Image variable is not openable by this profile"
                                                            : entry->diagnostics.front().message;
            return MakeError(ErrorCode::unsupported_data_type, message, image_name);
        }
        const auto make_image = [&](const ImageDescriptor& descriptor) {
            // The descriptor already carries the stored layout; the geometry is that layout permuted
            // into logical order, so it is derived here rather than read again.
            const StorageLayout layout = descriptor.storage ? *descriptor.storage : StorageLayout{};
            return Image{std::make_shared<Image::Impl>(_impl->context, _impl->location, _impl->profile, _impl->store,
                                                       descriptor, BuildChunkGeometry(descriptor, layout))};
        };
        const auto cached = _impl->image_descriptors.find(image_name);
        if (cached != _impl->image_descriptors.end()) {
            return make_image(cached->second);
        }
        auto image_descriptor = _impl->profile.DescribeVerified(*_impl->store, image_id);
        if (!image_descriptor) {
            return image_descriptor.error();
        }
        auto [inserted, _] = _impl->image_descriptors.emplace(image_name, std::move(image_descriptor.value()));
        return make_image(inserted->second);
    });
}

ProbeResult Probe(std::string_view location, const ProbeOptions&) {
    try {
        auto store_result = internal::OpenStore(location);
        if (!store_result) {
            ProbeResult result;
            result.kind = ProbeKind::not_zarr;
            if (store_result.error().code == ErrorCode::invalid_metadata ||
                store_result.error().code == ErrorCode::io_error) {
                result.kind = ProbeKind::invalid_dataset;
            }
            result.diagnostics.push_back(Diagnostic{internal::zarr::ErrorCodeName(store_result.error().code),
                                                    store_result.error().message, store_result.error().node_path});
            return result;
        }
        auto probe_result = internal::ProbeStore(store_result.value());
        if (!probe_result) {
            ProbeResult result;
            result.kind = ProbeKind::invalid_dataset;
            result.diagnostics.push_back(Diagnostic{internal::zarr::ErrorCodeName(probe_result.error().code),
                                                    probe_result.error().message, probe_result.error().node_path});
            return result;
        }
        return probe_result.value();
    } catch (const std::exception& error) {
        ProbeResult result;
        result.kind = ProbeKind::invalid_dataset;
        result.diagnostics.push_back(Diagnostic{"exception", error.what(), std::string(location)});
        return result;
    }
}

Result<SchemaProbeResult> ProbeSchema(std::string_view location, std::string_view schema_id) {
    try {
        // Resolve the profile before touching the store, so an unknown schema reports itself rather
        // than whatever happens to be wrong with the path.
        auto profile = internal::SchemaProfile::For(schema_id);
        if (!profile) {
            return profile.error();
        }
        auto store_result = internal::OpenStore(location);
        if (!store_result) {
            return store_result.error();
        }
        return profile.value().Probe(store_result.value());
    } catch (const std::exception& error) {
        return MakeError(ErrorCode::invalid_metadata, error.what(), std::string(location));
    }
}

Result<bool> IsXradioImage(std::string_view location) {
    auto result = ProbeSchema(location, kXradioImageSchema);
    if (!result) {
        return result.error();
    }
    if (result.value().kind == SchemaMatchKind::match) {
        return true;
    }
    if (result.value().kind == SchemaMatchKind::invalid) {
        return MakeError(ErrorCode::invalid_metadata,
                         RejectionMessage(result.value().diagnostics, "The requested schema did not match"),
                         std::string(location));
    }
    return false;
}

}  // namespace carta::zarr
