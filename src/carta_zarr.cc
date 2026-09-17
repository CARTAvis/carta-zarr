/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "carta-zarr/carta_zarr.h"

#include "read/pieces.h"
#include "readable_image.h"
#include "reduce/plane_histogram.h"
#include "reduce/spectral_reduce.h"
#include "schema/chunk_geometry.h"
#include "schema/profile.h"
#include "store.h"
#include "store_pixel_source.h"
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

// Every public entry point reports its failures as a Result, and a consumer that checks one should
// never have to catch as well. Underneath, though, metadata comes from a file and buffers are sized
// from what it says: nlohmann throws on a value that is not the type it is read as, and the standard
// library throws on a length it cannot allocate. This is where that becomes a report.
//
// It is the only try in this file. Three entry points used to write their own beside the twelve
// that went through here, with fallback codes of their own, so whether one caught was a question you
// answered by reading to the end of it.
template <typename Function, typename OnThrow>
auto GuardedWith(Function&& function, OnThrow&& on_throw) -> decltype(function()) {
    try {
        return function();
    } catch (const std::exception& error) {
        return on_throw(error);
    }
}

template <typename Function>
auto Guarded(ErrorCode code, std::string node, Function&& function) -> decltype(function()) {
    return GuardedWith(std::forward<Function>(function), [&](const std::exception& error) {
        return Error{code, error.what(), std::move(node)};
    });
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
          geometry(std::move(geometry)),
          source(*this->store, this->descriptor) {}

    std::shared_ptr<Context::Impl> context;
    std::string location;
    // The profile that described this image, rather than the name of one to look up again. It is a
    // pointer into a table that outlives every store, so holding it costs nothing and removes an
    // error path that could only fire if a descriptor named a profile the library does not have.
    internal::SchemaProfile profile;
    std::shared_ptr<internal::Store> store;
    ImageDescriptor descriptor;
    ChunkGeometry geometry;
    // Built once, from members rather than from the constructor's arguments, and declared last so
    // that both of those are already initialised. A PixelSource cannot be copied or moved, which is
    // why it lives here rather than being made at each entry point.
    internal::StorePixelSource source;

    // What every read and reduction is against. Built per call, so that an image whose axes cannot
    // be mapped fails the operation that needed them rather than the open -- which is where that
    // failure has always reached the caller.
    Result<internal::ReadableImage> Readable() const {
        return internal::ReadableImage::Of(source, descriptor, geometry, *context->workers);
    }
};

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
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        auto image = _impl->Readable();
        if (!image) {
            return image.error();
        }
        return internal::ReadInPieces(image.value(), request, destination, options);
    });
}

Result<std::size_t> Image::ReadPixelMask(const ReadRequest& request, BufferView<std::uint8_t> destination) const {
    return ReadPixelMask(request, destination, ReadOptions{});
}

Result<std::size_t> Image::ReadPixelMask(const ReadRequest& request, BufferView<std::uint8_t> destination,
                                         const ReadOptions& options) const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<std::size_t> {
        if (!_impl) {
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        auto image = _impl->Readable();
        if (!image) {
            return image.error();
        }
        return internal::ReadPixelMask(image.value(), request, destination, options);
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
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        auto image = _impl->Readable();
        if (!image) {
            return image.error();
        }
        return internal::ReduceSpectral(image.value(), request, sink, options);
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
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        auto image = _impl->Readable();
        if (!image) {
            return image.error();
        }
        return internal::ComputeHistogram(image.value(), request, sink, options);
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
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        auto image = _impl->Readable();
        if (!image) {
            return image.error();
        }
        return internal::ComputeCubeHistogram(image.value(), request, options);
    });
}

Result<std::vector<Beam>> Image::ReadBeams() const {
    const std::string node = _impl ? _impl->descriptor.id : std::string{};
    return Guarded(ErrorCode::invalid_metadata, node, [&]() -> Result<std::vector<Beam>> {
        if (!_impl) {
            return Error{ErrorCode::invalid_argument, "Image handle is empty"};
        }
        if (!_impl->store) {
            return Error{ErrorCode::invalid_argument, "Image store is unavailable"};
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
    return Guarded(ErrorCode::invalid_metadata, std::string(location), [&]() -> Result<Dataset> {
        if (!context._impl) {
            return Error{ErrorCode::invalid_argument, "Context handle is empty"};
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
        if (auto openable = internal::RequireOpenableDataset(probe, location); !openable) {
            return openable.error();
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
    });
}

const DatasetDescriptor& Dataset::descriptor() const noexcept {
    static const DatasetDescriptor empty_descriptor;
    return _impl ? _impl->descriptor : empty_descriptor;
}

Result<DatasetSize> Dataset::Size(std::chrono::milliseconds directory_size_timeout) const {
    const std::string node = _impl ? _impl->location : std::string{};
    return Guarded(ErrorCode::io_error, node, [&]() -> Result<DatasetSize> {
        if (!_impl) {
            return Error{ErrorCode::invalid_argument, "Dataset handle is empty"};
        }

        std::uint64_t physical_size = 0;
        if (TryComputeDirectorySize(_impl->location, directory_size_timeout, physical_size)) {
            return DatasetSize{physical_size, false};
        }

        auto logical_size = internal::TotalArraySizeBytes(*_impl->store);
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
            return Error{ErrorCode::invalid_argument, "Dataset handle is empty"};
        }
        std::scoped_lock const lock(_impl->mutex);
        const std::string image_name(image_id);
        // Asked of the listing this dataset kept, rather than by enumerating the store again. It is
        // the same question and the same answer the profile would give, which is why it is the
        // profile's function and not a second copy of the rule here.
        if (auto openable = internal::RequireOpenable(_impl->descriptor.images, image_name); !openable) {
            return openable.error();
        }
        const auto make_image = [&](const ImageDescriptor& descriptor) {
            // The descriptor already carries the stored layout; the geometry is that layout permuted
            // into logical order, so it is derived here rather than read again.
            const StorageLayout layout = descriptor.storage ? *descriptor.storage : StorageLayout{};
            return Image{std::make_shared<Image::Impl>(_impl->context, _impl->location, _impl->profile, _impl->store,
                                                       descriptor, internal::BuildChunkGeometry(descriptor, layout))};
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
    // The guard that does not report an Error, because a probe answers with a ProbeResult whatever
    // happens: "not a Zarr store" is an answer rather than a failure. A throw becomes a diagnostic,
    // and that is the only way this differs from every other entry point here.
    const auto as_diagnostic = [](const auto& failure) {
        return Diagnostic{internal::zarr::ErrorCodeName(failure.code), failure.message, failure.node_path};
    };
    return GuardedWith(
        [&]() -> ProbeResult {
            auto store_result = internal::OpenStore(location);
            if (!store_result) {
                const auto& failure = store_result.error();
                ProbeResult result;
                result.kind = failure.code == ErrorCode::invalid_metadata || failure.code == ErrorCode::io_error
                                  ? ProbeKind::invalid_dataset
                                  : ProbeKind::not_zarr;
                result.diagnostics.push_back(as_diagnostic(failure));
                return result;
            }
            auto probe_result = internal::ProbeStore(store_result.value());
            if (!probe_result) {
                ProbeResult result;
                result.kind = ProbeKind::invalid_dataset;
                result.diagnostics.push_back(as_diagnostic(probe_result.error()));
                return result;
            }
            return probe_result.value();
        },
        [&](const std::exception& error) {
            ProbeResult result;
            result.kind = ProbeKind::invalid_dataset;
            result.diagnostics.push_back(Diagnostic{"exception", error.what(), std::string(location)});
            return result;
        });
}

Result<SchemaProbeResult> ProbeSchema(std::string_view location, std::string_view schema_id) {
    return Guarded(ErrorCode::invalid_metadata, std::string(location), [&]() -> Result<SchemaProbeResult> {
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
    });
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
        return Error{ErrorCode::invalid_metadata,
                     internal::RejectionMessage(result.value().diagnostics,
                                                "The requested schema did not match"),
                     std::string(location)};
    }
    return false;
}

}  // namespace carta::zarr
