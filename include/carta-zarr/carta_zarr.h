/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_CARTA_ZARR_H_
#define CARTA_ZARR_CARTA_ZARR_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/export.h"
#include "carta-zarr/read.h"
#include "carta-zarr/reduce.h"
#include "carta-zarr/result.h"

#include <chrono>
#include <cstdint>
#include <memory>

namespace carta::zarr {

class CARTA_ZARR_EXPORT Context final {
public:
    Context(const Context&) = default;
    Context& operator=(const Context&) = default;
    Context(Context&&) noexcept = default;
    Context& operator=(Context&&) noexcept = default;
    ~Context();

    static Result<Context> Create(const OpenOptions& options = {});

private:
    class Impl;
    explicit Context(std::shared_ptr<Impl> impl);

    std::shared_ptr<Impl> _impl;

    friend class Dataset;
    friend class Image;
};

class CARTA_ZARR_EXPORT Image final {
public:
    Image(const Image&) = default;
    Image& operator=(const Image&) = default;
    Image(Image&&) noexcept = default;
    Image& operator=(Image&&) noexcept = default;
    ~Image();

    const ImageDescriptor& descriptor() const noexcept;

    // The read geometry, in the same axis order as descriptor().axes.
    const ChunkGeometry& chunk_geometry() const noexcept;

    // Reads a densely packed result in logical axis order, axis 0 fastest-varying. Returns the
    // number of elements written. Safe to call concurrently on one handle. On failure, the
    // destination may be unchanged, partially written, or fully written; callers must discard it.
    //
    // float is the only output this library produces. It is said in the destination's type rather
    // than asked for in the request, because a request that could name a type the buffer was not
    // shaped for is a mistake worth making unspellable.
    Result<std::size_t> Read(const ReadRequest& request, BufferView<float> destination) const;
    Result<std::size_t> Read(const ReadRequest& request, BufferView<float> destination,
                             const ReadOptions& options) const;
    // Watched as it advances, which also splits it into pieces. See ProgressCallback.
    Result<std::size_t> Read(const ReadRequest& request, BufferView<float> destination,
                             const ReadOptions& options, const ProgressCallback& progress) const;

    // Reads this image's pixel mask over the same region, one byte per pixel, true meaning a good
    // pixel. Reports not_found when the image has no mask.
    //
    // Takes a ReadControl rather than a ReadOptions, which is the whole of what it used to honour:
    // apply_pixel_mask means nothing here because this is the mask, and the read is issued in one
    // piece -- the destination is the caller's, so there is no temporary of ours for a ceiling to
    // bound and nowhere to report from. Those three used to be fields a caller could set and this
    // would quietly ignore; now they are not fields it can be handed.
    Result<std::size_t> ReadPixelMask(const ReadRequest& request, BufferView<std::uint8_t> destination) const;
    Result<std::size_t> ReadPixelMask(const ReadRequest& request, BufferView<std::uint8_t> destination,
                                      const ReadControl& control) const;

    // Reduces every region over the same channels in one pass over the pixels, handing results to
    // the sink block by block.
    //
    // The pass visits each covered chunk once and accumulates every region that touches it, which
    // is the whole point of taking N regions instead of being called N times: a position-velocity
    // cut along the diagonal of a 4096^2 image is 5,792 overlapping boxes, and reducing them one at
    // a time decompresses the same chunks thousands of times over.
    //
    // Results stream rather than accumulate: those 5,792 regions over 30,000 channels would be
    // 2.59 GiB returned at once. Each block is valid only inside the sink call.
    //
    // The image's pixel mask is applied when it has one and ReadOptions::apply_pixel_mask is left
    // on, exactly as it is for a read: a flagged pixel reaches the statistics as NaN, and declining
    // the mask means the flag is never read at all. ReadOptions also supplies cancellation, the
    // deadline, and a ceiling on the pixel buffer the pass may hold.
    Result<void> ReduceSpectral(const SpectralReduceRequest& request, const SpectralSink& sink) const;
    Result<void> ReduceSpectral(const SpectralReduceRequest& request, const SpectralSink& sink,
                                const ReadOptions& options) const;

    // Bins every pixel of each plane over a fixed range, handing counts to the sink block by block.
    //
    // Separate from ReduceSpectral because a histogram is not one of the statistics that reduction
    // accumulates, and because it needs none of that machinery: the region is always the whole
    // plane, so there is nothing to index and no region raster to consult. The image's own pixel
    // mask is a different thing and still applies -- see below.
    //
    // The image's pixel mask is applied when it has one and ReadOptions::apply_pixel_mask is left
    // on, so a flagged pixel is not counted -- the same thing that happens to a NaN. Declining the
    // mask counts every stored pixel, flagged or not.
    Result<void> ComputeHistogram(const HistogramRequest& request, const HistogramSink& sink) const;
    Result<void> ComputeHistogram(const HistogramRequest& request, const HistogramSink& sink,
                                  const ReadOptions& options) const;

    // One histogram for the whole selection in a single pass, when the range is not known in
    // advance. See CubeHistogramRequest for what that costs and what it keeps exact. The pixel mask
    // is applied on the same terms as ComputeHistogram, and on both passes it makes.
    Result<CubeHistogramResult> ComputeCubeHistogram(const CubeHistogramRequest& request) const;
    Result<CubeHistogramResult> ComputeCubeHistogram(const CubeHistogramRequest& request,
                                                     const ReadOptions& options) const;

    Result<std::vector<Beam>> ReadBeams() const;

private:
    class Impl;
    explicit Image(std::shared_ptr<Impl> impl);

    std::shared_ptr<Impl> _impl;

    friend class Dataset;
};

class CARTA_ZARR_EXPORT Dataset final {
public:
    Dataset(const Dataset&) = default;
    Dataset& operator=(const Dataset&) = default;
    Dataset(Dataset&&) noexcept = default;
    Dataset& operator=(Dataset&&) noexcept = default;
    ~Dataset();

    static Result<Dataset> Open(const Context& context, std::string_view location);

    const DatasetDescriptor& descriptor() const noexcept;
    // How much room this dataset takes where it is stored, and whether that number was measured or
    // inferred. Measured when the store can be sized within the timeout; otherwise the logical
    // uncompressed size of all arrays, marked as an upper bound.
    //
    // The timeout is not named after a directory because a dataset need not live in one: what can
    // be sized, and how quickly, is the transport's affair.
    Result<DatasetSize> Size(
        std::chrono::milliseconds stored_size_timeout = std::chrono::milliseconds(50)) const;
    Result<Image> OpenImage(std::string_view image_id) const;

private:
    class Impl;
    explicit Dataset(std::shared_ptr<Impl> impl);

    std::shared_ptr<Impl> _impl;
};

CARTA_ZARR_EXPORT ProbeResult Probe(std::string_view location, const ProbeOptions& options = {});

CARTA_ZARR_EXPORT Result<SchemaProbeResult> ProbeSchema(std::string_view location, std::string_view schema_id);

// Returns an error for an unreadable or malformed store; false is a valid non-match.
CARTA_ZARR_EXPORT Result<bool> IsXradioImage(std::string_view location);

}  // namespace carta::zarr

#endif  // CARTA_ZARR_CARTA_ZARR_H_
