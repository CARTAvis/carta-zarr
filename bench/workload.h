/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_BENCH_WORKLOAD_H_
#define CARTA_ZARR_BENCH_WORKLOAD_H_

// Where each operation of a trial reads, and the reading itself.
//
// Positions are a function of the seed, the mode, the trial and the cube's logical shape -- never of
// its layout -- so every layout of one cube is asked for the same pixels, and a difference in time is
// the layout's. They come from a generator written out here rather than from <random>, whose
// distributions are free to differ between standard libraries, so the Mac and a Linux server agree.

#include "options.h"

#include <carta-zarr/carta_zarr.h>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace carta::zarr::bench {

// The axes the operations address, as logical indices into ImageDescriptor::axes. An image with no
// polarization or time axis reads as one with a single plane of each.
struct CubeAxes {
    std::size_t rank = 0;
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t spectral = 0;
    std::optional<std::size_t> polarization;
    std::optional<std::size_t> time;

    std::uint64_t width = 1;
    std::uint64_t height = 1;
    std::uint64_t channels = 1;
    std::uint64_t polarizations = 1;

    // The image's axes, or an error naming the role it lacks.
    static Result<CubeAxes> Of(const ImageDescriptor& descriptor);
};

// One operation, as planned before the trial starts.
struct Operation {
    Mode mode = Mode::plane;
    std::uint64_t x = 0;
    std::uint64_t y = 0;
    std::uint64_t width = 0;
    std::uint64_t height = 0;
    std::uint64_t channel = 0;
    std::uint64_t channel_count = 0;
    std::uint64_t polarization = 0;
    // Whether another operation of this trial, in any process, reads the same position: there were
    // more operations than distinct positions to give them.
    bool overlap = false;

    // Where it reads, for the CSV: semicolon-separated, so that it stays one field.
    std::string Describe() const;
};

// The operations one process makes in one trial. Process p takes the p-th run of `ops` positions from
// one sequence drawn for the whole trial, so processes never share one while there are enough, and
// process 0 reads the same positions whatever the process count.
//
// cube-histogram is the exception, since its operation covers the plane: process p takes the p-th of
// `processes` contiguous runs of channels, and each operation a polarization of its own.
//
// A region box covers `region_fraction` of the plane, square in pixels' proportion to the plane, and
// sits in a cell of a grid of as many such boxes as fit along each side.
std::vector<Operation> PlanOperations(Mode mode, const CubeAxes& axes, std::uint64_t seed, unsigned trial,
                                      unsigned processes, unsigned process_index, unsigned ops,
                                      double region_fraction = 0.05);

// Runs operations against one image and remembers enough of the last result to fingerprint it.
//
// The fingerprint is taken after the clock stops, and only of what every layout of the same pixels
// must agree on exactly: the pixels a read returns, with every NaN one NaN; and for a reduction or a
// histogram the pixel counts and extremes, which are exact, but not sums, whose rounding depends on
// the order the chunks were visited in. An exact histogram's counts are in it too, since binning
// over fixed bounds is exact; a one-pass histogram's are not, since ComputeCubeHistogram says they
// depend on the thread count.
//
// A cube histogram reads through a cache pool that keeps nothing, as carta-backend's cube walks do,
// so that a scan of the whole cube neither fills nor finds the context's cache.
class Runner {
public:
    // For every mode but open, which brings its own context and image to each operation.
    Runner(const Context& context, Image image, HistogramMethod histogram = {});
    // For open.
    Runner(ContextOptions context, std::string dataset, std::string image_id);

    // How many elements of the cube the operation covered.
    Result<std::uint64_t> Run(const Operation& operation, const ReadOptions& options);

    std::uint64_t Fingerprint() const;

    // The image last read, if there was one: an open that failed leaves none.
    const std::optional<Image>& image() const noexcept {
        return _image;
    }

private:
    Result<std::uint64_t> Read(const Operation& operation, const ReadOptions& options);
    Result<std::uint64_t> Reduce(const Operation& operation, const ReadOptions& options);
    Result<std::uint64_t> Histogram(const Operation& operation, const ReadOptions& options);
    Result<std::uint64_t> ExactHistogram(const Operation& operation, const ReadOptions& options);
    Result<std::uint64_t> Open();

    std::optional<Image> _image;
    std::optional<CubeAxes> _axes;
    HistogramMethod _histogram;
    std::optional<CachePool> _keeping_nothing;
    ContextOptions _context;
    std::string _dataset;
    std::string _image_id;

    std::vector<float> _pixels;
    std::size_t _pixel_count = 0;
    // The exact statistics of a reduction or a cube histogram, in order.
    std::vector<double> _exact;
    std::vector<std::uint64_t> _counts;
};

// FNV-1a over a run of values, with every NaN hashed as the one quiet NaN.
std::uint64_t Fingerprint(const float* values, std::size_t count);
std::uint64_t Fingerprint(const double* values, std::size_t count);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_WORKLOAD_H_
