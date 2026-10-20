/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The threading promises docs/concepts.md and carta_zarr.h make, held through the public API: a Dataset
// opened from several threads at once, reductions and histograms running at once on one image and on
// two images of one context, and Read beside Prefetch on one handle.
//
// Each test computes its answers serially first, on a context of their own, and then has several
// threads ask the same questions at once of a fresh context, released together so that they overlap
// rather than take turns. A fresh context and dataset matter: what is computed lazily and memoised is
// then computed under contention, which is where a race would be.
//
// Agreement with the serial answer is what this can check. That no two threads touched the same
// memory unsynchronised is what ThreadSanitizer checks, and this test is written to be run under it,
// configured with CARTA_ZARR_THREAD_SANITIZER: see "Run the tests" in docs/installation.md.

#include "support/check.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include <carta-zarr/carta_zarr.h>

namespace {

using carta::zarr::testing::Require;

constexpr int kThreads = 8;

// Enough decode threads that the reductions on the wide fixture split their work, and few enough that
// eight callers sharing them have to queue.
carta::zarr::ContextOptions SharedPool() {
    carta::zarr::ContextOptions options;
    options.decode_threads = 4;
    options.cache_bytes = std::size_t{64} << 20;
    return options;
}

carta::zarr::Context NewContext() {
    auto context = carta::zarr::Context::Create(SharedPool());
    Require(static_cast<bool>(context), "Context::Create failed");
    return context.value();
}

carta::zarr::Dataset OpenDataset(const carta::zarr::Context& context, const char* fixture) {
    Require(std::filesystem::exists(fixture),
            std::string("a fixture is missing: ") + fixture + "; run tests/data/generate_zarr_fixtures.py");
    auto dataset = carta::zarr::Dataset::Open(context, fixture);
    Require(static_cast<bool>(dataset), std::string("Dataset::Open failed on ") + fixture);
    return dataset.value();
}

carta::zarr::Image OpenSky(const carta::zarr::Context& context, const char* fixture) {
    auto image = OpenDataset(context, fixture).OpenImage("SKY");
    Require(static_cast<bool>(image), std::string("SKY could not be opened in ") + fixture);
    return image.value();
}

// Runs body(t) on kThreads threads, released together, and fails with the first thread's failure.
//
// Each thread counts itself in and waits until all have, so that none has finished before the last
// has started; without that, eight short calls on eight new threads can run one after another and
// the test checks nothing a serial loop would not.
void RunTogether(const std::function<void(int)>& body) {
    std::atomic<int> arrived{0};
    std::mutex failure_mutex;
    std::string failure;
    std::vector<std::thread> threads;
    threads.reserve(kThreads);
    for (int t = 0; t < kThreads; ++t) {
        threads.emplace_back([&, t] {
            arrived.fetch_add(1);
            while (arrived.load() < kThreads) {
                std::this_thread::yield();
            }
            try {
                body(t);
            } catch (const std::exception& error) {
                const std::lock_guard<std::mutex> lock(failure_mutex);
                if (failure.empty()) {
                    failure = "thread " + std::to_string(t) + ": " + error.what();
                }
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    Require(failure.empty(), failure);
}

// Everything a descriptor says, written out, so that two can be compared whole. The library gives
// descriptors no operator==, and a test that compared a few fields would pass a race that tore any of
// the others.
void Write(std::ostream& out, const std::vector<carta::zarr::Diagnostic>& diagnostics) {
    for (const auto& diagnostic : diagnostics) {
        out << "diagnostic " << carta::zarr::DiagnosticCodeName(diagnostic.code) << " " << diagnostic.node_path << " "
            << diagnostic.message << "\n";
    }
}

void Write(std::ostream& out, const std::vector<double>& values) {
    for (const double value : values) {
        out << value << " ";
    }
    out << "\n";
}

std::string Fingerprint(const carta::zarr::ImageDescriptor& image) {
    std::ostringstream out;
    out << std::setprecision(17);
    out << image.id << " " << image.sub_type << " " << static_cast<int>(image.stored_type) << " " << image.unit << " "
        << image.pixel_mask_id << "\n";
    for (const auto& group : image.data_groups) {
        out << "group " << group << "\n";
    }
    for (const auto& axis : image.axes) {
        out << "axis " << axis.name << " " << static_cast<int>(axis.role) << " " << axis.length << " "
            << axis.storage_index << "\n";
    }
    if (const auto& direction = image.direction) {
        out << "direction " << direction->projection << " " << direction->reference_frame << " "
            << direction->equinox.value_or(-1.0);
        for (int i = 0; i < 2; ++i) {
            out << " " << direction->reference_pixel[i] << " " << direction->reference_value[i] << " "
                << direction->increment[i] << " " << direction->native_pole_direction[i] << " "
                << direction->transformation_matrix[i][0] << " " << direction->transformation_matrix[i][1];
        }
        out << "\n";
        Write(out, direction->projection_parameters);
    }
    if (const auto& spectral = image.spectral) {
        out << "spectral " << spectral->unit << " " << spectral->reference_frame << " "
            << spectral->reference_pixel.value_or(-1.0) << " " << spectral->reference_value.value_or(-1.0) << " "
            << spectral->increment.value_or(-1.0) << " " << spectral->rest_frequency.value_or(-1.0) << "\n";
        Write(out, spectral->values);
    }
    if (const auto& polarization = image.polarization) {
        for (const auto& label : polarization->labels) {
            out << "stokes " << label << "\n";
        }
    }
    if (const auto& temporal = image.temporal) {
        out << "time " << temporal->unit << " " << temporal->scale << " " << temporal->format << "\n";
        Write(out, temporal->values);
    }
    if (const auto& observation = image.observation) {
        out << "observation " << observation->object_name << " " << observation->observer << " "
            << observation->telescope_name << " " << observation->time_scale << " "
            << observation->observation_mjd.value_or(-1.0);
        if (const auto& position = observation->observatory_position) {
            out << " " << (*position)[0] << " " << (*position)[1] << " " << (*position)[2];
        }
        out << "\n";
    }
    Write(out, image.diagnostics);
    return out.str();
}

// What OpenImage answered for one id: the descriptor, or the error.
std::string Outcome(const carta::zarr::Result<carta::zarr::Image>& opened) {
    if (!opened) {
        return "error " + std::to_string(static_cast<int>(opened.error().code)) + " " + opened.error().node_path + " " +
               opened.error().message;
    }
    return Fingerprint(opened.value().descriptor());
}

// Every image the dataset lists, opened at once from one Dataset: first every thread the same id, then
// each a different one in its own order. The fixture lists images that open and images that do not,
// and a refusal has to come out the same as well.
void TestOpeningImagesAtOnce(const char* fixture) {
    std::vector<std::string> ids;
    std::map<std::string, std::string> serial;
    {
        const auto dataset = OpenDataset(NewContext(), fixture);
        for (const auto& entry : dataset.descriptor().images) {
            ids.push_back(entry.id);
            serial[entry.id] = Outcome(dataset.OpenImage(entry.id));
        }
    }
    Require(ids.size() >= 3, "the fixture should list several images, or different ids cannot overlap");

    for (const auto& id : ids) {
        const auto dataset = OpenDataset(NewContext(), fixture);
        RunTogether([&](int) {
            Require(Outcome(dataset.OpenImage(id)) == serial.at(id),
                    "opening " + id + " from several threads disagreed with opening it alone");
        });
    }

    const auto dataset = OpenDataset(NewContext(), fixture);
    RunTogether([&](int t) {
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const auto& id = ids[(i + static_cast<std::size_t>(t)) % ids.size()];
            Require(Outcome(dataset.OpenImage(id)) == serial.at(id),
                    "opening " + id + " beside other images disagreed with opening it alone");
        }
    });
}

// Every channel of one polarization, so that a reduction over it reads every chunk of that polarization.
carta::zarr::PlaneSelection EveryChannel(const carta::zarr::Image& image, std::uint64_t polarization) {
    const auto& axes = image.descriptor().axes;
    const auto spectral = carta::zarr::AxisIndex(axes, carta::zarr::AxisRole::spectral);
    Require(spectral.has_value(), "the fixture's image should have a spectral axis");
    carta::zarr::PlaneSelection planes;
    planes.spectral = {0, axes.at(*spectral).length, 1};
    planes.polarization = polarization;
    return planes;
}

std::uint64_t Length(const carta::zarr::Image& image, carta::zarr::AxisRole role) {
    const auto& axes = image.descriptor().axes;
    return axes.at(carta::zarr::AxisIndex(axes, role).value()).length;
}

// One reduction, every block that completes, as region -> statistic -> channel.
using Statistics = std::map<std::tuple<std::size_t, carta::zarr::Statistic, std::uint64_t>, double>;

const carta::zarr::Statistic kStatistics[]{
    carta::zarr::Statistic::num_pixels, carta::zarr::Statistic::nan_count, carta::zarr::Statistic::sum,
    carta::zarr::Statistic::sum_sq,     carta::zarr::Statistic::min,       carta::zarr::Statistic::max,
    carta::zarr::Statistic::sum_sq_dev,
};

Statistics Reduce(const carta::zarr::Image& image) {
    const std::uint64_t x = Length(image, carta::zarr::AxisRole::spatial_x);
    const std::uint64_t y = Length(image, carta::zarr::AxisRole::spatial_y);
    const std::vector<carta::zarr::RegionMask> regions{
        {0, 0, x, y},
        {x / 4, y / 5, x / 2, y / 2},
    };
    carta::zarr::SpectralReduceRequest request;
    request.planes = EveryChannel(image, 1);
    request.regions = {regions.data(), regions.size()};
    for (const auto statistic : kStatistics) {
        request.statistics |= statistic;
    }

    Statistics collected;
    const auto result = image.ReduceSpectral(request, [&](const carta::zarr::SpectralBlock& block) {
        if (!block.complete) {
            return true;
        }
        for (std::size_t r = 0; r < block.region_count; ++r) {
            for (const auto statistic : kStatistics) {
                const double* series = block.Series(r, statistic);
                Require(series != nullptr, "a block should carry every statistic asked for");
                for (std::uint64_t c = 0; c < block.channel_count; ++c) {
                    collected[{r, statistic, block.first_channel + c}] = series[c];
                }
            }
        }
        return true;
    });
    Require(static_cast<bool>(result),
            std::string("ReduceSpectral failed: ") + (result.has_value() ? "" : result.error().message));
    // Checked here so that two reductions which both lost their blocks cannot agree with each other.
    Require(collected.size() == regions.size() * std::size(kStatistics) * request.planes.spectral.count,
            "ReduceSpectral should report every region, statistic and channel once complete");
    return collected;
}

// The counts and extremes are exact whatever the order the work was done in. The sums are added in an
// order the reduction's workers decide, so they are held to a rounding: what a race would do to them
// is lose or double a whole chunk's worth, which is far outside it.
bool Exact(carta::zarr::Statistic statistic) {
    return statistic == carta::zarr::Statistic::num_pixels || statistic == carta::zarr::Statistic::nan_count ||
           statistic == carta::zarr::Statistic::min || statistic == carta::zarr::Statistic::max;
}

bool Agree(double actual, double expected, bool exact) {
    if (std::isnan(expected)) {
        return std::isnan(actual);
    }
    if (exact) {
        return actual == expected;
    }
    return std::abs(actual - expected) <= 1e-9 * (1.0 + std::abs(expected));
}

void RequireSame(const Statistics& actual, const Statistics& expected, const std::string& what) {
    Require(actual.size() == expected.size(), what + ": a different number of values came back");
    for (const auto& [key, value] : expected) {
        const auto found = actual.find(key);
        Require(found != actual.end(), what + ": a value went missing");
        Require(Agree(found->second, value, Exact(std::get<1>(key))),
                what + ": region " + std::to_string(std::get<0>(key)) + " statistic " +
                    std::to_string(static_cast<int>(std::get<1>(key))) + " channel " +
                    std::to_string(std::get<2>(key)) + " was " + std::to_string(found->second) + ", alone " +
                    std::to_string(value));
    }
}

// Fine enough that each channel's pixels spread over many bins. Over the range of the whole cube, as
// these are, a coarse grid gives each channel of the wide fixture one bin, and a comparison of those
// counts checks only how many pixels each channel has.
constexpr std::uint32_t kHistogramBins = 4096;

// Plane histograms over every channel, as channel -> counts, over the range a cube histogram found.
std::vector<std::vector<std::uint64_t>> Histogram(const carta::zarr::Image& image, double lower, double upper) {
    carta::zarr::HistogramRequest request;
    request.planes = EveryChannel(image, 0);
    request.bin_count = kHistogramBins;
    request.lower = lower;
    request.upper = upper;
    std::vector<std::vector<std::uint64_t>> counts(request.planes.spectral.count);
    const auto result = image.ComputeHistogram(request, [&](const carta::zarr::HistogramBlock& block) {
        if (!block.complete) {
            return true;
        }
        for (std::uint64_t c = 0; c < block.channel_count; ++c) {
            const auto* row = block.Counts(c);
            counts.at(block.first_channel + c).assign(row, row + block.bin_count);
        }
        return true;
    });
    Require(static_cast<bool>(result),
            std::string("ComputeHistogram failed: ") + (result.has_value() ? "" : result.error().message));
    for (const auto& row : counts) {
        Require(row.size() == kHistogramBins, "ComputeHistogram should report every channel once complete");
    }
    return counts;
}

carta::zarr::CubeHistogramResult CubeHistogram(const carta::zarr::Image& image) {
    carta::zarr::CubeHistogramRequest request;
    request.planes = EveryChannel(image, 0);
    request.bin_count = 128;
    auto result = image.ComputeCubeHistogram(request);
    Require(static_cast<bool>(result),
            std::string("ComputeCubeHistogram failed: ") + (result.has_value() ? "" : result.error().message));
    return result.value();
}

// The bins a cube histogram fills are fixed by how its walk is split, which a context's thread count
// decides and the order the tasks run in does not: each task bins into a histogram of its own. So two
// runs on contexts of one size agree bin for bin, at once or not.
void RequireSame(const carta::zarr::CubeHistogramResult& actual, const carta::zarr::CubeHistogramResult& expected,
                 const std::string& what) {
    const auto& a = actual.totals;
    const auto& e = expected.totals;
    Require(Agree(a.num_pixels, e.num_pixels, true) && Agree(a.nan_count, e.nan_count, true) &&
                Agree(a.min, e.min, true) && Agree(a.max, e.max, true),
            what + ": the cube histogram's counts or extremes differ from a run alone");
    Require(Agree(a.sum, e.sum, false) && Agree(a.sum_sq, e.sum_sq, false) && Agree(a.sum_sq_dev, e.sum_sq_dev, false),
            what + ": the cube histogram's sums differ from a run alone by more than a rounding");
    Require(actual.counts == expected.counts, what + ": the cube histogram's bins differ from a run alone");
    Require(actual.sampled == expected.sampled,
            what + ": the cube histogram says it was sampled when alone it was not");
}

// What each image answers alone, computed on a context of its own.
struct Answers {
    Statistics reduced;
    carta::zarr::CubeHistogramResult cube;
    std::vector<std::vector<std::uint64_t>> histogram;
};

Answers Alone(const char* fixture) {
    const auto image = OpenSky(NewContext(), fixture);
    Answers answers{Reduce(image), CubeHistogram(image), {}};
    answers.histogram = Histogram(image, answers.cube.totals.min, answers.cube.totals.max);
    for (const auto& row : answers.histogram) {
        const auto filled = std::count_if(row.begin(), row.end(), [](std::uint64_t count) { return count > 0; });
        Require(filled >= 4, std::string("a channel of ") + fixture + " fills " + std::to_string(filled) +
                                 " bins, too few for its counts to say where its pixels are");
    }
    return answers;
}

// Thread t runs operation t % 3 on image (t / 3) % images.size(), twice over, so every operation runs
// beside the other two and beside itself, on each image.
void RunReductionsTogether(const std::vector<carta::zarr::Image>& images, const std::vector<Answers>& alone,
                           const std::string& what) {
    RunTogether([&](int t) {
        const auto which = static_cast<std::size_t>(t / 3) % images.size();
        const auto& image = images.at(which);
        const auto& expected = alone.at(which);
        const std::string where = what + ", image " + std::to_string(which);
        for (int repeat = 0; repeat < 2; ++repeat) {
            switch (t % 3) {
                case 0:
                    RequireSame(Reduce(image), expected.reduced, where + ", ReduceSpectral");
                    break;
                case 1:
                    Require(Histogram(image, expected.cube.totals.min, expected.cube.totals.max) == expected.histogram,
                            where + ", ComputeHistogram: the counts differ from a run alone");
                    break;
                default:
                    RequireSame(CubeHistogram(image), expected.cube, where + ", ComputeCubeHistogram");
                    break;
            }
        }
    });
}

// ReduceSpectral, ComputeHistogram and ComputeCubeHistogram at once on one image, all sharing the
// context's worker pool. The wide fixture is large enough that each splits its work between the
// pool's workers, so the pool is running several callers' tasks at once rather than one at a time.
void TestReductionsAtOnceOnOneImage(const char* fixture) {
    const std::vector<Answers> alone{Alone(fixture)};
    const auto image = OpenSky(NewContext(), fixture);
    RunReductionsTogether({image}, alone, "on one image");
}

// The same on two images of one context: they share its pool and its cache, and nothing else. One
// of them has a pixel mask, so the flag is read beside the pixels on some threads and not on others.
void TestReductionsAtOnceOnTwoImages(const char* first, const char* second) {
    const std::vector<Answers> alone{Alone(first), Alone(second)};
    const auto context = NewContext();
    RunReductionsTogether({OpenSky(context, first), OpenSky(context, second)}, alone, "on two images");
}

carta::zarr::ReadRequest WholeImage(const carta::zarr::Image& image) {
    carta::zarr::ReadRequest request;
    for (const auto& axis : image.descriptor().axes) {
        request.axes.push_back(carta::zarr::Range{0, axis.length, 1});
    }
    return request;
}

std::vector<float> ReadWhole(const carta::zarr::Image& image) {
    const auto request = WholeImage(image);
    std::size_t elements = 1;
    for (const auto& range : request.axes) {
        elements *= range.count;
    }
    std::vector<float> pixels(elements);
    const auto read = image.Read(request, {pixels.data(), pixels.size()});
    Require(static_cast<bool>(read), std::string("Read failed: ") + (read.has_value() ? "" : read.error().message));
    Require(read.value() == elements, "Read reported the wrong element count");
    return pixels;
}

// A NaN is what a flagged pixel and an absent chunk read as, and it has to come back as one.
bool SamePixels(const std::vector<float>& actual, const std::vector<float>& expected) {
    if (actual.size() != expected.size()) {
        return false;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const bool same = std::isnan(expected[i]) ? std::isnan(actual[i]) : actual[i] == expected[i];
        if (!same) {
            return false;
        }
    }
    return true;
}

// Read and Prefetch at once on one handle, which carta_zarr.h declares safe. Half the threads prefetch
// the whole image while the other half read it, into one shared cache, so a read finds chunks that a
// prefetch has decoded, is decoding, or has not reached. Every read has to return what a read alone
// does, and every prefetch has to decode what one alone counts.
void TestReadAndPrefetchAtOnce(const char* fixture) {
    std::vector<float> serial_pixels;
    std::uint64_t serial_chunks = 0;
    {
        const auto image = OpenSky(NewContext(), fixture);
        const auto prefetched = image.Prefetch(WholeImage(image));
        Require(static_cast<bool>(prefetched), "Prefetch failed on its own");
        serial_chunks = prefetched.value();
        serial_pixels = ReadWhole(image);
    }
    Require(serial_chunks > 1, "the fixture should be more than one chunk, or there is nothing to overlap");

    const auto image = OpenSky(NewContext(), fixture);
    RunTogether([&](int t) {
        for (int repeat = 0; repeat < 3; ++repeat) {
            if (t % 2 == 0) {
                const auto prefetched = image.Prefetch(WholeImage(image));
                Require(static_cast<bool>(prefetched), "Prefetch failed beside reads");
                Require(prefetched.value() == serial_chunks, "Prefetch beside reads counted different chunks");
            } else {
                Require(SamePixels(ReadWhole(image), serial_pixels),
                        "a Read beside prefetches disagreed with a read alone");
            }
        }
    });
}

}  // namespace

int main() {
    const std::pair<const char*, std::function<void()>> tests[]{
        {"opening images at once", [] { TestOpeningImagesAtOnce(CARTA_ZARR_REFERENCE_FIXTURE); }},
        {"reductions at once on one image", [] { TestReductionsAtOnceOnOneImage(CARTA_ZARR_PIXEL_FIXTURE_WIDE); }},
        {"reductions at once on two images",
         [] { TestReductionsAtOnceOnTwoImages(CARTA_ZARR_PIXEL_FIXTURE_WIDE, CARTA_ZARR_PIXEL_FIXTURE); }},
        {"Read and Prefetch at once on the wide fixture",
         [] { TestReadAndPrefetchAtOnce(CARTA_ZARR_PIXEL_FIXTURE_WIDE); }},
        {"Read and Prefetch at once on a flagged fixture", [] { TestReadAndPrefetchAtOnce(CARTA_ZARR_PIXEL_FIXTURE); }},
    };
    for (const auto& [name, test] : tests) {
        try {
            test();
        } catch (const std::exception& error) {
            std::cerr << "concurrency test failed, " << name << ": " << error.what() << "\n";
            return 1;
        }
    }
    return 0;
}
