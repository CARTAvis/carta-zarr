/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// An ordinary read against pixels that were never written down.
//
// tests/piece_plan_test.cc covers what a read decides before it reads anything. This covers what it
// does once it starts, which until ReadInPieces took a PixelSource could only be asked of a
// directory tree -- and three of the four things below could not be asked at all, because a fixture
// on disk cannot be told to fail one read and not another.
//
// Links no Store and no TensorStore.

#include "read/pieces.h"

#include "support/check.h"

#include "support/synthetic_pixel_source.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

using carta::zarr::AxisRole;
using carta::zarr::BufferView;
using carta::zarr::ChunkGeometry;
using carta::zarr::ErrorCode;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ProgressCallback;
using carta::zarr::ReadControl;
using carta::zarr::ReadOptions;
using carta::zarr::ReadRequest;
using carta::zarr::internal::ReadableImage;
using carta::zarr::internal::ReadInPieces;
using carta::zarr::internal::ReadPixelMask;
using carta::zarr::internal::WorkPool;
using carta::zarr::testing::SyntheticPixelSource;

using carta::zarr::testing::Require;

// A ReadableImage names a pool because a reduction runs its arithmetic on one. An ordinary read
// never touches it, and a pool of one starts no threads at all, so this costs nothing here.
WorkPool& InlinePool() {
    static WorkPool pool(1);
    return pool;
}

ReadableImage Readable(const SyntheticPixelSource& source, const ImageDescriptor& image,
                       const ChunkGeometry& geometry) {
    auto readable = ReadableImage::Of(source, image, geometry, InlinePool());
    Require(static_cast<bool>(readable), "the synthetic image's axes could not be mapped");
    return readable.value();
}

constexpr std::uint64_t kX = 64;
constexpr std::uint64_t kY = 40;
constexpr std::uint64_t kZ = 6;

ImageDescriptor MakeImage(bool with_mask = false) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, kX},
             {"m", AxisRole::spatial_y, kY},
             {"frequency", AxisRole::spectral, kZ},
             {"polarization", AxisRole::polarization, 1},
             {"time", AxisRole::time, 1}};
    // The store varies m fastest, so the read is the transposing one -- which is the case worth
    // walking, because the destination's contract is logical order whatever the store did.
    const std::size_t storage[]{3, 4, 1, 2, 0};
    for (std::size_t i = 0; i < 5; ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = axes[i].name;
        axis.role = axes[i].role;
        axis.length = axes[i].length;
        axis.storage_index = storage[i];
        descriptor.axes.push_back(axis);
    }
    descriptor.has_pixel_mask = with_mask;
    if (with_mask) {
        descriptor.pixel_mask_id = "FLAG";
    }
    return descriptor;
}

ChunkGeometry MakeGeometry() {
    ChunkGeometry geometry;
    geometry.fastest_spatial_axis = AxisRole::spatial_y;
    geometry.chunk_shape = {16, 20, 2, 1, 1};
    return geometry;
}

// Distinct for every pixel of the cube, so a value landing at the wrong destination offset is a
// different number rather than a coincidence.
float Value(const std::vector<std::uint64_t>& logical) {
    return static_cast<float>((logical.at(2) * kX * kY) + (logical.at(1) * kX) + logical.at(0));
}

ReadRequest WholeCube() {
    ReadRequest request;
    request.axes = {Range{0, kX, 1}, Range{0, kY, 1}, Range{0, kZ, 1}, Range{0, 1, 1}, Range{0, 1, 1}};
    return request;
}

constexpr std::size_t kElements = static_cast<std::size_t>(kX * kY * kZ);
constexpr float kUntouched = -98765.0F;

// A ceiling that cuts this cube into three pieces of two channels each.
//
// Stated rather than left to the library, because the default budget has a 64 MiB floor and this
// whole cube is 61 KiB: under it, a watched read is one piece and reports once, which is correct
// and shows nothing. A ceiling is the other reason a read splits, and it is the one that fits in a
// test. One chunk row of this cube is eight chunks of 2560 bytes, so a budget of exactly that
// affords one row, and the split axis's chunk is two channels deep.
constexpr std::size_t kThreePieces = 20480;

// The destination is dense in logical order with axis 0 fastest, so the value at (x, y, z) is the
// formula's, at exactly this offset. A transposed read disagrees on the second pixel.
void RequireCubeMatchesTheFormula(const std::vector<float>& destination, const std::string& what) {
    for (std::uint64_t z = 0; z < kZ; ++z) {
        for (std::uint64_t y = 0; y < kY; ++y) {
            for (std::uint64_t x = 0; x < kX; ++x) {
                const auto at = static_cast<std::size_t>((z * kX * kY) + (y * kX) + x);
                Require(destination.at(at) == Value({x, y, z, 0, 0}),
                        what + ": the pixel at " + std::to_string(at) + " is not the one asked for");
            }
        }
    }
}

// The flag is read before the pixels so that a mask the read cannot get leaves the caller's
// destination alone. This is the contract ADR 0005 names as the opposite of the pass's, and the
// reason Image::Read stays outside one -- and it was checked nowhere, because showing it needs a
// source that fails the flag read while the pixel read would have succeeded.
void TestAFailedFlagLeavesTheDestinationAlone() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    source.fail_mask_read(1, ErrorCode::io_error);

    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, ReadOptions{}, ProgressCallback{});

    Require(!read && read.error().code == ErrorCode::io_error, "a failed flag read was not reported");
    Require(source.pixel_reads() == 0, "the pixels were read although the flag could not be");
    for (std::size_t i = 0; i < destination.size(); ++i) {
        Require(destination.at(i) == kUntouched,
                "the destination was written although the flag read failed, at " + std::to_string(i));
    }
}

// A ceiling no amount of splitting gets under is reported rather than allocated past. One plane is
// the smallest piece this request can be cut into, and one byte is less than that.
void TestACeilingTooLowToFitIsRefused() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);

    ReadOptions options;
    options.temporary_memory_limit_bytes = 1;
    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, options, ProgressCallback{});

    Require(!read && read.error().code == ErrorCode::buffer_too_small,
            "a flag buffer over the ceiling was allocated rather than refused");
    Require(source.pixel_reads() == 0, "a refused read still read pixels");
}

// Progress is counted in destination elements, not in chunks, and the finished part is a prefix --
// which is what lets a caller render or forward it as it arrives.
void TestProgressCountsElementsAndFinishesAtTheTotal() {
    const auto image = MakeImage();
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);

    std::vector<std::size_t> written;
    std::size_t reported_total = 0;
    ReadOptions options;
    options.temporary_memory_limit_bytes = kThreePieces;
    const ProgressCallback progress = [&](std::size_t elements_written, std::size_t elements_total) {
        written.push_back(elements_written);
        reported_total = elements_total;
        return true;
    };

    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, options, progress);
    Require(static_cast<bool>(read) && read.value() == kElements, "a watched read did not produce the whole cube");

    Require(reported_total == kElements, "progress reported a total that is not the destination's size");
    Require(source.pixel_reads() == 3, "the ceiling did not cut the cube into three pieces");
    Require(written.size() == 3, "a piece finished without reporting that it had");
    Require(written.back() == kElements, "the last progress report was not the whole read");
    for (std::size_t i = 1; i < written.size(); ++i) {
        Require(written.at(i) > written.at(i - 1), "progress went backwards or stood still");
    }
    RequireCubeMatchesTheFormula(destination, "a watched read");
}

void TestProgressCanStopTheRead() {
    const auto image = MakeImage();
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);

    ReadOptions options;
    options.temporary_memory_limit_bytes = kThreePieces;
    const ProgressCallback progress = [](std::size_t, std::size_t) { return false; };

    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, options, progress);
    Require(!read && read.error().code == ErrorCode::cancelled, "a progress callback returning false did not cancel");
    Require(source.pixel_reads() == 1, "the read carried on past the piece its caller stopped it at");
}

// Splitting is supposed to change how the read is issued and nothing else. The oracle is the
// formula, so this is not two implementations agreeing with each other.
void TestASplitReadAgreesWithAnUnsplitOne() {
    const auto image = MakeImage();
    const auto geometry = MakeGeometry();

    SyntheticPixelSource whole_source(image, geometry, Value);
    std::vector<float> whole(kElements, kUntouched);
    const auto unsplit = ReadInPieces(Readable(whole_source, image, geometry), WholeCube(),
                                      BufferView<float>{whole.data(), whole.size()}, ReadOptions{}, ProgressCallback{});
    Require(static_cast<bool>(unsplit), "the unsplit read failed");
    Require(whole_source.pixel_reads() == 1, "a read with no reason to split was issued in pieces");
    RequireCubeMatchesTheFormula(whole, "an unsplit read");

    SyntheticPixelSource split_source(image, geometry, Value);
    ReadOptions options;
    options.temporary_memory_limit_bytes = kThreePieces;
    std::vector<float> split(kElements, kUntouched);
    const auto in_pieces = ReadInPieces(Readable(split_source, image, geometry), WholeCube(),
                                        BufferView<float>{split.data(), split.size()}, options, ProgressCallback{});
    Require(static_cast<bool>(in_pieces), "the split read failed");
    Require(split_source.pixel_reads() == 3, "the split read was issued in one piece after all");
    Require(split == whole, "a split read and an unsplit one disagreed about the same cube");
}

// What crosses the seam is the caller's buffer from where a piece lands to its end, not the piece's
// size restated. The seam's one check -- that the selection fits what it is writing into -- is only a
// check if the length it is held to comes from the buffer: it used to be worked out from the very
// selection it was compared with, so a piece planned to run past the caller's buffer would have
// been written there. The buffer here is longer than the read, so the rest is visible, and so is
// whether anything was written into it.
void TestEachPieceIsHandedTheRestOfTheBuffer() {
    const auto image = MakeImage();
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);

    constexpr std::size_t kSpare = 7;
    constexpr std::size_t kPiece = static_cast<std::size_t>(kX * kY * 2);
    ReadOptions options;
    options.temporary_memory_limit_bytes = kThreePieces;
    std::vector<float> destination(kElements + kSpare, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, options,
                                   ProgressCallback{});
    Require(static_cast<bool>(read) && read.value() == kElements, "a read into a roomier buffer failed");

    const std::vector<std::size_t> handed{kElements + kSpare, kElements + kSpare - kPiece,
                                          kElements + kSpare - (2 * kPiece)};
    Require(source.pixel_destinations() == handed,
            "a piece was not handed the rest of the caller's buffer from where it lands");
    for (std::size_t i = kElements; i < destination.size(); ++i) {
        Require(destination.at(i) == kUntouched, "a read wrote past what it selected, at " + std::to_string(i));
    }
    destination.resize(kElements);
    RequireCubeMatchesTheFormula(destination, "a read into a roomier buffer");

    // The mask read is one read into the caller's buffer, so it is handed all of it.
    const auto masked = MakeImage(true);
    SyntheticPixelSource flags(masked, geometry, Value);
    std::vector<std::uint8_t> mask(kElements + kSpare, 9);
    const auto mask_read = ReadPixelMask(Readable(flags, masked, geometry), WholeCube(),
                                         BufferView<std::uint8_t>{mask.data(), mask.size()}, ReadControl{});
    Require(static_cast<bool>(mask_read), "a mask read into a roomier buffer failed");
    Require(flags.mask_destinations() == std::vector<std::size_t>{kElements + kSpare},
            "the mask read was not handed the caller's whole buffer");
    for (std::size_t i = kElements; i < mask.size(); ++i) {
        Require(mask.at(i) == 9, "a mask read wrote past what it selected, at " + std::to_string(i));
    }
}

// A flag the read can get is folded in, so a dropped pixel arrives as NaN rather than as a value
// the caller has to know to distrust.
void TestAFlaggedPixelArrivesAsNaN() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    source.set_flags([](const std::vector<std::uint64_t>& logical) { return logical.at(0) % 2 == 0; });

    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, ReadOptions{}, ProgressCallback{});
    Require(static_cast<bool>(read), "a masked read failed");
    Require(source.mask_reads() > 0, "the flag was never read");

    for (std::uint64_t z = 0; z < kZ; ++z) {
        for (std::uint64_t y = 0; y < kY; ++y) {
            for (std::uint64_t x = 0; x < kX; ++x) {
                const auto at = static_cast<std::size_t>((z * kX * kY) + (y * kX) + x);
                const float value = destination.at(at);
                if (x % 2 == 0) {
                    Require(value == Value({x, y, z, 0, 0}), "a good pixel was dropped");
                } else {
                    Require(value != value, "a flagged pixel did not arrive as NaN");
                }
            }
        }
    }
}

// Declining the mask means the flag is never read at all, not that it is read and ignored.
void TestDecliningTheMaskReadsNoFlag() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    source.fail_mask_read(1, ErrorCode::io_error);

    ReadOptions options;
    options.apply_pixel_mask = false;
    std::vector<float> destination(kElements, kUntouched);
    const auto read = ReadInPieces(Readable(source, image, geometry), WholeCube(),
                                   BufferView<float>{destination.data(), destination.size()}, options, ProgressCallback{});
    Require(static_cast<bool>(read), "declining the mask still failed on a flag that cannot be read");
    Require(source.mask_reads() == 0, "declining the mask still read the flag");
    RequireCubeMatchesTheFormula(destination, "a read that declined the mask");
}

// The mask read is the other entry this module has, and until it came here it was written inline in
// the facade and reached the pixel seam directly -- so the only way to its refusals was a fixture on
// disk, and to its cancellation branch there was no way at all.
void TestTheMaskReadIsRefusedForAnImageWithNoFlag() {
    const auto image = MakeImage(false);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);

    std::vector<std::uint8_t> destination(kElements, 7);
    const auto read = ReadPixelMask(Readable(source, image, geometry), WholeCube(),
                                    BufferView<std::uint8_t>{destination.data(), destination.size()},
                                    ReadControl{});
    Require(!read, "an image with no flag has no pixel mask to read");
    Require(read.error().code == ErrorCode::not_found, "and it should say so as not_found");
    Require(source.mask_reads() == 0, "without asking the source for one");
}

// The preamble the two entries share, asked of the one that used to have its own copy: a destination
// too small is refused before any storage work, and so is a request that is already cancelled.
void TestTheMaskReadChecksBeforeItReads() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    const auto readable = Readable(source, image, geometry);

    {
        std::vector<std::uint8_t> too_small(kElements - 1, 0);
        const auto read = ReadPixelMask(readable, WholeCube(),
                                        BufferView<std::uint8_t>{too_small.data(), too_small.size()},
                                        ReadControl{});
        Require(!read && read.error().code == ErrorCode::invalid_argument,
                "a destination one element short should be refused");
        Require(source.mask_reads() == 0, "before the source is asked for anything");
    }

    {
        std::vector<std::uint8_t> destination(kElements, 0);
        ReadControl control;
        control.cancellation_requested = []() { return true; };
        const auto read = ReadPixelMask(readable, WholeCube(),
                                        BufferView<std::uint8_t>{destination.data(), destination.size()},
                                        control);
        Require(!read && read.error().code == ErrorCode::cancelled,
                "a request that is already cancelled should not open an array");
        Require(source.mask_reads() == 0, "and should not reach the source");
    }
}

// What it reads, where it puts it, and what it reports. The flag's own formula is the oracle, and
// the destination is dense in logical order exactly as an ordinary read's is.
void TestTheMaskReadFillsTheDestinationInLogicalOrder() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    // A stripe pattern that depends on all three varying axes, so a transposed read disagrees.
    source.set_flags([](const std::vector<std::uint64_t>& logical) {
        return ((logical.at(0) + (2 * logical.at(1)) + (3 * logical.at(2))) % 5) != 0;
    });

    std::vector<std::uint8_t> destination(kElements, 9);
    const auto read = ReadPixelMask(Readable(source, image, geometry), WholeCube(),
                                    BufferView<std::uint8_t>{destination.data(), destination.size()},
                                    ReadControl{});
    Require(static_cast<bool>(read), "the mask read should succeed");
    Require(read.value() == kElements, "and report every selected element");
    Require(source.mask_reads() == 1, "in a single read, because a mask allocates nothing to bound");
    Require(source.pixel_reads() == 0, "and it should not read the pixels beside it");

    for (std::uint64_t z = 0; z < kZ; ++z) {
        for (std::uint64_t y = 0; y < kY; ++y) {
            for (std::uint64_t x = 0; x < kX; ++x) {
                const auto at = static_cast<std::size_t>((z * kX * kY) + (y * kX) + x);
                const std::uint8_t expected = ((x + (2 * y) + (3 * z)) % 5) != 0 ? 1 : 0;
                Require(destination.at(at) == expected,
                        "the flag at (" + std::to_string(x) + ", " + std::to_string(y) + ", " +
                            std::to_string(z) + ") should be " + std::to_string(expected));
            }
        }
    }
}

// A source that fails the read reports that failure rather than a count, which is the branch the
// facade's copy could only reach by breaking a file underneath it.
void TestAFailedMaskReadIsReported() {
    const auto image = MakeImage(true);
    const auto geometry = MakeGeometry();
    SyntheticPixelSource source(image, geometry, Value);
    source.fail_mask_read(1, ErrorCode::io_error);

    std::vector<std::uint8_t> destination(kElements, 0);
    const auto read = ReadPixelMask(Readable(source, image, geometry), WholeCube(),
                                    BufferView<std::uint8_t>{destination.data(), destination.size()},
                                    ReadControl{});
    Require(!read && read.error().code == ErrorCode::io_error, "a failed flag read should be reported");
}

}  // namespace

int main() {
    try {
        TestAFailedFlagLeavesTheDestinationAlone();
        TestACeilingTooLowToFitIsRefused();
        TestProgressCountsElementsAndFinishesAtTheTotal();
        TestProgressCanStopTheRead();
        TestASplitReadAgreesWithAnUnsplitOne();
        TestEachPieceIsHandedTheRestOfTheBuffer();
        TestAFlaggedPixelArrivesAsNaN();
        TestDecliningTheMaskReadsNoFlag();
        TestTheMaskReadIsRefusedForAnImageWithNoFlag();
        TestTheMaskReadChecksBeforeItReads();
        TestTheMaskReadFillsTheDestinationInLogicalOrder();
        TestAFailedMaskReadIsReported();
        std::cout << "carta-zarr read synthetic tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr read synthetic tests failed: " << error.what() << '\n';
        return 1;
    }
}
