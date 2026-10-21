/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// What an ordinary read decides before it reads anything: whether to cut itself into pieces, where,
// and how big each one should be.
//
// None of this needs a store, a transport or a fixture, so the decisions are checked directly
// rather than inferred from the progress callbacks of a read.

#include "chunk_blocks.h"
#include "read/pieces.h"
#include "support/check.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

// What a context gives every read here to aim its default budget at.
constexpr std::size_t kDecodeThreads = 4;

using carta::zarr::AxisRole;
using carta::zarr::ChunkGeometry;
using carta::zarr::ImageDescriptor;
using carta::zarr::Range;
using carta::zarr::ReadOptions;
using carta::zarr::ReadRequest;
using carta::zarr::internal::kHeldBytesPerDecodedByte;
using carta::zarr::internal::Piece;
using carta::zarr::internal::PlanPieces;

using carta::zarr::testing::Require;

// An image in the logical order XRADIO reports: x, y, spectral, polarization, time.
ImageDescriptor MakeImage(std::uint64_t x, std::uint64_t y, std::uint64_t channels, bool has_mask = false) {
    ImageDescriptor descriptor;
    descriptor.id = "SKY";
    descriptor.stored_type = carta::zarr::DataType::float32;
    descriptor.pixel_mask_id = has_mask ? "FLAG" : "";
    const struct {
        const char* name;
        AxisRole role;
        std::uint64_t length;
    } axes[]{{"l", AxisRole::spatial_x, x},
             {"m", AxisRole::spatial_y, y},
             {"frequency", AxisRole::spectral, channels},
             {"polarization", AxisRole::polarization, 1},
             {"time", AxisRole::time, 1}};
    for (std::size_t i = 0; i < 5; ++i) {
        carta::zarr::AxisDescriptor axis;
        axis.name = axes[i].name;
        axis.role = axes[i].role;
        axis.length = axes[i].length;
        axis.storage_index = std::vector<std::size_t>{3, 4, 1, 2, 0}.at(i);
        descriptor.axes.push_back(axis);
    }
    return descriptor;
}

ChunkGeometry MakeGeometry(std::uint64_t chunk_x, std::uint64_t chunk_y, std::uint64_t chunk_z) {
    ChunkGeometry geometry;
    geometry.chunk_shape = {chunk_x, chunk_y, chunk_z, 1, 1};
    return geometry;
}

ReadRequest WholeImage(const ImageDescriptor& descriptor) {
    ReadRequest request;
    for (const auto& axis : descriptor.axes) {
        request.axes.push_back(Range{0, axis.length, 1});
    }
    return request;
}

std::uint64_t ElementCount(const ReadRequest& request) {
    std::uint64_t elements = 1;
    for (const auto& range : request.axes) {
        elements *= range.count;
    }
    return elements;
}

// A read with a progress callback, which is cut to the chunks its budget affords; most of what is
// pinned here is how. See TestAReadWithoutProgressIsCutToItsBudgetButAPlane for a read without one.
std::vector<Piece> Pieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                          const ChunkGeometry& flag_geometry, const ReadRequest& request, const ReadOptions& options) {
    return PlanPieces(descriptor, geometry, flag_geometry, request, options, kDecodeThreads, true);
}

// A flag chunked as its pixels are, which is every image here but one.
std::vector<Piece> Pieces(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const ReadRequest& request,
                          const ReadOptions& options) {
    return Pieces(descriptor, geometry, geometry, request, options);
}

// A budget that affords `chunks` unmasked float32 chunks of `elements` each, whose pixels the caller
// holds: what each of them holds is what it decodes to, kHeldBytesPerDecodedByte over.
std::size_t BudgetFor(std::uint64_t chunks, std::uint64_t elements) {
    return static_cast<std::size_t>(chunks * elements * sizeof(float) * kHeldBytesPerDecodedByte);
}

bool SameRange(const Range& a, const Range& b) {
    return a.start == b.start && a.count == b.count && a.stride == b.stride;
}

// The axes on which a piece's request differs from the whole read's.
std::vector<std::size_t> CutAxes(const Piece& piece, const ReadRequest& whole) {
    std::vector<std::size_t> cut;
    for (std::size_t i = 0; i < whole.axes.size(); ++i) {
        if (!SameRange(piece.request.axes.at(i), whole.axes.at(i))) {
            cut.push_back(i);
        }
    }
    return cut;
}

// The pieces fill the destination end to end, each starting where the one before it stopped and the
// last stopping at the end -- which is what makes the finished part a prefix, and what a progress
// report counts.
void RequireTheyFillTheDestination(const std::vector<Piece>& pieces, const ReadRequest& whole,
                                   const std::string& what) {
    std::uint64_t filled = 0;
    for (const auto& piece : pieces) {
        Require(piece.first_element == filled, what + ": a piece did not start where the one before it stopped");
        filled += ElementCount(piece.request);
    }
    Require(filled == ElementCount(whole), what + ": the pieces did not fill the destination");
}

// A read that fits the library's own budget is issued exactly as it was asked for. The plan still
// describes it, as one piece covering everything, so the loop that reads it has one shape rather than
// two.
void TestAReadThatFitsIsOnePiece() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    const auto request = WholeImage(image);

    const auto pieces = Pieces(image, geometry, request, {});
    Require(pieces.size() == 1, "a read that fits its budget was cut");
    Require(CutAxes(pieces.front(), request).empty() && pieces.front().first_element == 0 &&
                pieces.front().segments.empty(),
            "an uncut read should be one piece covering everything");
}

// A read with a progress callback is cut to its budget. Without a stated budget it is the library's
// own, so it takes a read larger than that to show: 8192 planes of 64 x 64 float32 decode to 128 MiB
// and hold three times that. Planning allocates nothing, so the size costs nothing here.
void TestAReadWithProgressIsCutToItsBudget() {
    const auto large = MakeImage(64, 64, 8192);
    const auto large_geometry = MakeGeometry(32, 64, 1);
    const auto large_request = WholeImage(large);
    const auto unbounded = Pieces(large, large_geometry, large_request, ReadOptions{});
    Require(unbounded.size() > 1, "a read larger than the library's own budget was not cut");
    RequireTheyFillTheDestination(unbounded, large_request, "a cut to the library's budget");

    const auto image = MakeImage(64, 64, 8);
    ReadOptions bounded;
    bounded.read_budget_bytes = 4096;
    Require(Pieces(image, MakeGeometry(32, 64, 1), WholeImage(image), bounded).size() > 1,
            "a read with a memory ceiling was not cut");
}

// A read with a progress callback is cut to no more than kWatchedPieceChunks a piece however much
// its budget affords, so that a spectrum reports as often on many decode threads as on few -- but
// to no less than one chunk along the cut, which a plane of many chunks is, so that it is never
// gathered in segments.
void TestAPieceWithProgressIsSixteenChunksOrOneAlongTheCut() {
    constexpr std::size_t kManyThreads = 28;
    ReadOptions plenty;
    plenty.read_budget_bytes = std::size_t{2} << 30;

    // A spectrum of one pixel through 1024 channels of 16 x 16 x 1 chunks is a chunk a channel.
    const auto image = MakeImage(64, 64, 1024);
    const auto geometry = MakeGeometry(16, 16, 1);
    auto spectrum = WholeImage(image);
    spectrum.axes.at(0) = Range{40, 1, 1};
    spectrum.axes.at(1) = Range{20, 1, 1};
    const auto pieces = PlanPieces(image, geometry, geometry, spectrum, plenty, kManyThreads, true);
    Require(pieces.size() == 1024 / carta::zarr::internal::kWatchedPieceChunks,
            "a spectrum with progress was not cut into pieces of sixteen chunks");
    for (const auto& piece : pieces) {
        Require(piece.segments.empty(), "a piece of a spectrum with progress was read in segments");
    }
    RequireTheyFillTheDestination(pieces, spectrum, "a spectrum with progress");
    Require(PlanPieces(image, geometry, geometry, spectrum, plenty, kManyThreads, false).size() == 1,
            "a spectrum without progress, which its budget affords whole, was cut");

    // A plane of 8 x 8 x 1 chunks is sixty-four of them, more than sixteen: a piece a plane, whole.
    const auto planes = MakeImage(64, 64, 8);
    const auto fine = MakeGeometry(8, 8, 1);
    const auto plane_pieces = PlanPieces(planes, fine, fine, WholeImage(planes), plenty, kManyThreads, true);
    Require(plane_pieces.size() == 8, "a read with progress of planes of many chunks was not cut a plane a piece");
    for (const auto& piece : plane_pieces) {
        Require(piece.segments.empty(), "a plane of a read with progress was read in segments");
    }
}

// A read with no progress callback is cut to its budget as one with a callback is, into segments where
// a piece one chunk deep is still too much -- except a single plane, which is read whole once its
// budget affords a chunk a decode thread: cut, it reads slower one segment after another, and what it
// holds beyond its budget is bounded by the plane.
void TestAReadWithoutProgressIsCutToItsBudgetButAPlane() {
    const auto without_progress = [](const ImageDescriptor& image, const ChunkGeometry& geometry,
                                     const ReadRequest& request, const ReadOptions& options, std::size_t threads) {
        return PlanPieces(image, geometry, geometry, request, options, threads, false);
    };
    const auto no_segments = [](const std::vector<Piece>& pieces) {
        for (const auto& piece : pieces) {
            if (!piece.segments.empty()) {
                return false;
            }
        }
        return true;
    };

    // Chunks of 16 x 16 float32 decode to 1 KiB and hold 3 KiB; a plane of 64 x 64 is sixteen of them.
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(16, 16, 1);
    ReadOptions four;
    four.read_budget_bytes = 4 * 3 * 1024;
    const auto whole = without_progress(image, geometry, WholeImage(image), four, 4);
    Require(whole.size() == 8, "a read without progress was not cut a plane a piece, one chunk along the cut");
    RequireTheyFillTheDestination(whole, WholeImage(image), "a read without progress cut to its budget");
    for (const auto& piece : whole) {
        Require(piece.segments.size() == 4, "a plane of a read without progress of many was not read in four segments");
    }

    auto plane = WholeImage(image);
    plane.axes.at(2) = Range{3, 1, 1};
    const auto one = without_progress(image, geometry, plane, four, 4);
    Require(one.size() == 1 && no_segments(one),
            "a single plane without progress was cut, though its budget affords a chunk a thread");
    // Cut along m, a row of four chunks a piece, which the budget affords.
    Require(without_progress(image, geometry, plane, four, 5).size() == 4,
            "a single plane whose budget affords fewer chunks than threads was not cut to it");
    Require(PlanPieces(image, geometry, geometry, plane, four, 4, true).size() == 4,
            "a single plane with progress was not cut to its budget");

    ReadOptions all;
    all.read_budget_bytes = 128 * 3 * 1024;
    Require(without_progress(image, geometry, WholeImage(image), all, 4).size() == 1,
            "a read without progress its budget affords whole was cut");
    ReadOptions half;
    half.read_budget_bytes = 64 * 3 * 1024;
    const auto halves = without_progress(image, geometry, WholeImage(image), half, 4);
    Require(halves.size() == 2 && no_segments(halves),
            "a read without progress was not cut to what its budget affords");
}

// The cut goes on the slowest-varying axis that selects more than one element, because the
// destination is dense with axis 0 fastest and that is the only axis whose pieces extend a prefix.
void TestTheCutGoesOnTheSlowestSelectedAxis() {
    const auto image = MakeImage(64, 64, 8);
    // Chunks of 32 x 16, so one plane is eight of them, 16 KiB, and one row of chunks along m is two.
    const auto geometry = MakeGeometry(32, 16, 1);
    const auto request = WholeImage(image);
    // Axes 3 and 4 are degenerate here, so the spectrum at index 2 is the slowest one selected. A
    // ceiling of one plane makes each channel a piece.
    ReadOptions a_plane;
    a_plane.read_budget_bytes = BudgetFor(8, 32 * 16);
    const auto pieces = Pieces(image, geometry, request, a_plane);
    Require(pieces.size() == 8, "a ceiling of one plane did not make a piece of each channel");
    for (const auto& piece : pieces) {
        Require(CutAxes(piece, request) == std::vector<std::size_t>{2},
                "the cut should have gone on the spectral axis and nowhere else");
        // One channel is worth a plane of the destination.
        Require(piece.first_element == piece.request.axes.at(2).start * 64 * 64,
                "a piece of the spectrum did not land a plane per channel in");
    }
    RequireTheyFillTheDestination(pieces, request, "a spectrum cut");

    // Narrow the spectrum to one channel and the only axis left with more than one element is m. A
    // ceiling of one row of chunks makes each sixteen rows a piece.
    auto one_channel = request;
    one_channel.axes.at(2) = Range{0, 1, 1};
    ReadOptions a_row_of_chunks;
    a_row_of_chunks.read_budget_bytes = BudgetFor(2, 32 * 16);
    const auto narrowed = Pieces(image, geometry, one_channel, a_row_of_chunks);
    Require(narrowed.size() == 4, "a ceiling of one row of chunks did not cut the plane into four");
    for (const auto& piece : narrowed) {
        Require(CutAxes(piece, one_channel) == std::vector<std::size_t>{1},
                "with one channel the cut should move to m and nowhere else");
        Require(piece.first_element == piece.request.axes.at(1).start * 64, "one m is worth a row of the destination");
    }
    RequireTheyFillTheDestination(narrowed, one_channel, "a plane cut");
}

// A request that selects a single element of every axis has nowhere to be cut, so it is one piece
// whatever the caller asked for, holding the one chunk it cannot avoid.
void TestAReadWithNowhereToCutIsOnePiece() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(32, 64, 1);
    ReadRequest single;
    single.axes.assign(5, Range{0, 1, 1});

    ReadOptions bounded;
    bounded.read_budget_bytes = 1;
    Require(Pieces(image, geometry, single, bounded).size() == 1, "a single-element read has nowhere to be cut");
}

// A piece is measured in chunks, not in elements: the budget buys whole chunks along the cut axis
// because asking for part of one decodes all of it anyway -- so a piece ends on a chunk boundary.
void TestATighterCeilingBuysFewerChunks() {
    const auto image = MakeImage(64, 64, 32);
    const auto geometry = MakeGeometry(64, 64, 4);
    const auto request = WholeImage(image);

    ReadOptions generous;
    generous.read_budget_bytes = BudgetFor(8, 64 * 64 * 4);
    ReadOptions tight;
    tight.read_budget_bytes = BudgetFor(1, 64 * 64 * 4);

    const auto wide = Pieces(image, geometry, request, generous);
    const auto narrow = Pieces(image, geometry, request, tight);
    Require(narrow.size() > wide.size(), "a tighter ceiling should buy smaller pieces");
    for (const auto* pieces : {&wide, &narrow}) {
        for (const auto& piece : *pieces) {
            const auto& range = piece.request.axes.at(2);
            Require(range.start % 4 == 0 && (range.start + range.count == 32 || range.count % 4 == 0),
                    "a piece did not begin and end on a chunk boundary of the spectrum");
        }
        RequireTheyFillTheDestination(*pieces, request, "a budgeted cut");
    }
}

// Where a piece ends is the arithmetic that most needed moving: the end is rounded to a chunk
// boundary so that no chunk is decoded by two pieces -- the last element one piece reads and the
// first the next one reads never share a chunk.
//
// It takes a request that does not begin on a chunk boundary to show. Starting at zero, a piece a
// whole number of chunks long ends on a boundary without being moved there, and cutting at the naive
// end passes every other test in this repository.
void TestNoChunkIsReadByTwoPieces() {
    const auto image = MakeImage(64, 64, 32);
    const auto geometry = MakeGeometry(64, 64, 4);
    ReadOptions tight;
    tight.read_budget_bytes = BudgetFor(1, 64 * 64 * 4);

    const auto check = [&](const Range& spectrum, const std::string& what) {
        auto request = WholeImage(image);
        request.axes.at(2) = spectrum;
        const auto pieces = Pieces(image, geometry, request, tight);
        Require(pieces.size() > 1, what + " under a one-chunk budget was not cut");
        for (std::size_t i = 0; i + 1 < pieces.size(); ++i) {
            const auto& here = pieces.at(i).request.axes.at(2);
            const auto& next = pieces.at(i + 1).request.axes.at(2);
            Require(here.stride == spectrum.stride && next.stride == spectrum.stride,
                    what + ": a piece lost the request's stride");
            const auto last = here.start + ((here.count - 1) * here.stride);
            Require(last / 4 < next.start / 4, what + ": two pieces read from one chunk of the spectrum");
            Require(next.start == last + here.stride, what + ": a piece skipped or repeated a selected channel");
        }
        RequireTheyFillTheDestination(pieces, request, what);
    };

    check(Range{2, 28, 1}, "a read starting mid-chunk");
    check(Range{1, 15, 2}, "a strided read");  // channels 1, 3, ..., 29
}

// The flag is decoded beside the pixels, so a read that will apply it costs more per chunk and the
// same ceiling has to buy less. Counting it on one side of the division and not the other would
// size pieces against a cost the read does not have.
void TestApplyingTheFlagCostsTheBudget() {
    const auto image = MakeImage(64, 64, 32, /*has_mask=*/true);
    const auto geometry = MakeGeometry(64, 64, 1);
    const auto request = WholeImage(image);

    ReadOptions declined;
    declined.read_budget_bytes = 64 * 64 * sizeof(float) * 6;
    declined.apply_pixel_mask = false;
    ReadOptions applied = declined;
    applied.apply_pixel_mask = true;

    const auto without = Pieces(image, geometry, request, declined);
    const auto with = Pieces(image, geometry, request, applied);
    Require(with.size() > without.size(),
            "a read that also decodes the flag cannot afford as much of the spectrum per piece");
}

// A flag chunked coarser than its pixels decodes the whole of a flag chunk beside each pixel chunk,
// and a piece is sized against that. Two pixel chunks of 2 x 5 float32 cover the row, and the flag is
// one chunk of 4 x 5 x 2. A pixel chunk then decodes 40 bytes and the 40 of flag beside it, and holds
// three times that and its 10-byte share of the folded-in flag: 250, so a 400-byte budget buys one
// of them. Counting the flag in the pixels' chunks, 10 bytes each, it holds 160 and buys both.
void TestACoarseFlagCostsItsOwnChunks() {
    const auto image = MakeImage(4, 5, 2, /*has_mask=*/true);
    const auto geometry = MakeGeometry(2, 5, 1);
    ChunkGeometry flag;
    flag.chunk_shape = {4, 5, 2, 1, 1};
    auto request = WholeImage(image);
    request.axes.at(1) = Range{0, 1, 1};
    request.axes.at(2) = Range{0, 1, 1};

    ReadOptions budget;
    budget.read_budget_bytes = 400;
    Require(Pieces(image, geometry, geometry, request, budget).size() == 1,
            "a flag chunked alike costs a byte an element, so both pixel chunks fit");
    const auto pieces = Pieces(image, geometry, flag, request, budget);
    Require(pieces.size() == 2, "a piece was sized as if the flag beside it were chunked as its pixels");
    RequireTheyFillTheDestination(pieces, request, "a coarse flag");

    budget.apply_pixel_mask = false;
    Require(Pieces(image, geometry, flag, request, budget).size() == 1, "an unmasked read paid for a flag");
}

// A stride of a chunk or more along an axis the read is not cut on steps over whole chunks, and the
// chunks it steps over are not decoded -- so they are not what a piece's budget is spent on. Every
// other column of 16-wide chunks is two chunks across the image, not the three its first and last
// element span, and a budget of 16 chunks buys two channels a piece rather than one. (A piece of a
// read with a progress callback is no more than kWatchedPieceChunks anyway, so the budget stays
// within it.)
void TestChunksAStrideStepsOverCostNothing() {
    const auto image = MakeImage(64, 64, 32);
    const auto geometry = MakeGeometry(16, 16, 1);
    auto request = WholeImage(image);
    request.axes.at(0) = Range{0, 2, 32};  // columns 0 and 32: chunks 0 and 2
    ReadOptions budget;
    budget.read_budget_bytes = BudgetFor(16, 16 * 16);

    const auto pieces = Pieces(image, geometry, request, budget);
    for (std::size_t i = 0; i + 1 < pieces.size(); ++i) {
        Require(pieces.at(i).request.axes.at(2).count == 2,
                "a piece holds " + std::to_string(pieces.at(i).request.axes.at(2).count) +
                    " channels: two chunks across and four down is eight a channel, and 16 buys two");
    }
    RequireTheyFillTheDestination(pieces, request, "a read striding over whole chunks");
}

// The same along the axis the read is cut on: every eighth channel of a spectrum chunked four deep is
// one chunk per channel, so a budget of two chunks is two channels per piece, not the one that
// dividing the chunk depth by the stride would give.
void TestAStrideOverWholeChunksAlongTheCutIsAChunkAnElement() {
    const auto image = MakeImage(64, 64, 32);
    const auto geometry = MakeGeometry(64, 64, 4);
    auto request = WholeImage(image);
    request.axes.at(2) = Range{0, 4, 8};  // channels 0, 8, 16 and 24, a chunk each
    ReadOptions budget;
    budget.read_budget_bytes = BudgetFor(2, 64 * 64 * 4);

    const auto pieces = Pieces(image, geometry, request, budget);
    Require(pieces.size() == 2, "four channels a chunk apart under a two-chunk budget are two pieces, not " +
                                    std::to_string(pieces.size()));
    RequireTheyFillTheDestination(pieces, request, "a strided spectrum");
}

// The least a piece can be is one chunk along the cut and the whole of every other axis. A plane of
// chunks deep along the spectrum -- the 512 x 512 x 64 layout that held gigabytes for one plane -- has
// rows of chunks that hold more than a read affords, so each row is read in segments along the next
// axis, as many chunks as the budget affords, each straight into its place. The rows are still the pieces, so the
// finished part is still a prefix.
void TestAPieceTooLargeIsReadInSegments() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(16, 16, 8);
    auto request = WholeImage(image);
    request.axes.at(2) = Range{3, 1, 1};

    // Two chunks whose pixels the caller holds: less than a row of four. A segment is read straight into
    // its place, so the caller holds its pixels too, and it is two chunks.
    ReadOptions budget;
    budget.read_budget_bytes = BudgetFor(2, 16 * 16 * 8);
    const auto pieces = Pieces(image, geometry, request, budget);
    Require(pieces.size() == 4, "each row of chunks should be a piece, not " + std::to_string(pieces.size()));
    RequireTheyFillTheDestination(pieces, request, "a plane read in segments");
    for (const auto& piece : pieces) {
        Require(piece.segments.size() == 2, "a row of four chunks affording two at a time is two segments, not " +
                                                std::to_string(piece.segments.size()));
        std::uint64_t covered = 0;
        std::uint64_t next_x = 0;
        for (const auto& segment : piece.segments) {
            Require(SameRange(segment.axes.at(1), piece.request.axes.at(1)) &&
                        SameRange(segment.axes.at(2), piece.request.axes.at(2)),
                    "a segment should cut only the axis below the piece's");
            Require(segment.axes.at(0).start == next_x && segment.axes.at(0).count == 32,
                    "segments should cover the row two chunks at a time, in order");
            next_x += segment.axes.at(0).count;
            covered += ElementCount(segment);
        }
        Require(covered == ElementCount(piece.request), "the segments of a piece should cover it exactly");
    }

    // A budget that affords the row reads it whole.
    budget.read_budget_bytes = BudgetFor(4, 16 * 16 * 8);
    for (const auto& piece : Pieces(image, geometry, request, budget)) {
        Require(piece.segments.empty(), "a piece that fits was read in segments");
    }
}

// Segments go as deep as they must: a request over several channels of chunks that are a whole plane
// across has nothing to cut but the axes below the spectrum, and when a row along m is still too much
// they go on to x, down to a single chunk.
void TestSegmentsGoDownToOneChunk() {
    const auto image = MakeImage(64, 64, 8);
    const auto geometry = MakeGeometry(16, 16, 8);
    const auto request = WholeImage(image);

    ReadOptions budget;
    budget.read_budget_bytes = 1;
    const auto pieces = Pieces(image, geometry, request, budget);
    Require(pieces.size() == 1, "the spectrum is one chunk deep, so there is one piece along it");
    Require(pieces.front().segments.size() == 16, "sixteen chunks a single chunk at a time is sixteen segments, not " +
                                                      std::to_string(pieces.front().segments.size()));
    for (const auto& segment : pieces.front().segments) {
        Require(segment.axes.at(0).count == 16 && segment.axes.at(1).count == 16,
                "a segment past every cut is a single chunk");
    }
}

// A segment cut from a run that starts inside a chunk can reach the run's end and still span two of
// them, which cutting it the same way again would never undo: it is cut a chunk at a time instead.
void TestASegmentStartingInsideAChunkIsCutAtItsBoundary() {
    const auto image = MakeImage(16, 16, 8);
    const auto geometry = MakeGeometry(16, 16, 4);
    auto request = WholeImage(image);
    request.axes.at(0) = Range{5, 1, 1};
    request.axes.at(1) = Range{7, 1, 1};
    request.axes.at(2) = Range{1, 4, 1};

    ReadOptions budget;
    budget.read_budget_bytes = 1;
    const auto pieces = Pieces(image, geometry, request, budget);
    RequireTheyFillTheDestination(pieces, request, "a column across a chunk boundary");
    std::vector<Range> parts;
    for (const auto& piece : pieces) {
        if (piece.segments.empty()) {
            parts.push_back(piece.request.axes.at(2));
        }
        for (const auto& segment : piece.segments) {
            parts.push_back(segment.axes.at(2));
        }
    }
    Require(parts.size() == 2 && SameRange(parts.at(0), Range{1, 3, 1}) && SameRange(parts.at(1), Range{4, 1, 1}),
            "channels 1 to 4 over chunks of four should be read as 1 to 3 and then 4");
}

}  // namespace

int main() {
    try {
        TestAReadThatFitsIsOnePiece();
        TestAReadWithProgressIsCutToItsBudget();
        TestAPieceWithProgressIsSixteenChunksOrOneAlongTheCut();
        TestAReadWithoutProgressIsCutToItsBudgetButAPlane();
        TestTheCutGoesOnTheSlowestSelectedAxis();
        TestAReadWithNowhereToCutIsOnePiece();
        TestATighterCeilingBuysFewerChunks();
        TestNoChunkIsReadByTwoPieces();
        TestApplyingTheFlagCostsTheBudget();
        TestACoarseFlagCostsItsOwnChunks();
        TestChunksAStrideStepsOverCostNothing();
        TestAStrideOverWholeChunksAlongTheCutIsAChunkAnElement();
        TestAPieceTooLargeIsReadInSegments();
        TestSegmentsGoDownToOneChunk();
        TestASegmentStartingInsideAChunkIsCutAtItsBoundary();
    } catch (const std::exception& error) {
        std::cerr << "piece plan test failed: " << error.what() << "\n";
        return 1;
    }
    std::cout << "carta-zarr piece plan tests passed\n";
    return 0;
}
