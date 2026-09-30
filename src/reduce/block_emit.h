/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_REDUCE_BLOCK_EMIT_H_
#define CARTA_ZARR_SRC_REDUCE_BLOCK_EMIT_H_

#include "carta-zarr/result.h"

#include "reduce/pass_plan.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>

namespace carta::zarr::internal {

/**
 * One emitted block, while it is being filled.
 *
 * Nothing here is computed by the reduction. It hands `begin` and `end` straight to the walk and
 * indexes by the slab's own relative number, and the absolute one the sink reports comes from the
 * emitter.
 *
 * Which of those numbers is which used to be three paragraphs, here and at both of the types in
 * pass.h, because both were std::uint64_t and nothing else could say it. See ChannelIndex.
 */
struct EmitBlock {
    SelectionChannel begin;
    SelectionChannel end;
    // The chunks of this block the walk has finished, and the reads it has made. Both belong to the
    // block rather than to the reduction: they reset when a block does, which is the difference
    // between a reduction's walk and a whole-plane pass's. A walk made of several footprints carries
    // `reads_done` across them so that a footprint taken in a single read still reports once.
    std::uint64_t chunks_done = 0;
    std::uint64_t reads_done = 0;

    std::uint64_t length() const noexcept {
        return end - begin;
    }
};

/**
 * Cuts a selection into emitted blocks and hands each one over as it fills.
 *
 * ADR 0005 extracted the pass; what it left behind was the loop around the pass, and both
 * reductions wrote that out in full -- cut the block, reset the accumulator, count the block's
 * chunks, hand it over part-filled whenever the walk is about to spend another budget, walk it,
 * hand it over finished, advance. Thirty lines each, with the relative-versus-absolute hazard above
 * intact in both copies, and `pass.h` documenting that hazard in prose because the types could not
 * express it.
 *
 * It also holds an invariant the two used to restate. The chunks one spectral layer occupies is
 * both what the emit budget is spent against and the unit progress is counted in; they are the same
 * number and a reduction that used one for one and one for the other would read the right pixels
 * and report a bar that lies. It is given once, here.
 *
 * What it is not is a reduction base class. The visitor stays a template parameter all the way down
 * -- ADR 0005 measured a quarter of a reduction riding on that -- and so does the walk, because a
 * whole-plane pass and a walk over the chunk runs of a region set are genuinely different walks.
 * This is the envelope; the arithmetic inside it is still the reduction's own.
 */
class BlockEmitter {
public:
    /**
     * `layer_chunks` is what one spectral layer of whatever is being walked occupies: the plan's own
     * for a whole plane, the region set's for a reduction. `bytes_per_channel` and `hint` size the
     * block, as PassPlan::EmitChannels describes.
     */
    // `layer_chunks` is clamped for one of its two uses and not for the other, which is deliberate.
    // Zero is reachable -- a region set whose mask selects nothing occupies no chunks -- and the two
    // want opposite things about it. As the denominator of a part-filled block's completeness it
    // must never be zero, so it is clamped. As the layer the emit budget is spent against, zero is
    // the honest answer and gives the right one: there is nothing to read, so the whole selection is
    // handed over in a single block rather than cut into pieces sized for chunks nobody will decode.
    BlockEmitter(const PassPlan& plan, std::uint64_t layer_chunks, std::size_t bytes_per_channel,
                 std::uint32_t hint,
                 std::string cancelled = "The reduction was cancelled by its sink")
        : _plan(&plan),
          _layer_chunks(std::max<std::uint64_t>(1, layer_chunks)),
          _emit_channels(plan.EmitChannels(layer_chunks, bytes_per_channel, hint)),
          _cancelled(std::move(cancelled)) {}

    BlockEmitter(const BlockEmitter&) = delete;
    BlockEmitter& operator=(const BlockEmitter&) = delete;

    // How many channels one block may hold before chunk alignment has its say. Reported because a
    // test has something to say about it and because nothing else can ask.
    std::uint64_t emit_channels() const noexcept {
        return _emit_channels;
    }

    /**
     * Walk the whole selection, a block at a time.
     *
     * `reset(length)` prepares the accumulator for a block of `length` channels.
     *
     * `walk(block, report)` reads it. It hands `block.begin`, `block.end`, `block.chunks_done` and,
     * when its walk is made of more than one footprint, `block.reads_done` to the pass, and passes
     * `report` as the pass's own before-read callback.
     *
     * `hand_over(first_channel, length, complete, completeness)` fills the reduction's own block
     * struct and calls its sink, returning what the sink returned. It is called once for every read
     * after the first of a block, and once more when the block is finished. `first_channel` is a
     * SelectionChannel, and turning it into the plain number a public block carries is the one place
     * the type is left behind.
     *
     * All three are template parameters and none may become a std::function; `walk` carries the
     * visitor, which is the one ADR 0005 is about.
     */
    template <typename Reset, typename Walk, typename HandOver>
    Result<void> Over(Reset&& reset, Walk&& walk, HandOver&& hand_over) const {
        const SelectionChannel end_of_selection{_plan->planes.spectral.count};

        for (SelectionChannel begin{}; begin < end_of_selection;) {
            EmitBlock block;
            block.begin = begin;
            block.end = _plan->AlignedSlabEnd(begin, _emit_channels, end_of_selection);

            reset(block.length());

            // What the completeness of a part-filled block is a fraction of.
            const std::uint64_t chunks_total = _plan->ChunksCovering(_layer_chunks, block.begin, block.end);

            const auto deliver = [&](bool complete) -> Result<void> {
                const double completeness =
                    complete ? 1.0
                             : static_cast<double>(block.chunks_done) / static_cast<double>(chunks_total);
                if (!hand_over(block.begin, block.length(), complete, completeness)) {
                    return Error{ErrorCode::cancelled, _cancelled, _plan->descriptor->id};
                }
                return {};
            };

            if (auto walked = walk(block, [&](std::uint64_t) -> Result<void> { return deliver(false); });
                !walked) {
                return walked.error();
            }
            if (auto handed = deliver(true); !handed) {
                return handed.error();
            }
            begin = block.end;
        }
        return {};
    }

private:
    const PassPlan* _plan;
    std::uint64_t _layer_chunks;
    std::uint64_t _emit_channels;
    std::string _cancelled;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_BLOCK_EMIT_H_
