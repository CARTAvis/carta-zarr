/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_SRC_REDUCE_PASS_PLAN_H_
#define CARTA_ZARR_SRC_REDUCE_PASS_PLAN_H_

// What a pass decides before it reads anything, apart from the walk that carries it out.
//
// The plan is pure and the walk is not, so they are separate headers: code that needs to know how a
// pass is cut does not have to include the slab reader, the walk's templates and the pixel source
// seam. The walk is in pass.h, which includes this.

#include "carta-zarr/read.h"

#include "axis_map.h"
#include "chunk_blocks.h"
#include "reduce/plane_selection.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace carta::zarr::internal {

/**
 * An index into a run of channels, tagged by which run.
 *
 * Indices into different runs do not convert to each other, which is why this exists. Passing the
 * slab reader an index into the block being filled, where it takes an index into the pass's whole
 * spectral selection, would make every block after the first re-read the first block's channels,
 * and the totals would hide it: each block would still read the right number of pixels.
 *
 * An index plus a count is an index into the same run; two indices into one run are a count apart.
 * Writing the wrong one on purpose is still possible (`SelectionChannel{relative}` compiles), but it
 * cannot happen by plain assignment.
 */
template <typename Tag>
struct ChannelIndex {
    std::uint64_t index = 0;

    friend bool operator<(ChannelIndex a, ChannelIndex b) noexcept { return a.index < b.index; }
    friend bool operator==(ChannelIndex a, ChannelIndex b) noexcept { return a.index == b.index; }
    friend ChannelIndex operator+(ChannelIndex a, std::uint64_t channels) noexcept {
        return ChannelIndex{a.index + channels};
    }
    friend std::uint64_t operator-(ChannelIndex a, ChannelIndex b) noexcept { return a.index - b.index; }
};

// An index into the pass's own spectral selection: channel `index` of the image is
// `planes.spectral.start + index * planes.spectral.stride`.
using SelectionChannel = ChannelIndex<struct SelectionChannelTag>;

// An index into the channel range one walk of a pass covers (a block, or the whole run), which a
// visitor accumulating into its own block indexes by.
using BlockChannel = ChannelIndex<struct BlockChannelTag>;

/**
 * Everything a pass decides before a byte is read: which axis is contiguous, how wide a band is,
 * how deep a slab goes, and how much a read may decode.
 *
 * Pure, and separate from the pass for the same reason `PlanRowTasks` is separate from the pool: it
 * is the half with an answer worth testing, and testing it needs no store, no transport and no
 * fixture.
 */
class PassPlan;

PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry, const ChunkGeometry& flag_geometry,
                  const AxisMap& map, const CheckedPlanes& planes, std::uint64_t sample, const ReadOptions& options,
                  std::size_t decode_threads);

class PassPlan {
public:
    const ImageDescriptor* descriptor = nullptr;
    AxisMap map;
    // u is the spatial axis the store varies fastest; v is the other one.
    std::size_t axis_u = 0;
    std::size_t axis_v = 0;
    std::uint64_t u_length = 0;
    std::uint64_t v_length = 0;
    std::uint64_t chunk_u = 1;
    std::uint64_t chunk_v = 1;
    // What one chunk decodes to, counting the flag chunks with it when this read applies the mask (in
    // the flag's own chunks; see DecodedFlagBytes); what it holds while it is read, slab and mask
    // buffers included; and how much a single read may hold. ADR 0005 depends on the first and
    // ADR 0021 on the second, and the third is the caller's own limit when it gave one.
    std::uint64_t chunk_bytes = 1;
    std::uint64_t held_bytes = 1;
    std::size_t slab_budget_bytes = 0;
    std::uint64_t band_rows = 1;
    // The chunks of one spectral layer of the plane that the sample has a pixel in: every chunk,
    // unless the sample steps over some. A whole-plane walk decodes these for each layer, so its
    // reads are sized by them and its progress counted in them.
    std::uint64_t layer_chunks = 1;
    // The planes this pass is over, already checked against the descriptor above.
    PlaneSelection planes;
    // Take every nth pixel along both spatial axes. This does not reduce the chunks a read decodes,
    // since a chunk is decoded whole however few of its pixels are wanted, so it helps only when it
    // steps over whole chunks.
    std::uint64_t sample = 1;
    bool apply_mask = false;

    // Whether u is the image's y rather than its x, so a caller holding coordinates in the image's
    // own axes knows whether to swap them into the walk's. A reduction asks this rather than working
    // it out again from the geometry, so the rule lives only where the plan picks axis_u.
    bool SwapsSpatial() const noexcept { return axis_u == map.y; }

    // How many chunks along the spectrum the selected channels [begin, end) touch.
    //
    // A walk counts its progress in chunks, and on the spectral axis a selection's stride makes that
    // more than a division. It is computed from the run rather than its length: a run that starts
    // part-way into one chunk and ends part-way into another touches one more chunk than its length
    // divides into. The walk's slabs end on chunk boundaries and each counts what it touched, so a
    // total taken from the length alone would be too small, and progress would pass one before the
    // last read.
    std::uint64_t ChunksTouched(SelectionChannel begin, SelectionChannel end) const {
        return ::carta::zarr::internal::ChunksTouched(planes.spectral.start + (begin.index * planes.spectral.stride),
                                                      end - begin, planes.spectral.stride, _chunk_depth);
    }

    // The chunks a walk over the channels [begin, end) decodes, when each spectral layer of what it
    // walks occupies `layer_chunks`: what its progress is a fraction of. Never zero, so it can be
    // divided by: a selection holds at least one channel, and one that occupies nothing still reports
    // against a whole.
    //
    // A block's completeness and a cube histogram's progress are both this number.
    std::uint64_t ChunksCovering(std::uint64_t layer_chunks, SelectionChannel begin, SelectionChannel end) const {
        return std::max<std::uint64_t>(1, layer_chunks * ChunksTouched(begin, end));
    }

    // How many chunks one read may decode: the budget in the units everything else here counts in.
    // Zero when the budget is smaller than a single chunk, which UnitsAffordable floors at one.
    //
    // Public because it is all of the budget the occupancy needs, and passing one number rather than
    // a PassPlan keeps the occupancy testable with nothing linked behind it (ADR 0006).
    std::uint64_t ChunksPerRead() const noexcept { return slab_budget_bytes / std::max<std::uint64_t>(1, held_bytes); }

    // How many units of `chunks_per_unit` chunks one read's budget allows: the rule in chunk_blocks.h,
    // with this plan's budget.
    //
    // A unit is at least one chunk. The one caller that can mean none, a region set covering nothing,
    // handles that itself, in EmitChannels.
    std::uint64_t UnitsAffordable(std::uint64_t chunks_per_unit) const {
        return ::carta::zarr::internal::UnitsAffordable(ChunksPerRead(), chunks_per_unit);
    }

    // How many channels one slab may hold when its spatial footprint occupies `footprint_chunks`
    // chunks of each spectral chunk it touches.
    std::uint64_t SlabChannels(std::uint64_t footprint_chunks) const {
        return std::max<std::uint64_t>(1, UnitsAffordable(footprint_chunks) * _least_channels);
    }

    // The end of a slab that begins at `begin` and would like to be `desired` channels long, moved
    // onto a chunk boundary so that no decode serves two slabs.
    SelectionChannel AlignedSlabEnd(SelectionChannel begin, std::uint64_t desired, SelectionChannel end) const {
        return SelectionChannel{AlignedBlockEnd(begin.index, desired, end.index, planes.spectral.start,
                                                planes.spectral.stride, _chunk_depth)};
    }

    // How many channels one emitted block may hold.
    //
    // The hint is the caller's, the budgets are the library's, and the chunk alignment is the plan's;
    // the smallest wins, and the block reports what it used. Without a hint, a block costs one budget
    // of decoded bytes, as a piece of Read does, and spends whatever the spatial walk left over. A
    // small region leaves almost all of it, and the block spans many chunks along the spectrum; a
    // region covering the image spends the budget spatially, and the block is the single chunk layer
    // the walk is already reading.
    //
    // `layer_chunks` is the chunks one spectral layer of whatever the caller is walking occupies,
    // which is the plan's own for a whole plane and the region set's for a reduction.
    std::uint64_t EmitChannels(std::uint64_t layer_chunks, std::size_t bytes_per_channel, std::uint32_t hint) const;

private:
    friend PassPlan PlanPass(const ImageDescriptor& descriptor, const ChunkGeometry& geometry,
                             const ChunkGeometry& flag_geometry, const AxisMap& map, const CheckedPlanes& planes,
                             std::uint64_t sample, const ReadOptions& options, std::size_t decode_threads);

    // Intermediate values for the answers above, kept here so that no caller works them out itself.
    // Tests check them through SlabChannels and ChunksTouched, which are what callers ask.
    std::uint64_t _chunk_depth = 1;
    std::uint64_t _least_channels = 1;
};

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_REDUCE_PASS_PLAN_H_
