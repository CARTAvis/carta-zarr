/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_ZARR_PIXEL_SELECTION_H_
#define CARTA_ZARR_SRC_ZARR_PIXEL_SELECTION_H_

// What a pixel read asks for, as opposed to how it is served.
//
// Split from pixel_reader.h so that a caller deciding what to read does not take a dependency on
// the thing that reads it. A pass needs all of this and none of the Store, which is what lets a
// reduction run against pixels that were never on disk.

#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace carta::zarr::internal::zarr {

// Check cooperative cancellation and the deadline at a storage-operation boundary.
Result<void> CheckReadControl(const ReadControl& control, std::string_view node);

// The order a read's destination is laid out in. CONTEXT.md names both.
//
// An ordinary read wants logical order, because its destination is the caller's and its contract is
// a dense image in the order the library reports. A pass wants stored order, because a plane that
// arrives the way the store wrote it is never transposed, and the visitor pays for the difference in
// nothing but strides.
enum class DestinationOrder {
    logical,
    stored,
};

class PixelSelection;

// Translate a request over the logical axes into the stored axis order the array is written in,
// checking it against the descriptor on the way, with its destination laid out in `order`. Ranges are
// validated here rather than left to TensorStore so that an out-of-range request is an
// invalid_argument naming the axis, instead of an I/O error naming a domain.
//
// The only way to make a PixelSelection, so every one that exists is one this accepted.
Result<PixelSelection> BuildSelection(const ImageDescriptor& descriptor, const ReadRequest& request,
                                      DestinationOrder order);

/**
 * One hyperslab of an array, addressed in the array's own stored axis order, and the order its
 * destination is written in.
 *
 * The destination's order is carried alongside rather than applied by the caller because TensorStore
 * can fold it into the same copy that moves the decoded chunk into the destination: transposing here
 * costs a strided write, transposing afterwards costs a second full pass over the data.
 *
 * Only BuildSelection makes one, and it makes only well-formed ones: every vector has the array's
 * rank, the permutation is a permutation, no stride is zero, no count is zero, and the element count
 * fits. A reader used to check all of that again, because a struct of six public vectors could be
 * handed to it in any state; nothing could, and the check never fired. What a reader of pixels still
 * checks is whether the caller's destination holds this many elements, which is about the caller's
 * buffer rather than about the selection, and which a reader writing through a pointer must not take
 * on trust.
 */
class PixelSelection {
public:
    // All in stored axis order, one entry per stored dimension.
    const std::vector<std::uint64_t>& start() const noexcept { return _start; }
    const std::vector<std::uint64_t>& count() const noexcept { return _count; }
    const std::vector<std::uint64_t>& stride() const noexcept { return _stride; }

    // destination_to_stored()[i] is the stored dimension that the destination's axis i is. The
    // destination is dense with its axis 0 fastest-varying.
    //
    // It says how the destination is laid out and nothing else. The order is chosen when the
    // selection is built, and which logical axis a stored dimension is stays with the descriptor.
    const std::vector<std::size_t>& destination_to_stored() const noexcept { return _destination_to_stored; }

    // How many elements the selection produces. Never zero, and counted once, when it was built.
    std::uint64_t elements() const noexcept { return _elements; }

    // How far one step along each stored dimension moves in the destination, by stored dimension.
    std::vector<std::uint64_t> DestinationStrides() const;

private:
    PixelSelection() = default;
    friend Result<PixelSelection> BuildSelection(const ImageDescriptor& descriptor, const ReadRequest& request,
                                                 DestinationOrder order);

    std::vector<std::uint64_t> _start;
    std::vector<std::uint64_t> _count;
    std::vector<std::uint64_t> _stride;
    std::vector<std::size_t> _destination_to_stored;
    std::uint64_t _elements = 0;
};

}  // namespace carta::zarr::internal::zarr

#endif  // CARTA_ZARR_SRC_ZARR_PIXEL_SELECTION_H_
