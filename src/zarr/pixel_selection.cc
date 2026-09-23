/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// The part of a pixel read that decides what to ask for, as opposed to asking.
//
// Separate from pixel_reader.cc so that it links without TensorStore. Translating a request over
// the logical axes into the stored order, checking it against the descriptor, and counting what it
// selects are arithmetic over the descriptor; a build that reads no pixels still needs all three,
// and until now got them from a hand-written copy in tests/support that could disagree with this
// one about the very ranges it exists to reject.

#include "zarr/pixel_selection.h"

#include <chrono>
#include <limits>
#include <string>
#include <utility>

namespace carta::zarr::internal::zarr {

Result<PixelSelection> BuildSelection(const ImageDescriptor& descriptor, const ReadRequest& request,
                                      DestinationOrder order) {
    const auto rank = descriptor.axes.size();
    if (request.axes.size() != rank) {
        return Error{ErrorCode::invalid_argument,
                     "Request has " + std::to_string(request.axes.size()) + " axes but the image has " +
                         std::to_string(rank),
                     descriptor.id};
    }

    PixelSelection selection;
    selection.start.assign(rank, 0);
    selection.count.assign(rank, 0);
    selection.stride.assign(rank, 1);
    selection.destination_to_stored.resize(rank);

    for (std::size_t logical = 0; logical < rank; ++logical) {
        const auto& axis = descriptor.axes.at(logical);
        const auto& range = request.axes.at(logical);
        if (range.stride == 0) {
            return Error{ErrorCode::invalid_argument, "Axis '" + axis.name + "' has a zero stride",
                         descriptor.id};
        }
        if (range.count == 0) {
            return Error{ErrorCode::invalid_argument, "Axis '" + axis.name + "' selects no elements",
                         descriptor.id};
        }
        // The last selected index is what has to fall inside the axis. It is compared by dividing
        // the room that is left rather than by multiplying out the span, because the span
        // overflows: a count and a stride of about 2^32 each multiply to a small number, which
        // passed this check and went on to size a buffer and drive the loops.
        if (range.start >= axis.length || range.count - 1 > (axis.length - 1 - range.start) / range.stride) {
            return Error{ErrorCode::invalid_argument,
                         "Axis '" + axis.name + "' request exceeds its length of " +
                             std::to_string(axis.length),
                         descriptor.id};
        }
        const auto stored = axis.storage_index;
        if (stored >= rank) {
            return Error{ErrorCode::invalid_metadata, "Axis '" + axis.name + "' has an out-of-range storage index",
                         descriptor.id};
        }
        selection.start.at(stored) = range.start;
        selection.count.at(stored) = range.count;
        selection.stride.at(stored) = range.stride;
        selection.destination_to_stored.at(logical) = stored;
    }
    // In stored order the destination's axis 0 is the last stored dimension -- the one the array is
    // contiguous along -- so that it lands fastest, as it was written.
    if (order == DestinationOrder::stored) {
        for (std::size_t axis = 0; axis < rank; ++axis) {
            selection.destination_to_stored.at(axis) = rank - 1 - axis;
        }
    }
    return selection;
}

std::vector<std::uint64_t> PixelSelection::DestinationStrides() const {
    std::vector<std::uint64_t> strides(count.size(), 1);
    std::uint64_t running = 1;
    for (const auto stored : destination_to_stored) {
        strides.at(stored) = running;
        running *= count.at(stored);
    }
    return strides;
}

Result<void> CheckReadControl(const ReadControl& control, std::string_view node) {
    if (control.cancellation_requested && control.cancellation_requested()) {
        return Error{ErrorCode::cancelled, "Pixel read was cancelled", std::string(node)};
    }
    if (std::chrono::steady_clock::now() >= control.deadline) {
        return Error{ErrorCode::cancelled, "Pixel read deadline expired", std::string(node)};
    }
    return {};
}

std::uint64_t SelectionElementCount(const PixelSelection& selection) {
    if (selection.count.empty()) {
        return 0;
    }
    std::uint64_t elements = 1;
    for (const auto value : selection.count) {
        if (value == 0 || elements > std::numeric_limits<std::uint64_t>::max() / value) {
            return 0;
        }
        elements *= value;
    }
    return elements;
}

}  // namespace carta::zarr::internal::zarr
