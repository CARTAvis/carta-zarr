/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

// The part of a pixel read that decides what to ask for, as opposed to asking.
//
// Separate from pixel_reader.cc so that it links without TensorStore. Translating a request into
// stored order, checking it against the descriptor, and counting what it selects are arithmetic over
// the descriptor, and a build that reads no pixels, such as a test's, still needs all three.

#include "zarr/pixel_selection.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <string>
#include <utility>

namespace carta::zarr::internal::zarr {

Result<PixelSelection> BuildSelection(const ImageDescriptor& descriptor, const ReadRequest& request,
                                      DestinationOrder order) {
    const auto rank = descriptor.axes.size();
    if (request.axes.size() != rank) {
        return Error{
            ErrorCode::invalid_slice,
            "Request has " + std::to_string(request.axes.size()) + " axes but the image has " + std::to_string(rank),
            descriptor.id};
    }

    PixelSelection selection;
    selection._start.assign(rank, 0);
    selection._count.assign(rank, 0);
    selection._stride.assign(rank, 1);
    selection._destination_to_stored.resize(rank);

    for (std::size_t logical = 0; logical < rank; ++logical) {
        const auto& axis = descriptor.axes.at(logical);
        const auto& range = request.axes.at(logical);
        if (range.stride == 0) {
            return Error{ErrorCode::invalid_slice, "Axis '" + axis.name + "' has a zero stride", descriptor.id};
        }
        if (range.count == 0) {
            return Error{ErrorCode::invalid_slice, "Axis '" + axis.name + "' selects no elements", descriptor.id};
        }
        // The last selected index has to fall inside the axis. This divides the room that is left
        // rather than multiplying out the span, because the span can overflow: a count and a stride
        // of about 2^32 each multiply to a small number, which would pass the check and then size a
        // buffer and drive the loops.
        if (range.start >= axis.length || range.count - 1 > (axis.length - 1 - range.start) / range.stride) {
            return Error{ErrorCode::invalid_slice,
                         "Axis '" + axis.name + "' request exceeds its length of " + std::to_string(axis.length),
                         descriptor.id};
        }
        const auto stored = axis.storage_index;
        if (stored >= rank) {
            return Error{ErrorCode::invalid_metadata, "Axis '" + axis.name + "' has an out-of-range storage index",
                         descriptor.id};
        }
        selection._start.at(stored) = range.start;
        selection._count.at(stored) = range.count;
        // A stride over one element means nothing, so any positive one was accepted above. The readers
        // take a signed stride, and one of 2^63 or more would reach TensorStore as a negative number,
        // so it is set to 1. Over two or more elements, the check above already keeps it below the
        // axis's length.
        selection._stride.at(stored) = range.count == 1 ? 1 : range.stride;
        selection._destination_to_stored.at(logical) = stored;
    }
    // In stored order the destination's axis 0 is the last stored dimension -- the one the array is
    // contiguous along -- so that it lands fastest, as it was written.
    if (order == DestinationOrder::stored) {
        for (std::size_t axis = 0; axis < rank; ++axis) {
            selection._destination_to_stored.at(axis) = rank - 1 - axis;
        }
    }

    // Counted here, once. A count that does not fit is an error rather than a wrapped value, which a
    // read would report as a destination too small for the request.
    std::uint64_t elements = 1;
    for (const auto count : selection._count) {
        if (elements > std::numeric_limits<std::uint64_t>::max() / count) {
            return Error{ErrorCode::invalid_slice, "Request selects more elements than can be counted", descriptor.id};
        }
        elements *= count;
    }
    selection._elements = elements;
    selection._destination_extent = elements;
    return selection;
}

Result<PixelSelection> BuildSelectionWithin(const ImageDescriptor& descriptor, const ReadRequest& part,
                                            const ReadRequest& whole) {
    auto selection = BuildSelection(descriptor, part, DestinationOrder::logical);
    if (!selection) {
        return selection;
    }
    if (whole.axes.size() != part.axes.size()) {
        return Error{ErrorCode::invalid_argument, "A part has a different number of axes from the read it is part of",
                     descriptor.id};
    }
    auto& built = selection.value();
    const auto rank = part.axes.size();
    built._destination_strides.assign(rank, 1);
    std::uint64_t running = 1;
    std::uint64_t offset = 0;
    std::uint64_t last = 0;
    for (std::size_t axis = 0; axis < rank; ++axis) {
        const auto& outer = whole.axes.at(axis);
        const auto& inner = part.axes.at(axis);
        const auto outer_stride = std::max<std::uint64_t>(1, outer.stride);
        // Every element of the part is an element of the whole: it starts on one, steps as the whole
        // does or over a single element, and ends inside it.
        const bool on_whole = inner.start >= outer.start && (inner.start - outer.start) % outer_stride == 0;
        const auto first = on_whole ? (inner.start - outer.start) / outer_stride : 0;
        if (!on_whole || (inner.count > 1 && inner.stride != outer.stride) || first >= outer.count ||
            inner.count > outer.count - first) {
            return Error{ErrorCode::invalid_argument,
                         "Axis '" + descriptor.axes.at(axis).name + "' of a part is not within the read it is part of",
                         descriptor.id};
        }
        built._destination_strides.at(axis) = running;
        offset += first * running;
        last += (inner.count - 1) * running;
        running *= outer.count;
    }
    built._destination_offset = offset;
    built._destination_extent = offset + last + 1;
    return selection;
}

std::vector<std::uint64_t> PixelSelection::DestinationStrides() const {
    std::vector<std::uint64_t> strides(_count.size(), 1);
    std::uint64_t running = 1;
    for (std::size_t axis = 0; axis < _destination_to_stored.size(); ++axis) {
        const auto stored = _destination_to_stored.at(axis);
        strides.at(stored) = _destination_strides.empty() ? running : _destination_strides.at(axis);
        running *= _count.at(stored);
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

}  // namespace carta::zarr::internal::zarr
