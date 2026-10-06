/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_TESTS_SUPPORT_SYNTHETIC_SLAB_SOURCE_H_
#define CARTA_ZARR_TESTS_SUPPORT_SYNTHETIC_SLAB_SOURCE_H_

#include "pixel_source.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

namespace carta::zarr::testing {

/**
 * The pass's second adapter: pixels computed from their own coordinates, never written down.
 *
 * The filesystem adapter is the one production uses, and until this existed it was the only thing
 * that could answer a pass at all -- so every question about how the pass walks a cube had to be
 * asked of a directory tree. That put two limits on what could be asked. A fixture large enough to
 * reach the parallel paths had to be committed, so there is exactly one and it is 512x520x4; and
 * the pass's own reads could not be observed, so the tests count how many times a result was handed
 * over and reason backwards to what must have been read, and say in three places that they stop
 * testing anything if the chunk shape ever changes.
 *
 * This one states a cube of any size in a line, and counts what was asked of it.
 */
class SyntheticPixelSource final : public carta::zarr::internal::PixelSource {
public:
    // The value at one set of logical coordinates, indexed the way ImageDescriptor::axes is.
    using Formula = std::function<float(const std::vector<std::uint64_t>& logical)>;
    // Whether the pixel at those coordinates is good. Only consulted when the read applies a mask.
    using FlagFormula = std::function<bool(const std::vector<std::uint64_t>& logical)>;

    SyntheticPixelSource(const carta::zarr::ImageDescriptor& descriptor, const carta::zarr::ChunkGeometry& geometry,
                         Formula formula)
        : _descriptor(&descriptor), _geometry(&geometry), _formula(std::move(formula)) {}

    void set_flags(FlagFormula flags) { _flags = std::move(flags); }

    // Serve one value everywhere, without the per-pixel call the formula costs. For the tests that
    // are about how the pass splits a cube rather than about what is in it: at the sizes those use,
    // a std::function per pixel is most of the test's running time and none of its point.
    void set_constant(float value) {
        _constant = value;
        _has_constant = true;
    }

    // Fail the nth pixel read, counting from one, so that a caller's error path is reachable.
    void fail_read(std::uint64_t nth, carta::zarr::ErrorCode code) {
        _fail_at = nth;
        _fail_code = code;
    }

    // The same for the flag. Separate because the order the two are read in is a contract: an
    // ordinary read reads the flag first so that a mask it cannot get leaves the caller's
    // destination alone, and only a source that can fail one without the other can show it.
    void fail_mask_read(std::uint64_t nth, carta::zarr::ErrorCode code) {
        _mask_fail_at = nth;
        _mask_fail_code = code;
    }

    carta::zarr::Result<void> ReadPixels(const carta::zarr::internal::zarr::PixelSelection& selection,
                                         carta::zarr::BufferView<float> destination,
                                         const carta::zarr::ReadControl&) const override {
        ++_pixel_reads;
        _pixel_destinations.push_back(destination.size);
        if (_fail_at != 0 && _pixel_reads == _fail_at) {
            return carta::zarr::Error{_fail_code, "The synthetic source was told to fail here", "SKY"};
        }
        Record(selection);
        if (_has_constant) {
            const std::uint64_t total = selection.elements();
            if (auto fits = Fits(total, destination.size); !fits) {
                return fits;
            }
            _elements += total;
            std::fill(destination.data, destination.data + total, _constant);
            return carta::zarr::Result<void>{};
        }
        return Fill(selection, destination.size, [&](const std::vector<std::uint64_t>& logical, std::size_t at) {
            destination.data[at] = _formula(logical);
        });
    }

    carta::zarr::Result<void> ReadMask(const carta::zarr::internal::zarr::PixelSelection& selection,
                                       carta::zarr::BufferView<std::uint8_t> destination,
                                       const carta::zarr::ReadControl&) const override {
        ++_mask_reads;
        _mask_destinations.push_back(destination.size);
        _mask_selections.push_back(selection.elements());
        if (_mask_fail_at != 0 && _mask_reads == _mask_fail_at) {
            return carta::zarr::Error{_mask_fail_code, "The synthetic source was told to fail this flag read", "FLAG"};
        }
        return Fill(selection, destination.size, [&](const std::vector<std::uint64_t>& logical, std::size_t at) {
            // A flag byte is true for a flagged pixel, as XRADIO writes it: the opposite of the formula.
            destination.data[at] = _flags && !_flags(logical) ? 1 : 0;
        });
    }

    // What the pass asked for. The chunk counts are the assertion the fixture-driven tests cannot
    // make: a walk that decodes a chunk twice reads the same cell twice, whatever its answer.
    std::uint64_t pixel_reads() const { return _pixel_reads; }
    std::uint64_t mask_reads() const { return _mask_reads; }
    // How long a destination each read was handed, in the order they were asked. What a caller
    // hands across the seam is its own buffer, so this is what says whether it handed that or
    // restated the selection's size -- the one mistake a check at the seam cannot see from inside.
    const std::vector<std::size_t>& pixel_destinations() const { return _pixel_destinations; }
    const std::vector<std::size_t>& mask_destinations() const { return _mask_destinations; }
    // How many elements each flag read selected.
    const std::vector<std::uint64_t>& mask_selections() const { return _mask_selections; }
    std::uint64_t elements_read() const { return _elements; }
    // How many chunk cells were touched more than once, and the most any one of them was touched.
    std::size_t chunks_touched() const { return _chunk_hits.size(); }
    std::uint64_t most_hits_on_one_chunk() const {
        std::uint64_t most = 0;
        for (const auto& entry : _chunk_hits) {
            most = std::max(most, entry.second);
        }
        return most;
    }

private:
    // The destination is dense with its axis 0 fastest, laid out as
    // PixelSelection::destination_to_stored says -- logical order for a plain read, stored order for
    // a pass.
    //
    // Worked out here rather than asked of PixelSelection::DestinationStrides, deliberately. This
    // adapter writes pixels where these strides say and the pass reads them where the selection's
    // say; were both the same function, a mistake in it would be made twice, consistently, and
    // every test would pass while the real writer -- TensorStore, which lays the destination out
    // from the permutation and not from that function -- disagreed. Deriving it again is the
    // cross-check.
    std::vector<std::uint64_t> DestinationStrides(const carta::zarr::internal::zarr::PixelSelection& selection) const {
        const auto rank = selection.count().size();
        std::vector<std::uint64_t> strides(rank, 1);
        std::uint64_t running = 1;
        for (std::size_t axis = 0; axis < rank; ++axis) {
            const auto stored = selection.destination_to_stored().at(axis);
            strides.at(stored) = running;
            running *= selection.count().at(stored);
        }
        return strides;
    }

    // The seam's guard, as the TensorStore reader has one: the destination is the caller's, and
    // this writes through it as far as the selection reaches. Longer than the selection is the
    // ordinary case -- a piece is handed the rest of the caller's buffer -- so only shorter is
    // refused.
    static carta::zarr::Result<void> Fits(std::uint64_t total, std::size_t destination_elements) {
        if (total > destination_elements) {
            return carta::zarr::Error{carta::zarr::ErrorCode::invalid_argument, "Destination buffer is too small",
                                      "SKY"};
        }
        return carta::zarr::Result<void>{};
    }

    template <typename Write>
    carta::zarr::Result<void> Fill(const carta::zarr::internal::zarr::PixelSelection& selection,
                                   std::size_t destination_elements, Write&& write) const {
        const auto rank = selection.count().size();
        const std::uint64_t total = selection.elements();
        if (auto fits = Fits(total, destination_elements); !fits) {
            return fits;
        }
        _elements += total;

        // stored_to_logical comes from the descriptor, because the selection does not say it: its
        // permutation is the destination's layout, which for a pass is not logical order.
        std::vector<std::size_t> stored_to_logical(rank, 0);
        for (std::size_t logical = 0; logical < _descriptor->axes.size(); ++logical) {
            stored_to_logical.at(_descriptor->axes.at(logical).storage_index) = logical;
        }
        const auto strides = DestinationStrides(selection);

        std::vector<std::uint64_t> index(rank, 0);
        std::vector<std::uint64_t> logical(rank, 0);
        for (std::uint64_t n = 0; n < total; ++n) {
            std::size_t at = 0;
            for (std::size_t stored = 0; stored < rank; ++stored) {
                at += static_cast<std::size_t>(index.at(stored) * strides.at(stored));
                logical.at(stored_to_logical.at(stored)) =
                    selection.start().at(stored) + (index.at(stored) * selection.stride().at(stored));
            }
            write(logical, at);
            for (std::size_t stored = rank; stored-- > 0;) {
                if (++index.at(stored) < selection.count().at(stored)) {
                    break;
                }
                index.at(stored) = 0;
            }
        }
        return {};
    }

    void Record(const carta::zarr::internal::zarr::PixelSelection& selection) const {
        const auto rank = selection.count().size();
        // The chunk shape is in logical order; the selection is in stored order.
        std::vector<std::uint64_t> chunk(rank, 1);
        for (std::size_t logical = 0; logical < _descriptor->axes.size(); ++logical) {
            const auto stored = _descriptor->axes.at(logical).storage_index;
            chunk.at(stored) = std::max<std::uint64_t>(1, _geometry->chunk_shape.at(logical));
        }

        std::vector<std::uint64_t> first(rank, 0);
        std::vector<std::uint64_t> last(rank, 0);
        for (std::size_t stored = 0; stored < rank; ++stored) {
            const auto begin = selection.start().at(stored);
            const auto end = begin + ((selection.count().at(stored) - 1) * selection.stride().at(stored));
            first.at(stored) = begin / chunk.at(stored);
            last.at(stored) = end / chunk.at(stored);
        }

        std::vector<std::uint64_t> cell = first;
        while (true) {
            ++_chunk_hits[cell];
            std::size_t stored = rank;
            while (stored-- > 0) {
                if (++cell.at(stored) <= last.at(stored)) {
                    break;
                }
                cell.at(stored) = first.at(stored);
            }
            if (stored == static_cast<std::size_t>(-1)) {
                break;
            }
        }
    }

    const carta::zarr::ImageDescriptor* _descriptor;
    const carta::zarr::ChunkGeometry* _geometry;
    Formula _formula;
    FlagFormula _flags;
    float _constant = 0.0F;
    bool _has_constant = false;
    std::uint64_t _mask_fail_at = 0;
    carta::zarr::ErrorCode _mask_fail_code = carta::zarr::ErrorCode::io_error;
    std::uint64_t _fail_at = 0;
    carta::zarr::ErrorCode _fail_code = carta::zarr::ErrorCode::io_error;

    mutable std::uint64_t _pixel_reads = 0;
    mutable std::uint64_t _mask_reads = 0;
    mutable std::uint64_t _elements = 0;
    mutable std::vector<std::size_t> _pixel_destinations;
    mutable std::vector<std::size_t> _mask_destinations;
    mutable std::vector<std::uint64_t> _mask_selections;
    mutable std::map<std::vector<std::uint64_t>, std::uint64_t> _chunk_hits;
};

}  // namespace carta::zarr::testing

#endif  // CARTA_ZARR_TESTS_SUPPORT_SYNTHETIC_SLAB_SOURCE_H_
