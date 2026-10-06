/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
*/

#ifndef CARTA_ZARR_BENCH_MODES_CUBE_H_
#define CARTA_ZARR_BENCH_MODES_CUBE_H_

// What the modes that read a cube share, and no other file needs.

#include "../mode.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>

namespace carta::zarr::bench {

// A runner over an image opened for it, with the bytes one of its elements is stored in.
class CubeRunner : public Runner {
public:
    CubeRunner(Image image, CubeAxes axes)
        : _image(std::move(image)), _axes(axes), _item_size(ItemSize(_image.descriptor().stored_type)) {}

    void Record(Row& result) const override { result.logical_bytes = result.elements * _item_size; }

protected:
    const Image& image() const noexcept { return _image; }
    const CubeAxes& axes() const noexcept { return _axes; }

private:
    Image _image;
    CubeAxes _axes;
    std::size_t _item_size;
};

// How many pixels ReadPixels writes: a plane's, or a spectrum's. What a buffer for it is sized by --
// sized for the larger of the two in either mode, a spectrum of a 32768-square image held a 4 GiB
// buffer per process for the few kilobytes it wrote.
std::size_t PixelsRead(const CubeAxes& axes, bool plane);

// A whole plane at `operation`'s channel, or when `plane` is false every channel at its pixel, read
// into `pixels`, which has room for PixelsRead(axes, plane): how many pixels were written.
Result<std::size_t> ReadPixels(const Image& image, const CubeAxes& axes, const Operation& operation, bool plane,
                               float* pixels, const ReadOptions& options);

// How an operation is written, for the modes that read whole planes over a run of channels.
std::string DescribeChannels(const Operation& operation);

}  // namespace carta::zarr::bench

#endif  // CARTA_ZARR_BENCH_MODES_CUBE_H_
