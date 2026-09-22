/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// Companion to value_reader_stub.cc: the schema profile tests link without TensorStore. The
// in-memory transport holds no array data, so a pixel read has no answer to give and
// unsupported_transport is the honest one.
//
// Only the reads are stubbed. Building a selection, counting what it holds and checking the read
// controls are arithmetic rather than reads, so they come from src/zarr/pixel_selection.cc -- the
// same code the real build runs. They used to be copied out by hand here, which meant a build that
// reads no pixels could disagree with the one that does about which requests are legal.

#include "zarr/pixel_reader.h"

#include <string>

namespace carta::zarr::internal::zarr {

Result<void> ReadFloat32(const std::filesystem::path&, const StoreContextPtr&, std::string_view node,
                         std::string_view, const PixelSelection&, float*, std::size_t, const ReadControl&) {
    return Error{ErrorCode::unsupported_transport, "This build reads no pixels", std::string(node)};
}

Result<void> ReadMaskBytes(const std::filesystem::path&, const StoreContextPtr&, std::string_view node,
                           std::string_view, const PixelSelection&, std::uint8_t*, std::size_t, const ReadControl&) {
    return Error{ErrorCode::unsupported_transport, "This build reads no pixels", std::string(node)};
}

}  // namespace carta::zarr::internal::zarr
