/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <carta-zarr/carta_zarr.h>

int main() {
    // The build directory is not an image dataset. Whether it reads as no match or as a store that
    // cannot be read at all, the call has reached the installed library and answered.
    const auto probe = carta::zarr::ProbeSchema(".", carta::zarr::kXradioImageSchema);
    return !probe || probe.value().kind == carta::zarr::SchemaMatchKind::no_match ? 0 : 1;
}
