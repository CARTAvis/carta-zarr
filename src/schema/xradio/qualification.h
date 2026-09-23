/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_SCHEMA_XRADIO_QUALIFICATION_H_
#define CARTA_ZARR_SRC_SCHEMA_XRADIO_QUALIFICATION_H_

// Whether this profile will open a variable, and what it has to say about one it will not.
//
// One module rather than a rule per stage. Discovery decided what to list, and describing an image
// decided again whether it was an image at all, on a weaker rule -- so a listing could offer a
// variable that opening would refuse. Everything that needs the answer asks here, which is what
// makes an entry in the listing a statement about opening rather than about what a variable looks
// like.
//
// It takes a Store because that is the narrowest thing its tests can stand up: it reads array
// metadata and nothing else, and a Store over the in-memory transport costs a map entry. See ADR
// 0006.

#include "../../store.h"

#include "carta-zarr/descriptor.h"
#include "carta-zarr/result.h"

#include <array>
#include <optional>
#include <string_view>

namespace carta::zarr::internal::xradio {

// The axes an image of the sky plane carries. Every one of them is required: XRADIO writes optional
// coordinate arrays over the spatial pair alone -- right_ascension and declination are float64 over
// (l, m) and carry no type attribute -- so a rule keyed on "has l and m" offers those to a consumer
// as openable images. See ADR 0001.
//
// Shared with the probe, which checks that the dataset carries a well-formed coordinate array for
// each of them.
inline constexpr std::array<std::string_view, 5> kSkyAxes{"time", "frequency", "polarization", "l", "m"};

/**
 * What this profile decided about one node of a store.
 *
 * Three answers rather than two, because a node can fail to be an image in two different ways that
 * a consumer acts on differently. A coordinate array, a flag, and a beam table are not images and
 * there is nothing to say about them. A complex sky-plane variable is an image this profile will
 * not open, and a consumer that lists the dataset's images should see it with the reason attached.
 */
struct NodeQualification {
    // It belongs in the dataset's image listing, openable or not.
    bool listed = false;
    // This profile will open it. Never true unless it is listed.
    bool openable = false;
    // What there is to say: why a listed variable will not open, or why a node was passed over.
    // Absent for an ordinary image and for a node that is simply not an image.
    std::optional<Diagnostic> diagnostic;
    // Whether what there is to say is a defect in the store rather than a limit of this library.
    //
    // A complex sky-plane variable and an aperture-plane variable are well-formed things this
    // profile does not open; a variable whose extent disagrees with the coordinate it names, and a
    // node whose metadata will not parse, are the store being wrong about itself. The probe acts on
    // the difference: a dataset with nothing openable in it is a store this profile does not match,
    // unless something in it was malformed, which is a store it matched and found broken.
    bool malformed = false;
};

// Asked of every entry the store's inventory holds. What the node is, and whether an array's
// metadata parsed, the inventory has already decided; a node that will not parse is one of the
// answers here rather than a failure to arrive at one. The Store is for the coordinates an image
// names, which are other nodes.
NodeQualification QualifyNode(const Store& store, const NodeEntry& entry);

// The precondition of describing an image, as the error a caller should report when it is not met.
//
// Describing reads coordinate values and builds a descriptor; every one of those steps assumes the
// variable is an image this profile opens. Asking here rather than re-deciding is the whole point
// of the module: the answer is the one the listing gave.
Result<void> RequireQualified(const Store& store, std::string_view image_id);

}  // namespace carta::zarr::internal::xradio

#endif  // CARTA_ZARR_SRC_SCHEMA_XRADIO_QUALIFICATION_H_
