/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_SRC_STORE_H_
#define CARTA_ZARR_SRC_STORE_H_

#include "carta-zarr/descriptor.h"
#include "carta-zarr/read.h"
#include "carta-zarr/result.h"

#include "memo.h"
#include "zarr/array_metadata.h"
#include "zarr/transport.h"

#include <nlohmann/json.hpp>

#include <chrono>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace carta::zarr::internal {

// Defined in zarr/pixel_reader.h. Forward declared so that store.h stays free of the reader that
// pulls in TensorStore, the same way StoreContext is kept opaque below.
namespace zarr {
struct PixelSelection;
}

// Everything a Store remembers for the lifetime of its read-only view.
//
// Each table carries its own lock, and that is load-bearing rather than incidental: a Memo holds its
// lock across the computation, and these computations nest -- reading array metadata reads node
// metadata, listing nodes reads each node's metadata. One shared lock would deadlock. The nesting
// forms a DAG, and a reader added later must not add an edge back:
//
//     string_arrays ---> array_metadata ---> node_metadata ---> transport
//     listed_nodes ------------------------> node_metadata ---> transport
//     double_arrays -------------------------------------------> transport
//
// These locks are not the store's concurrency story. Descriptor construction is serialized a level
// up, by the mutex Dataset::OpenImage holds across the whole of DescribeSchema, so today only two
// concurrent Image::ReadBeams calls on one Dataset reach these tables at the same time. When pixel
// reads arrive that changes, and holding a lock across a TensorStore read becomes worth revisiting
// -- in here, rather than in five hand-written places as before.
struct StoreCaches {
    Memo<std::string, Result<nlohmann::json>> node_metadata;
    Memo<std::string, Result<zarr::ArrayMetadata>> array_metadata;
    Lazy<Result<std::vector<std::string>>> listed_nodes;
    Memo<std::string, Result<std::vector<double>>> double_arrays;
    Memo<std::string, Result<std::vector<std::string>>> string_arrays;
};

// Defined in zarr/store_context.h. Kept opaque here so that including store.h does not pull in
// TensorStore; a null pointer means "use TensorStore's default resources".
class StoreContext;
using StoreContextPtr = std::shared_ptr<const StoreContext>;

// A Store is the Dataset-scoped storage session that interprets a Zarr hierarchy over a Transport.
// It owns everything that turns bytes into meaning -- node path validation, JSON parsing,
// consolidated metadata, array metadata, storage layout, data locations, and the caches -- so that
// every read uses one consistent session. Nothing above this module learns where the bytes came from.
class Store {
public:
    // `consolidated_metadata` is the root's copy of its children's documents, which the store takes
    // over as the node metadata it has already read rather than keeping beside it.
    Store(TransportPtr transport, nlohmann::json root_attributes,
          std::map<std::string, nlohmann::json> consolidated_metadata, bool has_consolidated_metadata,
          StoreContextPtr context);

    // The root group's attributes, or an empty object when it declares none.
    const nlohmann::json& RootAttributes() const noexcept;

    // These hand back what the store holds, not a copy of it: the reference is good for as long as
    // the store is, and a caller that wants its own copy says so.
    const Result<nlohmann::json>& ReadNodeMetadata(std::string_view node) const;
    const Result<zarr::ArrayMetadata>& ReadArrayMetadata(std::string_view node) const;
    // Every node in the hierarchy, sorted, each one's metadata read and parsed. The names are what a
    // caller walks; what a node holds it asks for by name, which is already in hand by then.
    const Result<std::vector<std::string>>& ListNodes() const;
    // Values in C order, flattened. The rank is in the node's ArrayMetadata; ArrayView addresses
    // them by dimension name rather than by offset.
    Result<std::vector<double>> ReadNumericArray(std::string_view node) const;
    Result<std::vector<std::string>> ReadStringArray1D(std::string_view node) const;
    Result<StorageLayout> ReadStorageLayout(std::string_view node) const;

    // How many bytes this store occupies where it lives. Handed straight to the transport, which is
    // the only thing that knows; Store fronts it for the same reason it fronts every other
    // transport question, so that nothing above here learns where the bytes came from.
    Result<std::uint64_t> StoredSizeBytes(std::chrono::steady_clock::time_point deadline) const;

    // Pixel reads are deliberately uncached here: a slab is requested once and is far larger than
    // anything the metadata tables hold. Reuse belongs in TensorStore's chunk cache, which already
    // works at chunk granularity and is sized by the consumer's Context.
    //
    // One function rather than two, because reading an image's pixels and reading its flag differed
    // only in the element type: both find the array's metadata, ask the transport where its bytes
    // are, and hand both to the reader. `float` is pixels, converted from whatever the array holds;
    // `std::uint8_t` is a flag, one byte an element, true meaning a good pixel. Instantiated for
    // those two in store.cc and for nothing else.
    template <typename T>
    Result<void> ReadPixelsInto(std::string_view node, const zarr::PixelSelection& selection, T* destination,
                                std::size_t destination_elements, const ReadOptions& options) const;

private:
    Result<std::filesystem::path> ResolveArrayDirectory(std::string_view node) const;
    Result<std::vector<double>> ReadNumericArrayUncached(std::string_view node) const;
    Result<std::vector<std::string>> ReadStringArray1DUncached(std::string_view node) const;

    TransportPtr _transport;
    nlohmann::json _root_attributes;
    // The names the root's copy accounted for. The documents themselves went into the node metadata
    // table at construction, so there is one place a node's metadata is looked up rather than two.
    std::vector<std::string> _consolidated_nodes;
    bool _has_consolidated_metadata = false;
    StoreContextPtr _context;
    // Held indirectly so that Store stays movable: the tables own mutexes and cannot be moved.
    std::shared_ptr<StoreCaches> _caches;
};

// Open a store on the local filesystem.
Result<Store> OpenStore(std::string_view location, StoreContextPtr context = {});

// Open a store over an already-built transport. This is the seam tests enter through.
Result<Store> OpenStore(TransportPtr transport, StoreContextPtr context = {});

// How large a dataset is, and whether that number was measured or inferred.
//
// The store's own size when the transport can report it before the deadline, and otherwise the
// total uncompressed size of every array, marked as an upper bound. Which of the two a caller gets
// is this module's decision rather than the facade's: it is one question -- how much room does this
// take -- and answering half of it up there is what had the facade walking a directory of its own.
//
// Reports invalid_metadata for a store with no arrays at all, and for a size that overflows.
Result<DatasetSize> DatasetSizeBytes(const Store& store, std::chrono::steady_clock::time_point deadline);

}  // namespace carta::zarr::internal

#endif  // CARTA_ZARR_SRC_STORE_H_
