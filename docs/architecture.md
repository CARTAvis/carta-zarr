# Architecture and responsibility boundaries

carta-zarr is a library in its caller's process. The intended carta-backend
integration is under development; this repository does not contain a running
image server or frontend.

```{mermaid}
flowchart LR
    Input["FITS / CASA"] -->|"convert"| XR["XRADIO"]
    XR -->|"write"| Disk[("Local Zarr v3 dataset")]
    Caller["carta-backend / C++ application"] -->|"public API"| Lib["carta-zarr"]
    Lib -->|"read"| Disk
    Lib -->|"pixels, descriptors, totals"| Caller
```

The application owns world coordinates, region definitions, scientific derived
quantities and its protocol. This library owns storage access and operations
expressible in pixel coordinates: counts, sums, spread, extrema and histograms.

```{mermaid}
flowchart TB
    API["Context / Dataset / Image"] -->|"recognise and describe"| Profile["XRADIO schema profile"]
    Profile -->|"inspect metadata and coordinates"| Store["Store: hierarchy, metadata, verification"]
    API -->|"Read / Prefetch"| Pieces["Piece planner: bounded reads"]
    API -->|"statistics / histograms"| Pass["Pass: chunk traversal and reduction"]
    Pieces -->|"pixels and flags"| Source["PixelSource / StorePixelSource"]
    Pass -->|"pixels and flags"| Source
    Source -->|"read selected values"| Store
    Store -->|"metadata, listing, array paths"| Transport["Filesystem Transport"]
    Store -->|"numeric reads"| Reader["Zarr readers / TensorStore"]
    Transport -->|"filesystem access"| Disk[("Zarr metadata and chunks")]
    Reader -->|"decode chunks and cache"| Disk
```

- **Transport** supplies raw metadata, node names, local array paths and stored
  size. It is the extension boundary for where bytes live, though remote pixel
  transport is not implemented.
- **Store** interprets Zarr metadata and caches an inventory and parsed arrays.
  Consolidated metadata is checked against the array's own document before its
  values are read.
- **Schema profile** decides what constitutes an image, maps axes by role and
  describes flags, coordinates and observation information.
- **Read planning** cuts requests into pieces constrained by the read budget.
  **Pass** traverses chunks for reductions, with occupancy restricting work to
  chunks actually touched by the selected regions.
- **TensorStore** decodes numeric arrays, transforms selections and provides
  decoded-chunk caching. It is absent from public headers.
- **Context** shares execution resources. `ReadAhead` schedules background
  `Prefetch` calls for an animation's next distinct run of chunks.

Tests can replace `PixelSource` with synthetic pixels and Transport with an
in-memory implementation, exercising traversal and metadata without production
storage. The in-memory metadata transport cannot serve filesystem pixel reads.

See the [glossary](glossary.md) for domain terminology and the repository's
[architecture decisions](https://github.com/CARTAvis/carta-zarr/tree/dev/docs/adr)
for design alternatives and measured tradeoffs.
