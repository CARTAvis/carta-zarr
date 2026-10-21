# Architecture and responsibility boundaries

carta-zarr is a library in its caller's process. The intended carta-backend
integration is under development; this repository does not contain a running
image server or frontend.

## System context

Where carta-zarr sits in a CARTA deployment, and where its datasets come from:

```{mermaid}
flowchart LR
    User(["Astronomer"]) -->|"web browser"| FE["carta-frontend"]
    FE <-->|"WebSocket, protobuf"| BE
    subgraph BE["carta-backend process"]
        Lib["carta-zarr"]
    end
    Input["FITS / CASA image"] -->|"convert"| XR["XRADIO"]
    XR -->|"write"| Disk[("Zarr v3 image dataset")]
    Lib -->|"read"| Disk
```

carta-zarr runs inside each backend process rather than as a service of its
own; any other C++ application links it the same way.

## Components

The application owns world coordinates, region definitions, scientific derived
quantities and its protocol. This library owns storage access and operations
expressible in pixel coordinates: counts, sums, spread, extrema and histograms.
A reduction belongs here only if it can be stated without the application's
terms: a flux density needs the beam, and a position-velocity image needs world
coordinates, so both stay with the application. carta-backend, for example,
builds casacore coordinates and the `ImageInterface` CARTA reads from what the
descriptor reports.

Inside it:

```{mermaid}
flowchart TB
    subgraph G_api["Public API"]
        API["Context / Dataset / Image"]
    end
    subgraph G_schema["XRADIO schema"]
        Profile["XRADIO schema profile"]
    end
    subgraph G_compute["Reads and reductions"]
        Pieces["Piece planner: bounded reads"]
        Pass["Pass: chunk traversal and reduction"]
        Source["PixelSource / StorePixelSource"]
    end
    subgraph G_storage["Zarr storage"]
        Store["Store: hierarchy, metadata, verification"]
        Transport["Filesystem Transport"]
        Reader["Zarr readers / TensorStore"]
    end
    API -->|"recognize and describe"| Profile
    Profile -->|"inspect metadata and coordinates"| Store
    API -->|"Read / Prefetch"| Pieces
    API -->|"statistics / histograms"| Pass
    Pieces -->|"pixels and flags"| Source
    Pass -->|"pixels and flags"| Source
    Source -->|"read selected values"| Store
    Store -->|"metadata, listing, array paths"| Transport
    Store -->|"numeric reads"| Reader
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
  describes flags, coordinates and observation information. It reads XRADIO
  image schema 0.0.2 only, and refuses an image dataset that declares no
  version or another one, or that declares 0.0.2 without the fields 0.0.2
  requires, with a diagnostic saying why.
- **Piece planner** cuts requests into pieces constrained by the read budget.
  **Pass** traverses chunks for reductions, with occupancy restricting work to
  chunks actually touched by the selected regions.
- **TensorStore** decodes numeric arrays, transforms selections and provides
  decoded-chunk caching. It is absent from public headers.
- **Context** shares execution resources. `ReadAhead` schedules background
  `Prefetch` calls for an animation's next distinct run of chunks.

See the [glossary](glossary.md) for domain terminology and the repository's
[architecture decisions](https://github.com/CARTAvis/carta-zarr/tree/dev/docs/adr)
for design alternatives and measured tradeoffs.
