# carta-zarr

A read-only C++17 library for XRADIO images stored as Zarr v3. It is written as the Zarr reader of
[carta-backend](https://github.com/CARTAvis/carta-backend), CARTA's image server, and kept apart
from it so that CARTA and other C++ consumers can read these datasets without depending on
TensorStore. It depends on neither casacore nor CARTA's protobuf either.

**Status:** not yet released. carta-backend's Zarr support, the consumer this library is shaped by,
is under development and not yet in carta-backend's main branch or any of its releases; the
carta-backend this README describes is that work.

TensorStore is a private implementation detail: it is linked into this library and never appears in
its public headers, which need only the C++ standard library.

## What it reads

- **Image datasets** written by [XRADIO](https://github.com/casangi/xradio) 1.2.x: one Zarr v3 group
  holding sky-plane image variables and the coordinates they share. Schema id `xradio.image`,
  schema version `1.2` -- the layout this library was checked against, not a schema version XRADIO
  publishes. The conformance fixture is written by XRADIO 1.2.3 itself (see
  [tests/data/README.md](tests/data/README.md)), and following a later XRADIO means pinning that
  generator to it and seeing which assumption breaks.
- **Images**: real-valued variables carrying all five of `time`, `frequency`, `polarization`, `l`
  and `m`. Complex and aperture-plane (`u`, `v`) variables are listed with a diagnostic saying why
  they cannot be opened.
- **Local filesystem stores.** Every node carries its own metadata document; a root
  `consolidated_metadata` copy, which is what zarr-python writes for a consolidated dataset, is used
  to answer for them without reading each one. An array's own document is still read, and held to
  the copy, before any value is read from it
  ([ADR 0017](docs/adr/0017-an-array-is-held-to-its-own-document-once.md)).

Writing is out of scope, and so is anything but the `xradio.image` profile.

## Getting an image dataset

XRADIO converts FITS and CASA images itself, and what it writes is what this library reads:

```python
from xradio.image import open_image, write_image

xds = open_image("cube.fits")      # a FITS or CASA image, or a dict naming several
write_image(xds, "cube.zarr", out_format="zarr")
```

XRADIO's FITS reader is strict: every axis needs `CTYPE` and `CUNIT`, the Stokes axis included, and
the header must carry `LONPOLE`, `LATPOLE`, the `PC` matrix and `DATE-OBS`. `open_image` imports
the casacore layer even for FITS, and on macOS that means casatools.
[tests/data/generate_conformance_fixtures.py](tests/data/generate_conformance_fixtures.py) does
exactly this with XRADIO 1.2.3 pinned, and runs as it is under `uv`.

To store a cube in another chunk layout, or to make one of a given size when there is none to hand,
[tools/zarr-bench/generate.py](tools/zarr-bench/README.md) rewrites a real cube or synthesizes one.

## Where the library ends

The library decides how to reach the pixels; carta-backend decides what the numbers mean.

Inside, from the bottom up: a **transport** gives raw access to a store's nodes -- the local
filesystem, or memory in the tests, with a remote store as one more transport ([ADR 0004](docs/adr/0004-store-seam-at-the-transport.md)); the
**store** reads the Zarr hierarchy and array metadata and knows nothing of images; a **schema
profile** recognises an XRADIO dataset and describes its images by the role of each axis
([ADR 0007](docs/adr/0007-the-schema-profile-table-keeps-its-dispatch.md), [ADR 0012](docs/adr/0012-an-axis-is-reached-by-its-role.md)); and **reads and reductions** walk the chunks those images are stored in,
through TensorStore, which goes no further than this library.

carta-backend keeps everything that gives the numbers a meaning: casacore coordinates and the
`ImageInterface` it shows CARTA ([ADR 0002](docs/adr/0002-casacore-coordinates-without-fits.md)), regions and world coordinates, the statistics derived
from totals, position-velocity images, moments and the protocol. A reduction belongs in this library
only if it can be stated without any of those: counts, sums, extrema and histograms can be; a flux
density, which needs the beam, or a position-velocity image, which needs world coordinates, cannot.

## Requirements

Install these before the first configure:

| Tool | Version | Why |
|---|---|---|
| CMake | 3.24 or newer | TensorStore's own minimum. |
| C++ compiler | GCC 10, Clang 8 or Apple Xcode 11.3.1, or newer | C++17, and the compilers TensorStore supports. |
| Python 3 | 3.10 or newer | TensorStore's build turns its Bazel rules into CMake with a Python script at configure time, and the script uses `match`. Only the interpreter is needed. |
| NASM | any | TensorStore's configure enables the NASM assembler for codecs it declares, although carta-zarr compiles none of them. |
| patch | any | TensorStore patches some of its dependencies as it unpacks them. |

The first configure also needs network access. It downloads TensorStore 0.1.84 and the 42 libraries
TensorStore declares -- 169 MB of archives, 1.4 GB unpacked -- although only ten of them are compiled:
TensorStore itself, Abseil, riegeli, re2, zstd, zlib, blosc, snappy, lz4 and nlohmann_json. Nothing
else needs to be installed: those ten are linked into the library statically, and it exports none of
their symbols, so they cannot collide with copies a consumer links itself. A release build tree takes
about 1.3 GB.

Built and tested with:

- macOS 26 on Apple silicon: Apple clang 21, CMake 4.0, Python 3.13, NASM 3.01, from Homebrew
  (`brew install cmake nasm python`) and the Xcode command line tools.
- Ubuntu 24.04: GCC 13.3, CMake 3.28, Python 3.12, NASM 2.16
  (`sudo apt-get install build-essential cmake python3 nasm patch`).
- AlmaLinux 8.10 with the distribution's gcc-toolset-13 (GCC 13.3), CMake 3.26 and Python 3.12, with
  PowerTools enabled for NASM (`sudo dnf install gcc-toolset-13-gcc-c++ cmake python3.12 nasm patch`).
  Its own GCC 8 and Python 3.6 are too old: build in `scl enable gcc-toolset-13 bash`, and configure
  with `-DPython3_EXECUTABLE=/usr/bin/python3.12`, since Python 3.6 is installed too.

The storage benchmark, `-DCARTA_ZARR_BUILD_BENCH=ON`, also needs [`uv`](https://docs.astral.sh/uv/);
see [tools/zarr-bench/README.md](tools/zarr-bench/README.md).

## Build, test, install

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

The first configure builds TensorStore, which takes a while; later ones do not.

### Without network access

Build farms such as Launchpad and most RPM builders have no network. Configure once on a machine
that has, pack the 43 archives that configure downloaded, and point the offline configure at them:

```bash
# Where there is network, after a configure into build/:
(cd build/_deps && find . -maxdepth 4 -path '*-subbuild/*-populate-prefix/src/*' -type f \
    | tar czf ../../carta-zarr-deps.tar.gz -T -)

# Where there is not:
mkdir deps && tar xzf carta-zarr-deps.tar.gz -C deps
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFETCHCONTENT_BASE_DIR=$PWD/deps
```

Each archive is checked against the hash TensorStore pins and then not downloaded again. It is the
archives that move, not the unpacked sources: TensorStore patches each source as it unpacks it, and
the patches name the build tree they were made for, so neither copying another tree's `_deps` nor
`FETCHCONTENT_SOURCE_DIR_<name>` gives a tree that configures.

## Licence of what it links in

The ten libraries linked into `libcarta-zarr`, and the header-only `half`, are listed with their
licences in `THIRD_PARTY_NOTICES`, which the build writes from their sources and installs beside
`LICENSE`. A binary distribution of the library carries that file.

## Using it from CMake

```cmake
find_package(CartaZarr 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

`tests/consumer` is exactly this, built against the install tree by the `carta-zarr-package-consumer`
test, so the packaging stays consumable.

Before 1.0 a minor release may change the ABI, so `0.1` finds 0.1.x and nothing else, and the soname
is `libcarta-zarr.so.0.1`. A binary is rebuilt when the minor moves. From 1.0 on that boundary is the
major ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)).

## Using the API

```c++
#include <carta-zarr/carta_zarr.h>

auto context = carta::zarr::Context::Create();          // one per process is enough
auto dataset = carta::zarr::Dataset::Open(context.value(), "/path/to/image.zarr");
auto image = dataset.value().OpenImage(*dataset.value().descriptor().default_image_id);

const auto& axes = image.value().descriptor().axes;
carta::zarr::ReadRequest request;                        // one Range per axis, each defaulting to
request.axes.resize(axes.size(), {0, 1, 1});             // the first plane of that axis
request.axes[*carta::zarr::AxisIndex(axes, carta::zarr::AxisRole::spatial_x)] = {0, 512, 1};
request.axes[*carta::zarr::AxisIndex(axes, carta::zarr::AxisRole::spatial_y)] = {0, 512, 1};
std::vector<float> pixels(512 * 512);
auto written = image.value().Read(request, {pixels.data(), pixels.size()});
```

`carta::zarr::ProbeSchema(path, carta::zarr::kXradioImageSchema)` answers "is this an image dataset
this library reads?" without opening an image: a `SchemaMatchKind` of match, no match, or a match
that is malformed. `Dataset::Open` lists the images in one, and `ImageEntry::openable` says which of
them will open; an openable entry also carries the axes it will open with, so a consumer can choose
among them without opening any.

### What an image can do

| Call | What it does |
|---|---|
| `Image::Read` | Pixels of any range along each axis, as above. |
| `Image::ReduceSpectral` | Per-channel statistics of many regions in one pass -- the pixel count, NaN count, sum, sum of squares, spread, minimum and maximum -- streamed to a callback block by block. What spectral profiles and region statistics are made from. |
| `Image::ComputeHistogram` | A histogram of each selected plane, over a range and bin count the caller gives. |
| `Image::ComputeCubeHistogram` | One histogram of the whole selection when its range is not known beforehand, with progress. |
| `Image::ReadBeams` | The restoring beam fitted on each plane. |
| `Image::Prefetch`, `Image::DecodedChunkBytes` | Decode a read's chunks into the cache without reading them out, and say how much cache one chunk takes. |
| `ReadAhead::For` | Read ahead of a playing animation, so that the frame entering a new run of chunks does not stall ([ADR 0016](docs/adr/0016-reading-ahead-is-work-a-caller-holds.md)). |
| `Dataset::Size` | The store's size on disk within a timeout, or else the size its metadata declares, and which of the two it is. |

The headers under [include/carta-zarr](include/carta-zarr) are the reference: every type, field and
call is documented where it is declared -- `carta_zarr.h` for the handles, `read.h` for a read's
request and options, `reduce.h` for reductions and histograms, `read_ahead.h` for read-ahead, and
`descriptor.h` for what a dataset and an image describe.

### What the API promises

- **Everything is a `Result<T>`.** A failure carries an `ErrorCode`, a message and the node it came
  from. No public entry point throws: a malformed value in the metadata or a length no allocator can
  serve is reported as an error like any other. `Result` is `[[nodiscard]]` and reads like
  `std::expected`: `has_value`, `value`, `value_or`, `*`, `->` and `error`.
- **A handle always refers to something.** `Context`, `Dataset`, `Image` and `CachePool` come only
  from their factories, and moving one copies it, so the handle moved from still works. A consumer
  that fills one in later holds a `std::optional` of it.
- **Axes are found by role, not by position.** A descriptor lists an image's axes in a logical order
  that its schema profile chooses, and a request's ranges follow it; `AxisIndex(axes, role)` finds
  one. The XRADIO profile's order is `l`, `m`, `frequency`, `polarization`, `time`, whatever order
  the store holds them in, but that is the profile's choice rather than a promise. A read returns
  the axes densely packed with axis 0 varying fastest; `AxisDescriptor::storage_index` says where
  each one lives on disk.
- **Pixels come back as `float32`**, converted during the read rather than materialised in their
  stored type first. Reductions accumulate and report in `double`.
- **Pixel masks are applied by default.** An image with a flag reads a flagged pixel as NaN, so a
  consumer that wants the mask reads it as finiteness; `ReadOptions::apply_pixel_mask` turns that
  off. A declared flag must be boolean -- a Zarr `bool`, or the `int8` xarray writes for one -- and
  carry the image's own dimensions in order, or the image does not open.
- **Reads are cancellable.** `ReadOptions::control` carries a cancellation callback, a deadline and
  the cache pool the read keeps what it decodes in -- the context's shared one unless it brings its
  own from `Context::NewCachePool` -- which every read and reduction honours; `ReadOptions` adds a
  ceiling on the temporary memory a request may hold. Progress is an argument of `Image::Read` and of
  `Image::ComputeCubeHistogram`, and a reduction's arrives with its blocks; a request is only ever
  data.
- **Threading**: `Context`, `Dataset` and `Image` handles may be shared and read from several
  threads at once. Reductions running concurrently share one worker pool, which runs one read's
  arithmetic at a time: two reductions interleave read by read rather than running side by side.

## Documentation

- [include/carta-zarr](include/carta-zarr) — the API reference, in the headers.
- [CONTEXT.md](CONTEXT.md) — the vocabulary this codebase uses, and the words it avoids.
- [CHANGELOG.md](CHANGELOG.md) — what each release changed.
- [docs/storage-tuning.md](docs/storage-tuning.md) — choosing a Zarr layout and carta-backend's reader
  settings for Lustre or BeeGFS, what measurements found, and how to measure your own.
- [docs/adr](docs/adr) — decisions and their alternatives.

## Licence

Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA), Associated
Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA).
GPL-3.0-or-later. See [LICENSE](LICENSE).
