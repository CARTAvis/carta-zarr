# carta-zarr

A read-only C++17 library for XRADIO images stored as Zarr v3, extracted from
[carta-backend](https://github.com/CARTAvis/carta-backend) so that CARTA and other C++ consumers can
read these datasets without depending on TensorStore, casacore, or CARTA protobuf.

TensorStore is a private implementation detail: it is linked into this library and never appears in
its public headers, which need only the C++ standard library.

## What it reads

- **Image datasets** written by [XRADIO](https://github.com/casangi/xradio) 1.2.x: one Zarr v3 group
  holding sky-plane image variables and the coordinates they share. Schema id `xradio.image`,
  schema version `1.2`.
- **Images**: real-valued variables carrying all five of `time`, `frequency`, `polarization`, `l`
  and `m`. Complex and aperture-plane (`u`, `v`) variables are listed with a diagnostic saying why
  they cannot be opened.
- **Local filesystem stores.** Metadata may be per-node or consolidated into the root.

Writing is out of scope, and so is anything but the `xradio.image` profile.

## Requirements

CMake 3.24, a C++17 compiler, and `nlohmann_json` 3.11 or newer. TensorStore 0.1.84 is fetched and
built as part of the build, and brings its own zlib, zstd, blosc and protobuf.

## Build, test, install

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

The first configure builds TensorStore, which takes a while; later ones do not.

## Using it from CMake

```cmake
find_package(CartaZarr 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

`tests/consumer` is exactly this, built against the install tree by the `carta-zarr-package-consumer`
test, so the packaging stays consumable.

## Using the API

```c++
#include <carta-zarr/carta_zarr.h>

auto context = carta::zarr::Context::Create();          // one per process is enough
auto dataset = carta::zarr::Dataset::Open(context.value(), "/path/to/image.zarr");
auto image = dataset.value().OpenImage(*dataset.value().descriptor().default_image_id);

carta::zarr::ReadRequest request;                        // one Range per logical axis
request.axes = {{0, 512, 1}, {0, 512, 1}, {0, 1, 1}, {0, 1, 1}, {0, 1, 1}};
std::vector<float> pixels(512 * 512);
auto written = image.value().Read(request, {pixels.data(), pixels.size() * sizeof(float)});
```

`carta::zarr::IsXradioImage(path)` answers the "can this library open it?" question on its own, and
`Probe` reports what a path is without opening anything.

### What the API promises

- **Everything is a `Result<T>`.** A failure carries an `ErrorCode`, a message and the node it came
  from. No public entry point throws: a malformed value in the metadata or a length no allocator can
  serve is reported as an error like any other.
- **Logical axis order is `l`, `m`, `frequency`, `polarization`, `time`** (`kXradioImageAxisOrder`),
  whatever order the store holds them in. A read returns them densely packed with axis 0 varying
  fastest; `AxisDescriptor::storage_index` says where each one lives on disk.
- **Pixels come back as `float32`**, converted during the read rather than materialised in their
  stored type first. Reductions accumulate and report in `double`.
- **Pixel masks are applied by default.** An image with a flag variable reads a masked pixel as NaN;
  `ReadOptions::apply_pixel_mask` turns that off, and `ReadPixelMask` returns the mask itself as one
  byte per pixel. A declared mask must be boolean and carry the image's own dimensions in order, or
  the image does not open.
- **Reads are cancellable.** `ReadOptions` carries a cancellation callback, a deadline, a progress
  callback, and a ceiling on the temporary memory a request may hold.
- **Threading**: `Context`, `Dataset` and `Image` handles may be shared and read from several
  threads at once. Reductions running concurrently share one worker pool and are serialised inside
  it, so they are safe but do not run in parallel with each other.

## Documentation

- [CONTEXT.md](CONTEXT.md) — the vocabulary this codebase uses, and the words it avoids.
- [docs/design.md](docs/design.md) — why the library is split from the backend the way it is.
- [docs/read-api.md](docs/read-api.md) — the read and reduction API, and the reasoning behind its
  output types and streaming shape.
- [docs/adr](docs/adr) — decisions and their alternatives.
- [AGENTS.md](AGENTS.md) — building and testing without drowning in TensorStore's output.

## Licence

GPL-3.0-or-later. See [LICENSE](LICENSE).
