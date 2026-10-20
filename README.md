# carta-zarr

[![CI](https://github.com/CARTAvis/carta-zarr/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/CARTAvis/carta-zarr/actions/workflows/ci.yml)
[![Documentation](https://github.com/CARTAvis/carta-zarr/actions/workflows/docs.yml/badge.svg?branch=dev)](https://cartavis.org/carta-zarr/)
[![Release](https://img.shields.io/github/v/release/CARTAvis/carta-zarr)](https://github.com/CARTAvis/carta-zarr/releases)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)](LICENSE)

> [!WARNING]
> carta-zarr is still under heavy development, and its API may change often and without notice.
> Before 1.0 any minor release may break source and binary compatibility
> ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)), so build against an exact
> version and expect to change your code when you move to another.

A read-only C++17 library for [XRADIO](https://github.com/casangi/xradio) images stored as Zarr v3.
It is written as the Zarr reader of [carta-backend](https://github.com/CARTAvis/carta-backend),
CARTA's image server, and kept apart from it so that CARTA and other C++ consumers can read these
datasets without depending on TensorStore, casacore or CARTA's protobuf. TensorStore is linked in
privately; the public headers need only the C++ standard library.

**[Documentation](https://cartavis.org/carta-zarr/)** -- installation, a first read, concepts, task
guides and the API reference.

## What it reads

- XRADIO image schema 0.0.2, as XRADIO 1.2.5 and later write it, from a local filesystem.
- Real-valued images on the sky plane (`l`, `m`) or the aperture plane (`u`, `v`), with
  `frequency`, `polarization` and `time`.
- Pixels as float32, with the image's flag applied as NaN by default; statistics accumulated in
  double, and histogram counts as integers.

What is refused, and why, is in
[schema versions](https://cartavis.org/carta-zarr/reference/schema-versions.html).

## Build and use

CMake 3.24, a C++17 compiler, Python 3.10, NASM and `patch`; the first configure downloads
TensorStore. The [installation guide](https://cartavis.org/carta-zarr/installation.html) has the details, the platforms it is
tested on, an offline build, and the ThreadSanitizer build.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
cmake --install build --prefix "$HOME/.local"
```

```cmake
find_package(CartaZarr 0.2 CONFIG REQUIRED)  # before 1.0, finds 0.2.x and nothing else
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

```c++
#include <carta-zarr/carta_zarr.h>

auto context = carta::zarr::Context::Create();
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

Every call returns a `Result` to check before use; the [first read](https://cartavis.org/carta-zarr/quickstart.html) is the same
example, complete and with its checks.

## More

- [include/carta-zarr](include/carta-zarr) -- every type, field and call, documented where it is
  declared.
- [GLOSSARY.md](GLOSSARY.md) -- the vocabulary this codebase uses, and the words it avoids.
- [CHANGELOG.md](CHANGELOG.md) -- what each release changed.
- [Storage tuning](https://cartavis.org/carta-zarr/storage-tuning.html) -- choosing a Zarr layout and carta-backend's
  reader settings for Lustre or BeeGFS.
- [docs/adr](docs/adr) -- decisions and their alternatives.
- [docs/contributing-docs.md](docs/contributing-docs.md) -- building the documentation locally.

## Licence

Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA), Associated
Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA).
GPL-3.0-or-later; see [LICENSE](LICENSE). A binary distribution also carries the
`THIRD_PARTY_NOTICES` the build installs.
