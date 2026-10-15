# Build and install

## Prerequisites

Use CMake 3.24 or newer, a C++17 compiler, Python 3.10 or newer, NASM and `patch`.
The first configure downloads TensorStore and its pinned dependencies. Python
is needed for TensorStore's configure step; Python packages are not required
to use the C++ library.

On macOS, install the Xcode command line tools, then:

```sh
brew install cmake nasm python
```

On Ubuntu 24.04:

```sh
sudo apt-get install build-essential cmake python3 nasm patch
```

For AlmaLinux 8 and offline dependency archives, see the repository's
[build instructions](https://github.com/CARTAvis/carta-zarr/blob/dev/README.md#requirements).

## Build and install into a prefix you own

```sh
git clone https://github.com/CARTAvis/carta-zarr.git
cd carta-zarr
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
cmake --install build --prefix "$HOME/.local"
```

The README installs into `/opt/carta-zarr` with `sudo` instead; use whichever
prefix suits you, and pass it as `CMAKE_PREFIX_PATH` when building against it.
Choose a parallel job count appropriate to your machine. TensorStore makes the
first build substantially larger than subsequent builds.

## Run the tests (optional)

The tests are for checking a build on a new platform or a change to the
library; installing does not need them. The build above already compiled them;
run them with:

```sh
ctest --test-dir build --output-on-failure
```

The tests read the fixtures under `tests/data`, which a clone includes.

## Link from another CMake project

```cmake
find_package(CartaZarr 0.2 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

Configure your application with `-DCMAKE_PREFIX_PATH="$HOME/.local"` when that
prefix is not otherwise searched. Public headers require only the C++ standard
library; TensorStore is a private library dependency.

Before 1.0, `0.2` accepts 0.2.x rather than another minor version. Record the
exact library version used by your application and rebuild when changing minor
versions. Binary distributions must carry `LICENSE` and the installed
`THIRD_PARTY_NOTICES`.

The current development version is 0.2.0. Applications built for 0.1 must update
their CMake version requirement and rebuild against 0.2: public descriptor
changes affect the ABI. See [upgrading to 0.2](upgrading/0.2.md)
for the metadata changes to check in your application.

Continue with [your first image read](quickstart.md).
