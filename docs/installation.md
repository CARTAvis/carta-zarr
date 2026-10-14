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

## Build, test and install into a prefix you own

```sh
git clone https://github.com/CARTAvis/carta-zarr.git
cd carta-zarr
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build --output-on-failure
cmake --install build --prefix "$HOME/.local"
```

Choose a parallel job count appropriate to your machine. TensorStore makes the
first build substantially larger than subsequent builds.

## Link from another CMake project

```cmake
find_package(CartaZarr 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

Configure your application with `-DCMAKE_PREFIX_PATH="$HOME/.local"` when that
prefix is not otherwise searched. Public headers require only the C++ standard
library; TensorStore is a private library dependency.

Before 1.0, `0.1` accepts 0.1.x rather than another minor version. Record the
exact library version used by your application and rebuild when changing minor
versions. Binary distributions must carry `LICENSE` and the installed
`THIRD_PARTY_NOTICES`.

Continue with [your first image read](quickstart.md).
