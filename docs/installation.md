# Build and install

## Prerequisites

Use `cmake` 3.24 or newer, a C++17 compiler (GCC 10, Clang 8 or Apple Xcode
11.3.1, or newer, as TensorStore supports), `python` 3.10 or newer, `nasm` and
`patch`. Python is needed only for TensorStore's configure step, NASM because
that step enables the assembler, and `patch` because TensorStore patches some
of its dependencies as it unpacks them.

The first configure downloads TensorStore 0.1.84 and the 42 libraries it
declares, 43 archives of 169 MiB that unpack to about 0.8 GiB of source, so it
needs network access unless you follow the [offline
build](#build-without-network-access) steps. Only ten are compiled --
TensorStore, Abseil, riegeli, re2, zstd, zlib, blosc, snappy, lz4 and
nlohmann_json -- and they are linked into the library statically, with none of
their symbols exported, so they cannot collide with copies your application
links. A release build tree takes about 1.3 GB. The storage benchmark,
`-DCARTA_ZARR_BUILD_BENCH=ON`, also needs [`uv`](https://docs.astral.sh/uv/).

The library is built and tested on macOS 26 on Apple silicon (Apple clang 21,
CMake 4.0, Python 3.13, NASM 3.01), Ubuntu 24.04 (GCC 13.3, CMake 3.28,
Python 3.12, NASM 2.16) and AlmaLinux 8.10 (gcc-toolset-13's GCC 13.3,
CMake 3.26, Python 3.12). CI
builds and tests every change on macOS 15 on Apple silicon, Ubuntu 24.04 and
AlmaLinux 8.

On macOS, install the Xcode command line tools, which provide the Apple clang C++
compiler and the macOS SDK, then:

```sh
brew install cmake nasm python
```

On Ubuntu 24.04, the Python package is named `python3`:

```sh
sudo apt-get install build-essential cmake python3 nasm patch
```

On AlmaLinux 8, the system GCC 8 and Python 3.6 are too old. Install the
gcc-toolset-13 compiler and Python 3.12, enable PowerTools for NASM, and enter
the toolset's shell before building:

```sh
sudo dnf config-manager --set-enabled powertools
sudo dnf install gcc-toolset-13-gcc-c++ cmake python3.12 nasm patch
scl enable gcc-toolset-13 bash
```

Then configure as in [Build from source](#build-from-source), naming Python 3.12
explicitly, since Python 3.6 is also installed.

## Build from source

```sh
git clone https://github.com/CARTAvis/carta-zarr.git
cd carta-zarr

# configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# build
cmake --build build --parallel
# install
cmake --install build --prefix "$HOME/.local"
```

`--parallel` builds several source files at once. Without a number, it uses all the
cores on the machine; give a number, such as `--parallel 4`, to leave some free.

The first build takes substantially longer than later builds, because it also
builds TensorStore.

On AlmaLinux 8, configure in the gcc-toolset shell with Python 3.12 named
explicitly, instead of the configure line above:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPython3_EXECUTABLE=/usr/bin/python3.12
```

To install into a system prefix such as `/opt/carta-zarr`, run
`sudo cmake --install build --prefix /opt/carta-zarr`, then pass
`-DCMAKE_PREFIX_PATH=/opt/carta-zarr` when building against it.

## Build without network access

Some build systems, such as Launchpad, have no network access during builds.
Configure once on a machine that has network access, pack the archives that the
configure downloaded, copy the archive to the offline machine, and point the
offline configure at it:

```sh
# Where there is network, after a configure into build/; writes the archive to the repository root:
(cd build/_deps && find . -maxdepth 4 -path '*-subbuild/*-populate-prefix/src/*' -type f \
    | tar czf ../../carta-zarr-deps.tar.gz -T -)

# Where there is not:
mkdir deps && tar xzf carta-zarr-deps.tar.gz -C deps
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFETCHCONTENT_BASE_DIR=$PWD/deps
```

Each archive is checked against the hash TensorStore pins and then not
downloaded again. Move the archives, not the unpacked sources: TensorStore
patches each source as it unpacks it, and the patches name the build tree they
were made for. Copying another tree's `_deps`, or setting
`FETCHCONTENT_SOURCE_DIR_<name>`, does not give a tree that configures.

The archives match the dependency versions pinned in this version of the
repository. After a change that updates them, configure again on a machine with
network access and pack a new archive.

## Run the tests (optional)

The tests are for checking a build on a new platform or a change to the
library; installing does not need them. A build from source compiles them by
default; run them with:

```sh
ctest --test-dir build --output-on-failure
```

The tests read the fixtures under `tests/data`, which a clone includes.

A change to how the library shares anything between threads is worth running
under ThreadSanitizer as well. It needs a build tree of its own, because
TensorStore has to be built with it too; `CARTA_ZARR_THREAD_SANITIZER` does
that:

```sh
cmake -S . -B build-tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCARTA_ZARR_THREAD_SANITIZER=ON
cmake --build build-tsan -j 8
ctest --test-dir build-tsan --output-on-failure -R 'concurrency|pixel-read|work-pool|memo|read-ahead'
```

These are the tests that share handles between threads. The pool-release and
allocation-failure tests count allocations by replacing `operator new`, which
ThreadSanitizer's own allocator replaces in turn, so they fail under it whatever
the threads do.

## Link from another CMake project

```cmake
find_package(CartaZarr 0.3 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

Configure your application with `-DCMAKE_PREFIX_PATH="$HOME/.local"` when that
prefix is not otherwise searched. Public headers require only the C++ standard
library; TensorStore is a private library dependency.

## Versions and distribution

Before 1.0, `find_package(CartaZarr 0.3)` accepts only 0.3.x: a different minor
version, such as 0.4, is not compatible. Record the exact library version your
application uses, and rebuild when you change minor versions. Binary distributions must carry `LICENSE` and
`THIRD_PARTY_NOTICES`, which the build writes from the sources of the ten
libraries linked in and of the header-only `half`, and installs beside
`LICENSE`.

The current development version is 0.3.0. Applications built for 0.2 must update
their CMake version requirement and rebuild against 0.3, because public descriptor
changes affect the ABI (the binary interface). The [changelog](changelog.md) lists
the changes to check in your application.

Continue with [your first image read](quickstart.md).
