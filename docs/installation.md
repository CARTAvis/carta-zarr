# Build and install

## Prerequisites

Use `cmake` 3.24 or newer, a C++17 compiler, `python` 3.10 or newer, `nasm` and `patch`.
The first configure downloads TensorStore and its pinned dependencies, so it needs
network access unless you follow the [offline build](#build-without-network-access)
steps. Python is needed for TensorStore's configure step.

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

## Run the tests (optional)

The tests are for checking a build on a new platform or a change to the
library; installing does not need them. The build above already compiled them;
run them with:

```sh
ctest --test-dir build --output-on-failure
```

The tests read the fixtures under `tests/data`, which a clone includes.

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

## Link from another CMake project

```cmake
find_package(CartaZarr 0.2 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE CARTA::zarr)
```

Configure your application with `-DCMAKE_PREFIX_PATH="$HOME/.local"` when that
prefix is not otherwise searched. Public headers require only the C++ standard
library; TensorStore is a private library dependency.

## Versions and distribution

Before 1.0, `find_package(CartaZarr 0.2)` accepts only 0.2.x: a different minor
version, such as 0.3, is not compatible. Record the exact library version your
application uses, and rebuild when you change minor versions. Binary distributions must carry `LICENSE` and the installed
`THIRD_PARTY_NOTICES`.

The current development version is 0.2.0. Applications built for 0.1 must update
their CMake version requirement and rebuild against 0.2, because public descriptor
changes affect the ABI (the binary interface). The [changelog](changelog.md) lists
the changes to check in your application.

Continue with [your first image read](quickstart.md).
