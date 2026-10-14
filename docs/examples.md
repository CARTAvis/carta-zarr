# Run the C++ examples

These are complete command-line applications linked against the installed
public C++ API. They use a supported XRADIO dataset passed as their only argument.
Each handles errors and returns a nonzero exit status on failure.

## Build and run

After [installing the library](installation.md), run from the repository root:

```sh
cmake -S examples -B build-examples -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build-examples

./build-examples/carta_zarr_first_read tests/data/images/zarr/xradio/pixels
./build-examples/carta_zarr_statistics tests/data/images/zarr/xradio/pixels
./build-examples/carta_zarr_read_control tests/data/images/zarr/xradio/pixels
./build-examples/carta_zarr_animation tests/data/images/zarr/xradio/pixels
```

Replace the fixture path with your own supported dataset when ready. The C++
examples require the library; the [documentation website build](contributing-docs.md)
requires only documentation tools.

| Example | What it demonstrates | Source and explanation |
|---|---|---|
| `first_read` | Discover candidates, find axes by role and read a bounded spatial patch | [Quickstart](quickstart.md) |
| `statistics` | Reduce two regions in one pass, consume final blocks and compute plane/cube histograms | [Statistics](guides/statistics.md) |
| `read_control` | Select a channel, use progress and cancellation, isolate a scan's cache and prefetch | [Read controls](guides/resources.md#complete-read-control-example) |
| `animation` | Create a ReadAhead controller, read frames, report upcoming frames and cancel | [Animation](guides/resources.md#complete-animation-example) |

The statistics example examines at most four channels at the first polarization
and time. Its region reduction uses a patch up to 64 × 64 pixels; histograms
cover **whole planes**. Reading and animation keep the destination bounded to a
64 × 64 patch, but decoding still operates on whole chunks. The animation cache
is 64 MiB: on large datasets, creation can fail with the required cache capacity
in its error message. The animation CLI demonstrates API usage, not frame-rate
or read-ahead performance.

## Shared setup

The advanced programs use `support.h` to open the first openable image, validate
required axes, create a bounded patch and calculate its destination size. The
quickstart remains self-contained so it can be copied into a new application.
An image retains its storage resources after the helper's Dataset goes out of
scope.

```{literalinclude} ../examples/support.h
:language: cpp
```

## CI validation

The existing package-consumer test compiles all examples using the installed
`CartaZarr` CMake package. It runs each on the committed synthetic pixel fixture,
XRADIO 1.2.5's 0.0.2 cube and continuum fixtures, and XRADIO 1.2.3's unversioned
cube. It also checks that missing arguments report failure. The statistics consumer checks
that every requested channel receives a final block; the read-control consumer
checks final progress and the expected `cancelled` error. No assertion depends
on background prefetch timing.

To run this check after building the library:

```sh
ctest --test-dir build -R '^carta-zarr-package-consumer$' --output-on-failure
```
