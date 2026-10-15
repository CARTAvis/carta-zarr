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

Replace the fixture path with your own supported dataset when ready.

| Example | What it demonstrates | Explanation |
|---|---|---|
| `first_read` | Discover candidates, find axes by role and read a bounded spatial patch | [Quickstart](quickstart.md) |
| `statistics` | Reduce two regions in one pass, consume final blocks and compute plane/cube histograms | [Statistics](guides/statistics.md) |
| `read_control` | Select a channel, report progress, cancel, give a scan its own cache pool and prefetch | [Read controls](guides/resources.md#complete-read-control-example) |
| `animation` | Read ahead during an animation with `ReadAhead`: read frames, report upcoming ones and cancel | [Animation](guides/resources.md#complete-animation-example) |

Each linked page explains what its program does and its limits. All of them
bound the destination to a 64 × 64 patch, but decoding still operates on whole
chunks.

## Shared setup

The advanced programs use `support.h` to open the first openable image, validate
required axes, create a bounded patch and calculate its destination size. The
quickstart remains self-contained so it can be copied into a new application.
The `Image` the helper returns keeps the dataset's resources alive after the
helper's `Dataset` goes out of scope.

```{literalinclude} ../examples/support.h
:language: cpp
```

CI builds and runs every example; see
[contributing to the documentation](contributing-docs.md#sources-and-checks).
