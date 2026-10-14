# Your first image read

This example lists the dataset's image candidates, selects an openable image,
and reads up to 64 × 64 pixels from its first channel, polarization and time.
It checks each `Result` before accessing its value and finds spatial axes by
role rather than by position.

The example's source is also its documentation: CI compiles and runs the file
included below against a committed Zarr fixture.

## Build the example

After [installing the library](installation.md), run from the repository root:

```sh
cmake -S examples -B build-examples -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build-examples
./build-examples/carta_zarr_first_read tests/data/images/zarr/xradio/pixels
```

Use a path to your own supported XRADIO dataset instead of the fixture when
ready. The [dataset guide](guides/datasets.md) explains how to produce one.

## Complete source

```{literalinclude} ../examples/first_read.cc
:language: cpp
```

## What each stage does

1. `Context::Create` creates shared cache and execution resources.
2. `Dataset::Open` opens the Zarr hierarchy, recognises the schema and lists
   image candidates. It does not load every image's pixels.
3. `dataset.descriptor().images` exposes a read-only list of `ImageEntry`
   values. Each has an ID, openability, role, axes and diagnostics.
4. `Dataset::OpenImage(id)` describes the chosen image and verifies the
   metadata needed to read it. Coordinate values may be read at this stage.
5. `Image::Read` reads into the caller's `float` buffer. Axis 0 varies fastest
   in the densely packed output; flagged pixels become NaN by default.

Successful `Open` does not mean every candidate is openable. Use the listing
to choose one, and still check the result of `OpenImage` and `Read`: files can
be inaccessible or metadata can be inconsistent.

Next, learn the [core concepts](concepts.md) or [read another channel](guides/reading.md).
For complete statistics, cancellation and animation programs, see
[the runnable examples](examples.md).
