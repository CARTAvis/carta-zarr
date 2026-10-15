# Your first image read

This example lists the dataset's image candidates, selects an openable image,
and reads up to 64 × 64 pixels from its first channel, polarization and time.
It checks each `Result` before accessing its value and finds spatial axes by
role rather than by position.

The source below is the complete example. CI builds it and runs it against a
committed Zarr fixture, so the code on this page is tested.

## Build the example

After [installing the library](installation.md), run from the repository root:

```sh
cmake -S examples -B build-examples -DCMAKE_PREFIX_PATH="$HOME/.local"
cmake --build build-examples
./build-examples/carta_zarr_first_read tests/data/images/zarr/xradio/pixels
```

The output for that fixture looks like this:

```text
schema=xradio.image version=0.0.2
SKY openable=1 role=sky
spectral frame=LSRK increment=1e+06 Hz
Read 20 pixels from SKY; 13 are finite
```

To read your own XRADIO dataset, pass its path instead of the fixture. It must
be XRADIO image schema 0.0.2, which XRADIO 1.2.5 and later write; an older
store is refused with `unsupported_schema_version`. The [dataset guide](guides/datasets.md)
explains how to produce one.

## Complete source

```{literalinclude} ../examples/first_read.cc
:language: cpp
```

## What each stage does

1. `Context::Create` creates shared cache and execution resources.
2. `Dataset::Open` opens the Zarr hierarchy, recognises the schema and lists
   image candidates. The example prints its declared schema version and
   diagnostics. It does not load every
   image's pixels.
3. `dataset.descriptor().images` exposes a read-only list of `ImageEntry`
   values. Each has an ID, openability, role, axes and diagnostics.
4. `Dataset::OpenImage(id)` describes the chosen image and verifies the
   metadata needed to read it. Coordinate values may be read at this stage.
   The example prints the image's subtype when it has one, its spectral frame,
   and its increment when there is one; [schema versions](reference/schema-versions.md#metadata-differences-that-matter-to-applications)
   says where each comes from.
5. `Image::Read` reads into the caller's `float` buffer, with axis 0 varying
   fastest (see [axes and output order](concepts.md#axes-and-output-order));
   flagged pixels become NaN by default.

Successful `Open` does not mean every candidate is openable. Use the listing
to choose one, and still check the result of `OpenImage` and `Read`: files can
be inaccessible or metadata can be inconsistent.

Next, learn the [core concepts](concepts.md) or [read another channel](guides/reading.md).
For complete statistics, cancellation and animation programs, see
[the runnable examples](examples.md).
