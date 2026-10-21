# Prepare an XRADIO dataset

carta-zarr reads datasets written by XRADIO. This page covers producing one and
checking it; which schema versions the library accepts, and how the schema,
XRADIO package and library versions differ, is set out in
[schema versions](../reference/schema-versions.md).

## Write a dataset

Install the XRADIO version you mean to produce datasets with, and keep it
pinned for a reproducible workflow; this page follows 1.3.0. The `casacore`
extra is needed to read CASA images:

```sh
pip install "xradio[zarr,casacore]==1.3.0"
```

XRADIO converts a FITS or CASA image:

```python
from xradio.image import open_image, write_image

xds = open_image("cube.fits")
write_image(xds, "cube", out_format="zarr")  # writes cube.img.zarr
```

Open the resulting `cube.img.zarr` with carta-zarr. XRADIO's FITS reader
needs a `CTYPE` on every axis and refuses axis types it does not support,
such as `LINEAR`. Without `DATE-OBS` or `MJD-OBS`, the image's time is MJD 0.

The repository's
[fixture generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/generate_conformance_fixtures.py)
is a reproducible example of producing conforming data with XRADIO 1.3.0. Its
fixtures with a flag are converted from CASA images, so it also needs casatools.

## Check the dataset

`ProbeSchema(path, kXradioImageSchema)` tells whether a path holds a dataset
carta-zarr reads. A successful `Result` alone is not a match: inspect its `kind`,
which is `match`, `no_match` or `invalid`, and its diagnostics.
`schema_version` is the version the dataset declares, or empty. Only 0.0.2 is
read: a dataset that declares none, as XRADIO 1.2.4 and earlier wrote, is
refused, and has to be written again with XRADIO 1.2.5 or later. Opening the
dataset then lists its image candidates. What is refused, and why, is in
[schema versions](../reference/schema-versions.md).

`Dataset::Size(timeout)` reports how large the dataset is, and which of two
answers it gives: `SizeBasis::measured`, the files under the store added up, or
`SizeBasis::declared`, every array's shape times its element size, uncompressed,
when the store cannot be measured: the timeout passed part-way through, a
directory could not be listed, or the total overflowed. The second is not a
bound on the first, since compression makes a store smaller and metadata makes
it larger, so show which one you got.

## Next steps

- [List images and read a channel](reading.md) from the dataset.
- Rewrite a real cube or synthesize one in other layouts with the
  [benchmark generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tools/zarr-bench/README.md).
- Choose chunk and shard sizes for a workload with [storage tuning](../storage-tuning.md).
