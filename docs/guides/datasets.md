# Prepare an XRADIO dataset

carta-zarr reads datasets written by XRADIO. This page covers producing one and
checking it; which schema versions the library accepts, and how the schema,
XRADIO package and library versions differ, is set out in
[schema versions](../reference/schema-versions.md).

## Write a dataset

Install the XRADIO version you mean to produce datasets with, and keep it
pinned for a reproducible workflow; this page follows 1.2.5:

```sh
pip install "xradio==1.2.5"
```

XRADIO converts a FITS or CASA image:

```python
from xradio.image import open_image, write_image

xds = open_image("cube.fits")
write_image(xds, "cube", out_format="zarr")  # writes cube.img.zarr
```

Open the resulting `cube.img.zarr` with carta-zarr. XRADIO's FITS reader needs
both celestial axes in the header and refuses axis types it does not support.

The repository's
[fixture generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/generate_conformance_fixtures.py)
is a reproducible example of producing conforming data with XRADIO 1.2.5. Its
fixture with a flag is converted from a CASA image, so it also needs casatools.

## Check the dataset

`ProbeSchema(path, kXradioImageSchema)` tells whether a path holds a dataset
carta-zarr reads. A successful `Result` alone is not a match: inspect its `kind`,
which is `match`, `no_match` or `invalid`, and its diagnostics.
`schema_version` is the version the dataset declares, or empty. Only 0.0.2 is
read: a dataset that declares none, as XRADIO 1.2.4 and earlier wrote, is
refused, and has to be written again with XRADIO 1.2.5 or later. Opening the
dataset then lists its image candidates. What is refused, and why, is in
[schema versions](../reference/schema-versions.md).

## Next steps

- [List images and read a channel](reading.md) from the dataset.
- Rewrite a real cube or synthesise one in other layouts with the
  [benchmark generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tools/zarr-bench/README.md).
- Choose chunk and shard sizes for a workload with [storage tuning](../storage-tuning.md).
