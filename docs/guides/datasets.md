# Prepare an XRADIO dataset

The supported profile was checked against XRADIO 1.2.x, with a conformance fixture
written by XRADIO 1.2.3. The profile's `xradio.image` version `1.2` is this
library's compatibility label, not a schema version XRADIO publishes.

XRADIO converts a FITS or CASA image:

```python
from xradio.image import open_image, write_image

xds = open_image("cube.fits")
write_image(xds, "cube.zarr", out_format="zarr")
```

XRADIO's FITS reader requires complete coordinate metadata. Its casacore layer
is imported even for FITS; on macOS that involves casatools. The repository's
[pinned fixture generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/generate_conformance_fixtures.py)
is a reproducible example of producing conforming data.

Use `ProbeSchema(path, kXradioImageSchema)` to distinguish `match`, `no_match`
and `invalid`. A successful `Result` alone is not evidence of a match: inspect
its `kind`. Opening the dataset then provides the list of image candidates.

For layout experiments, the
[benchmark generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tools/zarr-bench/README.md)
can rewrite a real cube or synthesise one. Consult [storage tuning](../storage-tuning.md)
before choosing chunk and shard sizes for a particular workload.
