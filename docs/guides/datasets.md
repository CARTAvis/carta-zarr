# Prepare an XRADIO dataset

The `xradio.image` profile supports image schema `0.0.2`, written by XRADIO 1.2.5,
and older unversioned stores. Conformance fixtures from both XRADIO 1.2.5 and
1.2.3 are committed and tested. Schema version, XRADIO package version and
carta-zarr library version are separate values; see
[schema compatibility](compatibility.md) before interpreting them.

XRADIO converts a FITS or CASA image:

```python
from xradio.image import open_image, write_image

xds = open_image("cube.fits")
write_image(xds, "cube", out_format="zarr")  # XRADIO 1.2.5 writes cube.img.zarr
```

Open the resulting `cube.img.zarr` with carta-zarr. Pin XRADIO when producing
datasets for a reproducible workflow; the example above follows 1.2.5.

XRADIO's FITS reader requires complete coordinate metadata. The repository's
[pinned fixture generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/generate_conformance_fixtures.py)
is a reproducible example of producing conforming data, pinned to XRADIO 1.2.5.
Its flagged fixture uses CASA images and requires casatools; this is separate
from the dependencies needed to build the C++ reader or the documentation site.

Use `ProbeSchema(path, kXradioImageSchema)` to distinguish `match`, `no_match`
and `invalid`. A successful `Result` alone is not evidence of a match: inspect
its `kind`. Opening the dataset then provides the list of image candidates.
The returned `schema_version` is the version declared by the store, or empty
when absent. Check diagnostics too: a store declaring 0.0.2 may be read using
the unversioned rules if its required version-specific fields are missing.

For layout experiments, the
[benchmark generator](https://github.com/CARTAvis/carta-zarr/blob/dev/tools/zarr-bench/README.md)
can rewrite a real cube or synthesise one. Consult [storage tuning](../storage-tuning.md)
before choosing chunk and shard sizes for a particular workload.
