# Schema versions and reading rules

Three version numbers describe different things:

| Version | Meaning |
|---|---|
| carta-zarr 0.2.0 | The C++ library and its API/ABI |
| XRADIO 1.2.5 | The Python package used to write the current conformance fixtures |
| Image schema 0.0.2 | The layout declared in the dataset's root `schema_version` attribute |

The schema profile that reads it is `kXradioImageSchema` (`xradio.image`). Do not compare the dataset's schema version to the installed
C++ library's version.

## How a store is read

carta-zarr reads XRADIO image schema 0.0.2 only, which XRADIO 1.2.5 and later
write.

| Dataset metadata | Result | First diagnostic, error code |
|---|---|---|
| Declares 0.0.2 and meets the 0.0.2 checks | Read | None |
| Declares 0.0.2 but fails the 0.0.2 checks | Refused | `invalid_metadata`, with the reason |
| An image dataset that declares no version, as XRADIO 1.2.4 and earlier wrote | Refused | `unsupported_schema_version` |
| An image dataset that declares another version, an empty string or a non-string value | Refused | `unsupported_schema_version` |

The 0.0.2 checks require data groups naming a sky array typed `sky` and, when
the frequency coordinate is readable, its `frame` and `units`. They apply even
to a store that declares 0.0.2: XRADIO 1.2.5 stamps 0.0.2 on an older store it
upgrades on read, whether or not the upgrade reached those fields. Ordinary
image and coordinate validation applies after them.

For a refused image dataset, `ProbeSchema` returns `kind == invalid` with the
reason as its first diagnostic, and `Dataset::Open` returns an error with the
code above. `unsupported_schema_version` means a different build, or a store
written again with XRADIO 1.2.5 or later, may open it; `invalid_metadata` means
the store is malformed. A non-image Zarr dataset is a non-match whether or not
it has a root `schema_version`.

`SchemaProbeResult::schema_version` reports what the store declares: `"0.0.2"`
for every store that opens, and empty for one that declares nothing.

## Metadata differences that matter to applications

- **Image roles:** a variable typed `sky` that a data group other than `base`
  names as its sky reports that group's name as its role, such as `model`,
  `residual` or `dirty`. Other images take their role from their `type`, such
  as `point_spread_function` or `mask`. Use the reported role instead of
  inferring it from the ID.
- **Default image:** the listing puts XRADIO's own image names first (`SKY`,
  `SKY_MODEL`, `SKY_RESIDUAL`, `SKY_DIRTY`, …), and the default is the first
  openable one, so a dataset without `SKY` defaults to, for example,
  `SKY_MODEL`.
- **Image subtype:** `ImageDescriptor::sub_type` describes what is measured,
  for example `Intensity`, `SpectralIndex` or `OpticalDepth`. It is separate
  from `image_role`, and is empty when the image does not name one.
  This field is on the opened image descriptor, not `ImageEntry`.
- **Spectral frame and units:** read from `frequency`'s `frame` and `units`.
  `SpectralCoordinate::system` is the frame in casacore's names, as 0.0.2
  writes it, such as `LSRK` or `GEO`. The reference frequency's `observer` is
  not read.
- **Single-channel images:** a finite, nonzero `channel_width`
  can supply `SpectralCoordinate::increment`, with reference pixel 1 and the
  single frequency as reference value. The width's unit must be absent or match
  the spectral unit. Multi-channel increments still come from their samples;
  a missing or invalid width does not invent a spacing.
- **Observation dates:** use `ObservationInfo::mjd_obs` (Modified Julian Date, MJD) and `timesys`.
  `date_obs` has been removed; string-valued observation dates and numeric
  equinox forms not written by XRADIO are not interpreted, and neither is a
  numeric date without a format.

- **Integer sky images:** XRADIO's FITS reader keeps an integer FITS image's
  type, although the schema allows only float32 and float64. Such an image
  opens, `ImageDescriptor::stored_type` reports its type, and its pixels are
  read as float32 like any other.

The tests cover XRADIO's own 0.0.2 conversions, including a single-channel
continuum image. See [ADR 0023 (architecture decision record)](https://github.com/CARTAvis/carta-zarr/blob/dev/docs/adr/0023-only-xradio-image-schema-0-0-2-is-read.md)
for the compatibility policy and the
[fixture notes](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/README.md)
for how these datasets were produced.

To move an application from carta-zarr 0.1, see the [changelog](../changelog.md).
