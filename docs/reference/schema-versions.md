# Schema versions and reading rules

Three version numbers describe different things:

| Version | Meaning |
|---|---|
| carta-zarr 0.2.0 | The C++ library and its API/ABI |
| XRADIO 1.2.5 | The Python package used to write the current conformance fixtures |
| Image schema 0.0.2 | The layout declared in the dataset's root `schema_version` attribute |

Every version is read by the same schema profile, `kXradioImageSchema`
(`xradio.image`). Do not compare the dataset's schema version to the installed
C++ library's version.

## How a store is read

| Dataset metadata | Reading rules | Reported diagnostic |
|---|---|---|
| Declares 0.0.2 and meets the version-specific checks | Schema 0.0.2 | No fallback diagnostic |
| Declares no schema version | Unversioned, as XRADIO 1.2.4 and earlier wrote | `unversioned_schema` |
| Declares 0.0.2 but fails the version-specific checks | Unversioned, as XRADIO 1.2.4 and earlier wrote | `unversioned_schema`, with the reason |
| An image dataset declares another version, an empty string or a non-string value | Refused | `unsupported_schema_version` |

The 0.0.2 checks require data groups naming a sky array typed `sky` and, when
the frequency coordinate is readable, its `frame` and `units`. These select
the version-specific rules; ordinary image and coordinate validation still
applies: falling back to the unversioned rules does not make malformed data
readable.

`SchemaProbeResult::schema_version` and `DatasetDescriptor::schema_version`
report **what the store declares**, even when 0.0.2 falls back to older rules.
They are empty for a store with no declared version. Inspect the diagnostics
to distinguish the fallback; do not infer the reading rules from the version
string alone.

For an unsupported image schema, `ProbeSchema` returns `kind == invalid` with
an `unsupported_schema_version` diagnostic. `Dataset::Open` returns an error
with that code. A non-image Zarr dataset is still a non-match even if it happens
to contain a root `schema_version`.

## Metadata differences that matter to applications

- **Image roles:** a variable typed `sky` that a data group other than `base`
  names as its sky reports that group's name as its role, such as `model`,
  `residual` or `dirty`. Other images take their role from their `type`. Use
  the reported role instead of inferring it from the ID.
- **Default image:** the listing puts XRADIO's own image names first (`SKY`,
  `SKY_MODEL`, `SKY_RESIDUAL`, `SKY_DIRTY`, …), and the default is the first
  openable one, so a dataset without `SKY` defaults to, for example,
  `SKY_MODEL`.
- **Image subtype:** `ImageDescriptor::sub_type` describes what is measured,
  for example `Intensity`, `SpectralIndex` or `OpticalDepth`. It is separate
  from `image_role`, is read under 0.0.2 rules and is empty when not available.
  This field is on the opened image descriptor, not `ImageEntry`.
- **Spectral frame and units:** 0.0.2 reads `frequency`'s `frame` and `units`.
  The unversioned rules use the reference frequency's observer and can use its
  unit when the coordinate has none. `SpectralCoordinate::system` uses names
  suitable for casacore/FITS: the older observer `gcrs` is reported as `GEO`.
- **Single-channel images:** under 0.0.2 rules, a finite, nonzero `channel_width`
  can supply `SpectralCoordinate::increment`, with reference pixel 1 and the
  single frequency as reference value. The width's unit must be absent or match
  the spectral unit. Multi-channel increments still come from their samples;
  a missing or invalid width does not invent a spacing.
- **Observation dates:** use `ObservationInfo::mjd_obs` (Modified Julian Date, MJD) and `timesys`.
  `date_obs` has been removed; string-valued observation dates and numeric
  equinox forms not written by XRADIO are no longer interpreted. A numeric
  date without a format is treated as MJD only under the unversioned rules.

The tests cover both schema families and a 0.0.2 single-channel continuum image.
See [ADR 0022 (architecture decision record)](https://github.com/CARTAvis/carta-zarr/blob/dev/docs/adr/0022-a-store-is-read-by-the-schema-version-it-names.md)
for the compatibility policy and the
[fixture notes](https://github.com/CARTAvis/carta-zarr/blob/dev/tests/data/README.md)
for how these datasets were produced.

To move an application from carta-zarr 0.1, see the [changelog](../changelog.md).
