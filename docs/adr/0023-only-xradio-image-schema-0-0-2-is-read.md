# Only XRADIO image schema 0.0.2 is read

*Supersedes [ADR 0022](0022-a-store-is-read-by-the-schema-version-it-names.md)'s rows: its checks,
its refusal of an unknown version and what the schema is held to stand.*

ADR 0022 read a store by one of two rows: schema 0.0.2's, and one for stores that name no version or
name 0.0.2 without being one, read as XRADIO 1.2 wrote them. The second row was meant to be
transitional. It kept every pre-0.0.2 form alive: the spectral frame from the reference frequency's
observer, with astropy's `gcrs` read as casacore's `GEO`; a unit from the reference frequency; an
observation date with no format read as an MJD; an image's part from a `type` of `model` or
`residual`; and the names `MODEL`, `RESIDUAL`, `DIRTY` and `MASK_DECONVOLVE`. Every rule that
differed between the rows was written, tested and documented twice, and the only stores that needed
the second row were not schema-conformant: XRADIO 1.2.4 and earlier, hand-written fixtures, and the
site converter's output.

## Decision

**One version.** Only stores that name `schema_version: "0.0.2"` and are one are read. The row table
goes with the second row (`src/schema/xradio/schema_version.h` holds the version and XRADIO's image
names): a row held only what differed between versions, and with one version nothing does.

**What is refused, and as what.**

- A store that names no `schema_version` is an XRADIO image dataset of a version this library does
  not read. It is refused with `DiagnosticCode::unsupported_schema_version` first among the probe's
  diagnostics and `ErrorCode::unsupported_schema_version` from `Dataset::Open`, as any other
  version is, and the message says to write it again with XRADIO 1.2.5 or later.
- A store that names 0.0.2 but fails ADR 0022's checks -- it has data groups, the sky image they
  name, if they name one, is typed `"sky"`, and its frequency coordinate, if it has one, gives a frame
  and units -- is malformed: `invalid`, with `invalid_metadata` and the reason. The checks stay
  because the version alone does not say a store conforms: XRADIO 1.2.5 stamps 0.0.2 on any store it
  upgrades on read, including stores its own `check_image` rejects.
- A store whose data groups name no sky image is not malformed for that
  ([#40](https://github.com/CARTAvis/carta-zarr/issues/40)). The schema does not require one: XRADIO
  writes data groups of an aperture plane, and of a point spread function alone. Such a store is
  decided by ADR 0009's rule as any other is: one of complex aperture-plane variables answers
  `no_match`, saying they are not openable, and one with an image this library opens matches, a real
  aperture-plane image among them (ADR 0024). Its frequency
  coordinate is still checked: an older store with no sky is stamped 0.0.2 by XRADIO 1.2.5's upgrade
  as readily as one with a sky.
- Zarr that is not an image dataset is not refused at all, whether it names a version or not, as
  before.

`DiagnosticCode::unversioned_schema` is removed. A dataset that opens always reports
`schema_version` `"0.0.2"`; a probe still reports what a refused store names, or nothing.

**Rules that are not about versions stay** (ADR 0022's list): a flag as a Zarr `bool` or an `int8`
declared `dtype: bool`; ADR 0019's flag search, which a point spread function's flag still needs; a
store without consolidated metadata; beam labels found by name; an untyped boolean not taken for an
image; and a sky image of any real type. XRADIO's FITS reader still writes an integer FITS image as
an integer sky image, although `schema.py` allows only float, so refusing one would refuse what
XRADIO's own conversion wrote. Its pixels are read as float32, like any other.

## Why now, and in 0.2.0

ADR 0022 said retiring the second row would wait until the site converter wrote stores that pass
`check_image`, and the data converted with it had been converted again. That wait is given up: the
converter moves to XRADIO's own `write_image` once XRADIO 1.2.6 gives it chunks, shards and a
compressor, and everyone testing carta-zarr regenerates their data with it then. Until then the
converter's current output, the converted cubes on the test machines and the read-path test set do
not open.

0.2.0 is not tagged, and it already changes the ABI, so removing a `DiagnosticCode` goes into it
rather than into 0.3.

## Consequences

ADR 0003, which described how a store naming no version was read, describes nothing this library
reads any more.

The fixtures are 0.0.2: `minimal`, `uniform_beam`, `time_axis` and the `pixels` fixtures name the
version and are written as 0.0.2 writes them; XRADIO 1.2.3's conformance fixtures and the hand-written
`legacy` store are deleted. carta-backend-test-data's copies are regenerated with them.

A new schema version is refused until this library reads it. Reading it beside 0.0.2 brings back a
table of what differs between the two, and the rule that picks one.
