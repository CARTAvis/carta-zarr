# A store is read by the schema version it names

*Its two rows are superseded by [ADR 0023](0023-only-xradio-image-schema-0-0-2-is-read.md): only
schema 0.0.2 is read, and a store naming no version is refused. Its checks, its refusal of an unknown
version and what the schema is held to stand. Its checks apply to a store whose data groups name a
sky image ([#40](https://github.com/CARTAvis/carta-zarr/issues/40)): the schema requires none.*

XRADIO 1.2.5 gives its image format a version. `IMAGE_SCHEMA_VERSION = "0.0.2"` in
`xradio/image/schema.py` is written to every store's root as `schema_version`, `schema.py` defines
the format, and `check_image` validates a dataset against it. XRADIO 1.2.4 and earlier wrote no
version, and their stores are what ADR 0003 was written against: a draft schema document and code
that disagreed, with the code taken as the authority.

This library read every store by one set of rules, written against XRADIO 1.2.3, with fallbacks for
older and hand-written forms among them, and reported the constant `"1.2"` as every store's schema
version. Nothing said which rule belonged to which XRADIO, so none could be retired, and a field
0.0.2 added could not be read without guessing whether an older store had it.

## What a store is read by

**A row per version.** The XRADIO profile reads a store by one of a table of rows
(`src/schema/xradio/schema_row.h`), one per schema version it reads: 0.0.2, and the stores that
name none. A row holds only what differs between versions:

- where the spectral frame is: `frequency.attrs.frame` for 0.0.2, and the reference frequency's
  observer otherwise, upper-cased with astropy's GCRS read as casacore's GEO;
- whether a frequency with no units takes its reference frequency's;
- whether an observation date with no format is an MJD;
- the fields only 0.0.2 writes: a channel's width, the increment of an axis of one channel, and a sky
  image's `sub_type`;
- XRADIO's names for a dataset's images, in the order they are listed: XRADIO 1.2 named a
  deconvolution's products both `MODEL` and `SKY_MODEL`, by whether it was handed a list of them or a
  dict, and 0.0.2 only the second.

A rule that is the same for every version is code, written once where it is applied. A sky image's
part (residual, model, dirty) is one: it is the name of the data group whose sky it is when the image
is typed `"sky"`, which 0.0.2 always does and XRADIO 1.2 did for a dict of images, and the image's own
`type` otherwise, which is where XRADIO 1.2 put it for a list.

**Which row.** A store naming 0.0.2 is read by 0.0.2's row if it is one: its data groups name a sky
image, that image is typed `"sky"`, and its frequency coordinate, if it has one, gives a frame and
units. Those three, rather than the whole schema, because they are what the rows read differently;
and checked at all because the version alone does not say a store conforms. XRADIO 1.2.5 upgrades an
older store when it reads one, and stamps 0.0.2 on what it writes back whether or not the upgrade
could reach those fields. A store naming 0.0.2 that is not one, and a store naming no version, is read
by the unversioned row, and a `DiagnosticCode::unversioned_schema` says why.

Any other version is refused, with `DiagnosticCode::unsupported_schema_version` first among the
probe's diagnostics and `ErrorCode::unsupported_schema_version` from `Dataset::Open`: a file a newer
build may open, which is neither a broken file nor one this library is not for. Only an image dataset
is refused for its version. XRADIO's measurement sets name a `schema_version` too, and Zarr of
something else is not called an image of an unknown schema for carrying one.

A probe reports the version the store names, or none: `"0.0.2"` for a store read as 0.0.2, for one
read as XRADIO 1.2's despite naming 0.0.2, and an empty string for a store that names nothing.

## What the schema is held to

For a store of a version XRADIO publishes, `schema.py` is the authority, and this supersedes ADR 0003
for those stores. A store `check_image` accepts must open. Where XRADIO's own writer disagrees with
`schema.py`, that is reported to XRADIO, and accommodated here only where refusing would refuse what
XRADIO's own conversion writes, as an integer sky image below. ADR 0003 still describes the
unversioned row, which reads what the code of XRADIO 1.2 wrote, until that row is retired.

Some rules are kept for every version, because they defend against stores found in practice rather
than against any version's layout:

- a boolean on the sky plane that is not typed `flag` is not offered as an image;
- a flag may be a Zarr `bool` or an `int8` declared `dtype: bool`;
- an image's flag is found by ADR 0019's three steps, including the guess by dimensions, which
  XRADIO 1.2.5's `FLAG_POINT_SPREAD_FUNCTION` needs because nothing links it to its image;
- a store need not have consolidated metadata;
- beam labels are found by name, in any order;
- a time's scale and format are read in any case;
- a sky image may be of any real type. `schema.py` allows float32 and float64, but XRADIO 1.2.5's FITS
  reader keeps an integer FITS image's type (`BITPIX` 16 is int16), and its writer stores it as it
  is. That disagreement is XRADIO's to settle, and is reported to it; until it is, refusing such a
  store would refuse what XRADIO's own conversion wrote.

Forms no XRADIO writes are no longer read: the v2 draft's `image_type` as an image's role, an equinox
given as a number, and an observation date given as a string (and `ObservationInfo::date_obs` with
it).

## Why rows of a profile and not of the profile table

ADR 0007 kept `SchemaProfile::BuiltIn()` a table of one entry for a second profile, and said that a
version of the XRADIO profile handled inside it rather than as a row of the table would reopen that.
The versions are inside it. Two entries of the table would each have to decide, of every XRADIO
store, whether it was theirs, and that decision is one rule over the store's `schema_version` and the
three checks above; split between two entries it would be stated twice and could disagree. A version
is not a second schema, and the table's own reason, that a schema id is public and a consumer may name
one this library does not have, is unchanged. The table is kept, and ADR 0007's last paragraph is
amended to say so.

## Consequences

Retiring XRADIO 1.2's stores is deleting a row, and the rules only it reads with it. That waits until
the site converter writes stores that pass `check_image`, the converted cubes and `tools/testset`
have been converted again, and the fixtures here and in carta-backend are 0.0.2. Until then README
keeps its promise of XRADIO 1.2.x.

A new schema version is a new row, and a store of it is refused until there is one. The conformance
fixtures are what tells the two apart: they are written by the pinned XRADIO, and XRADIO 1.2.3's are
kept beside them, as `conformance_1_2` and `conformance_flagged_1_2`, for as long as the unversioned
row is.
