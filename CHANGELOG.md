# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Before 1.0 a
minor release may change the ABI, and the soname says so ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)).

## [Unreleased]

### Added
* A store naming XRADIO image schema 0.0.2, which XRADIO 1.2.5 writes, is read as that schema defines it: the spectral frame from `frequency`'s `frame`, and the unit from its `units`. A store naming no version, or naming 0.0.2 without being one, is read as XRADIO 1.2 wrote it, with a `DiagnosticCode::unversioned_schema` saying so; any other version is refused with `ErrorCode::unsupported_schema_version` and a `DiagnosticCode` of that name ([ADR 0022](docs/adr/0022-a-store-is-read-by-the-schema-version-it-names.md), [#20](https://github.com/CARTAvis/carta-zarr/issues/20)).
* `carta-zarr-bench pixels` writes a run of channels as the library reads them, with the pixel mask applied or not, as raw float32 on stdout ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).
* `tools/testset`: builds the standard read-path test set -- a pancake (7763 x 4742 x 256) and a cigar (512 x 512 x 30000), each a FITS cube and the Zarr layouts the site's xradio converter makes from it -- and checks that each Zarr is the layout asked for, with its flag, holding its FITS cube's pixels ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).
* `tools/zarr-bench/generate.py --synthetic --fits PATH` writes the synthetic cube as a FITS file, NaN where it is not data ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).

### Changed
* A probe and a dataset descriptor report the schema version the store names, or an empty string when it names none, rather than `"1.2"` for every store ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).
* A sky image typed `"sky"` that a data group other than `base` names as its sky reports that group's name as its role (`model`, `residual`, `dirty`), as XRADIO writes a dict of images and every schema 0.0.2 store, rather than `"sky"`. Images are listed with XRADIO's own names first, `SKY_MODEL`, `SKY_RESIDUAL` and `SKY_DIRTY` among them, so a dataset without `SKY` opens on one of those rather than on its point spread function ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).
* `ReadOptions::read_budget_bytes` bounds the memory a read holds beyond the caller's destination, not the bytes it decodes: a chunk is charged about three times what it decodes to, for its compressed bytes and the codec's buffer, plus the buffers the library allocates for it. A budget a caller states buys about a third of the chunks it did. The default aims at two chunks for every decode thread of the context, held between 256 MiB and 2 GiB ([ADR 0021](docs/adr/0021-a-read-budget-bounds-what-a-read-holds.md), [#1](https://github.com/CARTAvis/carta-zarr/issues/1)).
* `Image::Read` keeps to the budget whether or not it is given a progress callback, and reads a piece too large at one chunk deep in parts, so a plane of an image in large chunks no longer holds all of them at once. A read of many planes no longer holds a flag for its whole extent or lets the compressed chunks queue without bound when the decode threads fall behind. A read with a progress callback is cut into pieces of no more than sixteen chunks, unless one chunk along the cut is more, so that a cursor spectrum reports as often on many decode threads as on few ([#1](https://github.com/CARTAvis/carta-zarr/issues/1)).
* The parts a piece too large at one chunk deep is read in go straight into their place in the caller's destination rather than being gathered in a buffer of the library's and copied, so a read with no progress callback is cut to its budget too, but for a single plane whose budget affords a chunk for each decode thread, which is read whole as before. A read of 128 planes of the 7763 x 4742 test cube on four decode threads holds about 90 MiB rather than about 950 MiB, or 2.3 GiB with its flag, for a tenth or so more time; a read of many planes with a progress callback is faster ([#11](https://github.com/CARTAvis/carta-zarr/issues/11)).
* `tools/zarr-bench/generate.py --synthetic` makes a continuum-subtracted HI cube calibrated against ASKAP's, on ASKAP's frequency axis, instead of continuum point sources inside a circle. `--sources` and `--nan-radius` are replaced by `--line-sources`, `--point-sources`, `--extended-sources`, `--footprint-fill` and `--flagged-channels`, and `--flagged-range` flags a run of channels whole; `--frequency-start`, `--channel-width` and `--rest-frequency` set the frequency axis; synthetic datasets written before are not reused ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).

### Removed
* Forms no XRADIO writes are no longer read: an image's role from the v2 draft's `image_type`, an equinox given as a number, and an observation date given as a string. `ObservationInfo::date_obs`, which only a string date filled, is gone. This changes the ABI ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).

### Fixed
* A cube in the geocentric spectral frame reports `SpectralCoordinate::system` as `GEO` rather than `GCRS`, which XRADIO writes as the observer and casacore does not accept, so carta-backend could not open it ([#18](https://github.com/CARTAvis/carta-zarr/issues/18)).
* A masked read under a budget smaller than one chunk reads a chunk at a time instead of failing with `buffer_too_small` ([#1](https://github.com/CARTAvis/carta-zarr/issues/1)).
* A boolean variable on the sky plane that is not typed `flag` is listed as not openable because it is boolean, rather than because it is complex. The `int8` that xarray writes for a bool is treated the same, rather than listed as an openable image ([#9](https://github.com/CARTAvis/carta-zarr/issues/9)).
* Plane histograms of a cube chunked deep and narrow no longer take twice as long: a read of many channels is split among the workers by channel, rather than each channel by rows, when its busiest worker bins no more rows that way ([#6](https://github.com/CARTAvis/carta-zarr/issues/6)).

## [0.1.0]

The first release.

### Added
* Reading XRADIO 1.2 image datasets stored as Zarr v3 on a local filesystem, with or without consolidated metadata.
* Probing a path for an image dataset, and listing the images in one with what each will open with or why it will not.
* Descriptors of a dataset and its images: axes by role, coordinates and their linear descriptions, the observation, and the chunk geometry.
* Pixel reads of any range along each axis as `float32`, with the image's flag applied as its pixel mask.
* Spectral reductions of many regions in one pass: pixel and NaN counts, sum, sum of squares, spread, minimum and maximum.
* Plane histograms, and cube histograms whose range is not known in advance.
* Beam tables, dataset sizes, chunk prefetch, and read-ahead for animations.
* Cancellation, deadlines, progress, and shareable decoded-chunk cache pools for every read and reduction.
* CMake package `CartaZarr` with target `CARTA::zarr`; TensorStore and its dependencies are linked in statically and none of their symbols are exported.
* Third-party licence notices installed beside the library, and an offline build from packed dependency archives.
* Storage benchmark and layout tools for Lustre and BeeGFS under `tools/zarr-bench`.
