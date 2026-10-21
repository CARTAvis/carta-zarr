# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Before 1.0 a
minor release may change the ABI, and the soname says so ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)).

## [Unreleased]

### Added
* Real-valued images of the aperture plane open, as CASA's `ia.fft` writes an amplitude or a phase and XRADIO converts it: `u` and `v` play `spatial_x` and `spatial_y`, so reads and reductions reach them as they reach `l` and `m`, and `ImageDescriptor::aperture`, a new `ApertureCoordinate`, gives their unit, reference pixel, reference value (`crval`) and increment in place of a direction. A dataset of them alone needs no `coordinate_system_info`. Their flag is found by its shape, and they have no beams. A complex variable on either plane is listed and not opened, with `unsupported_data_type`, and `DiagnosticCode::unsupported_coordinate_plane` is no longer emitted. An aperture-plane variable without its own `u` or `v` coordinate is listed as malformed and does not open, and the sky images beside it do. This changes the ABI ([ADR 0024](docs/adr/0024-the-spatial-roles-are-the-images-plane.md), [#41](https://github.com/CARTAvis/carta-zarr/issues/41)).
* `ObservationInfo::pointing_center`: where the telescope pointed, FITS `OBSRA`/`OBSDEC`, as a `SkyDirection` in degrees with its frame named as `DirectionCoordinate::reference_frame` names one. Reported whenever the image has one, including when it equals the reference direction, which is what XRADIO writes for a FITS file without `OBSRA`. `ImageDescriptor::user_keywords`: the dataset's and the image's `user` keywords, the image's winning, each a boolean, an integer, a number or a string; any other value is left out with a diagnostic. This changes the ABI ([#34](https://github.com/CARTAvis/carta-zarr/issues/34)).
* `CARTA_ZARR_THREAD_SANITIZER` builds the library and everything it links, TensorStore included, with ThreadSanitizer, for running the threading tests under it. It patches one memory order in TensorStore's task queue that ThreadSanitizer cannot follow; no other build is patched ([#26](https://github.com/CARTAvis/carta-zarr/issues/26)).

### Changed
* A read, a prefetch or a reduction whose request does not fit the image's axes reports `ErrorCode::invalid_slice`, and a read whose destination is smaller than its selection reports `ErrorCode::buffer_too_small`, as `docs/troubleshooting.md` already said; both were `invalid_argument`. That covers a range past an axis's end, a zero stride or count, the wrong number of ranges, a spectral, polarization or time index outside the image, and a reduction's region that is empty or falls off the image. `invalid_argument` stays for other caller mistakes, such as a missing sink, no regions or too many, or a region's mask whose length is not its box's ([#36](https://github.com/CARTAvis/carta-zarr/issues/36)).
* `Dataset::Open` reports a dataset a schema profile recognised and found malformed, the probe's `SchemaMatchKind::invalid`, with the new `ErrorCode::invalid_schema` rather than `invalid_metadata`. `invalid_metadata` stays for Zarr's own problems, such as a root `zarr.json` that will not parse (one that cannot be read is `io_error`, as before), so a consumer can tell a broken image dataset from unreadable Zarr without probing the store again. A dataset a profile matched with no image variables is `invalid_schema` too. The probe's diagnostics are unchanged. Adding the enumerator renumbers the codes after `ambiguous_schema`, which changes the ABI ([#37](https://github.com/CARTAvis/carta-zarr/issues/37)).
* `Image::Read` promises what it already did: when a read fails or is cancelled, the elements its progress callback last reported as written are final, so a caller can keep them and discard only the rest. The documentation said to discard the whole destination. A read without a progress callback has no such prefix ([#35](https://github.com/CARTAvis/carta-zarr/issues/35)).
* Public names made consistent. This changes the API and the ABI ([#39](https://github.com/CARTAvis/carta-zarr/issues/39)):
  * `SpectralCoordinate::channel_frequencies` is `values`, as `TemporalCoordinate` names its own, and for what it is: one value a channel, in `unit`, which nothing holds to Hz.
  * `SpectralCoordinate::system` is `reference_frame`, as `DirectionCoordinate` names its frame.
  * `ObservationInfo::timesys` is `time_scale` and `mjd_obs` is `observation_mjd`, words rather than FITS keywords, as the struct's other fields are.
  * `HistogramRequest::bins` and `CubeHistogramRequest::bins` are `bin_count`, and `provisional_bins` is `provisional_bin_count`, as `HistogramBlock` reports the number.
  * `Context::NewCachePool` is `CreateCachePool`, as every other factory is `Create` or `Open`.
  * `CubeHistogramProgress::progress` is `completeness`, as the blocks report a fraction done.
  * `Dataset::Size`'s parameter is `measure_timeout`, after `SizeBasis::measured`; a parameter name, with no effect on the ABI.

### Removed
* Public fields that repeated what another field already says, or that nothing set any more. This changes the API and the ABI ([#38](https://github.com/CARTAvis/carta-zarr/issues/38)):
  * `ImageEntry::image_role` and `ImageDescriptor::image_role`. The variable's `id` and the image's `data_groups` say the same, and no consumer read it. `carta-zarr-bench`'s JSON no longer has a `role` key.
  * `ImageDescriptor::has_pixel_mask`: an image has a pixel mask exactly when `pixel_mask_id` is not empty. `carta-zarr-bench`'s JSON still reports `has_pixel_mask`, derived the same way.
  * `AxisDescriptor::unit`. It was the same `units` attribute the coordinates report: `SpectralCoordinate::unit`, `TemporalCoordinate::unit` and `ApertureCoordinate::unit`. It was empty for `l`, `m` and `polarization`.
  * `ChunkGeometry::grid_shape`, which is each axis's length over `chunk_shape`, rounded up. It counted inner chunks, not shards or files.
  * `DiagnosticCode::unsupported_coordinate_plane`, no longer emitted since aperture-plane images open. The codes after it are renumbered.

### Fixed
* A store naming schema 0.0.2 whose data groups name no sky image is no longer refused as malformed. The schema does not require one, and XRADIO writes such stores: an aperture plane alone (`APERTURE`, `FLAG_APERTURE`), or a point spread function alone. A store of complex aperture-plane variables answers `no_match`, saying they are not openable, and one with an image this library opens matches. A sky a data group names must still be typed `"sky"`, and a frequency coordinate must still have `frame` and `units` whatever the groups name ([#40](https://github.com/CARTAvis/carta-zarr/issues/40)).

## [0.2.0]

### Added
* `ImageDescriptor::sub_type`: what a sky image of a schema 0.0.2 store measures (`Intensity`, `SpectralIndex`, ...), as XRADIO names it from casacore's image type or FITS `BTYPE`. Empty for an image naming none. This changes the ABI ([#21](https://github.com/CARTAvis/carta-zarr/issues/21)).
* A spectral axis of one channel in a schema 0.0.2 store, a continuum image's, is given its `channel_width` as its increment, measured from that channel, where before it had no linear description ([#21](https://github.com/CARTAvis/carta-zarr/issues/21)).
* A store naming XRADIO image schema 0.0.2, which XRADIO 1.2.5 and later write, is read as that schema defines it: the spectral frame from `frequency`'s `frame`, and the unit from its `units`. It is the only schema version read: any other is refused with `ErrorCode::unsupported_schema_version` and a `DiagnosticCode` of that name ([ADR 0022](docs/adr/0022-a-store-is-read-by-the-schema-version-it-names.md), [ADR 0023](docs/adr/0023-only-xradio-image-schema-0-0-2-is-read.md), [#20](https://github.com/CARTAvis/carta-zarr/issues/20), [#29](https://github.com/CARTAvis/carta-zarr/issues/29)).
* `carta-zarr-bench pixels` writes a run of channels as the library reads them, with the pixel mask applied or not, as raw float32 on stdout ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).
* `tools/testset`: builds the standard read-path test set -- a pancake (7763 x 4742 x 256) and a cigar (512 x 512 x 30000), each a FITS cube and the Zarr layouts the site's xradio converter makes from it -- and checks that each Zarr is the layout asked for, with its flag, holding its FITS cube's pixels ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).
* `tools/zarr-bench/generate.py --synthetic --fits PATH` writes the synthetic cube as a FITS file, NaN where it is not data ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).

### Changed
* The version is 0.2.0, and the soname and symbol versions 0.2: the changes to `DiagnosticCode`, `ImageDescriptor` and `ObservationInfo` change the ABI ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)). A consumer asks for `find_package(CartaZarr 0.2 ...)`.
* A probe reports the schema version the store names, or an empty string when it names none, rather than `"1.2"` for every store; a dataset that opens names `"0.0.2"` ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).
* A sky image typed `"sky"` that a data group other than `base` names as its sky reports that group's name as its role (`model`, `residual`, `dirty`), as schema 0.0.2 writes it, rather than `"sky"`. Images are listed with XRADIO's own names first, `SKY_MODEL`, `SKY_RESIDUAL` and `SKY_DIRTY` among them, so a dataset without `SKY` opens on one of those rather than on its point spread function ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).
* `ReadOptions::read_budget_bytes` bounds the memory a read holds beyond the caller's destination, not the bytes it decodes: a chunk is charged about three times what it decodes to, for its compressed bytes and the codec's buffer, plus the buffers the library allocates for it. A budget a caller states buys about a third of the chunks it did. The default aims at two chunks for every decode thread of the context, held between 256 MiB and 2 GiB ([ADR 0021](docs/adr/0021-a-read-budget-bounds-what-a-read-holds.md), [#1](https://github.com/CARTAvis/carta-zarr/issues/1)).
* `Image::Read` keeps to the budget whether or not it is given a progress callback, and reads a piece too large at one chunk deep in parts, so a plane of an image in large chunks no longer holds all of them at once. A read of many planes no longer holds a flag for its whole extent or lets the compressed chunks queue without bound when the decode threads fall behind. A read with a progress callback is cut into pieces of no more than sixteen chunks, unless one chunk along the cut is more, so that a cursor spectrum reports as often on many decode threads as on few ([#1](https://github.com/CARTAvis/carta-zarr/issues/1)).
* The parts a piece too large at one chunk deep is read in go straight into their place in the caller's destination rather than being gathered in a buffer of the library's and copied, so a read with no progress callback is cut to its budget too, but for a single plane whose budget affords a chunk for each decode thread, which is read whole as before. A read of 128 planes of the 7763 x 4742 test cube on four decode threads holds about 90 MiB rather than about 950 MiB, or 2.3 GiB with its flag, for a tenth or so more time; a read of many planes with a progress callback is faster ([#11](https://github.com/CARTAvis/carta-zarr/issues/11)).
* `tools/zarr-bench/generate.py --synthetic` makes a continuum-subtracted HI cube calibrated against ASKAP's, on ASKAP's frequency axis, instead of continuum point sources inside a circle. `--sources` and `--nan-radius` are replaced by `--line-sources`, `--point-sources`, `--extended-sources`, `--footprint-fill` and `--flagged-channels`, and `--flagged-range` flags a run of channels whole; `--frequency-start`, `--channel-width` and `--rest-frequency` set the frequency axis; synthetic datasets written before are not reused ([#2](https://github.com/CARTAvis/carta-zarr/issues/2)).

### Removed
* Image stores naming no schema version, which XRADIO 1.2.4 and earlier wrote and 0.1.0 read, are refused with `ErrorCode::unsupported_schema_version`, and the message says to write them again with XRADIO 1.2.5 or later. A store naming 0.0.2 that has no data group naming a sky image, a sky image not typed `"sky"`, or a frequency coordinate without `frame` and `units` is refused with `ErrorCode::invalid_metadata`: XRADIO 1.2.5 stamps 0.0.2 on an older store it upgrades on read whether or not the upgrade reached them. With those stores go the rules only they needed: the spectral frame from the reference frequency's observer, a unit from the reference frequency, an observation date with no format read as an MJD, an image's part read from a `type` of `model` or `residual`, and the names `MODEL`, `RESIDUAL`, `DIRTY` and `MASK_DECONVOLVE` listed among XRADIO's own. Integer sky images, which XRADIO's FITS reader still writes, are read as before ([ADR 0023](docs/adr/0023-only-xradio-image-schema-0-0-2-is-read.md), [#29](https://github.com/CARTAvis/carta-zarr/issues/29)).
* Forms no XRADIO writes are no longer read: an image's role from the v2 draft's `image_type`, an equinox given as a number, and an observation date given as a string. `ObservationInfo::date_obs`, which only a string date filled, is gone. This changes the ABI ([#20](https://github.com/CARTAvis/carta-zarr/issues/20)).

### Fixed
* A cube in the geocentric spectral frame reports `SpectralCoordinate::system` as `GEO`, casacore's name, rather than `GCRS`, which casacore does not accept, so carta-backend could not open it ([#18](https://github.com/CARTAvis/carta-zarr/issues/18)).
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
