# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Before 1.0 a
minor release may change the ABI, and the soname says so ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)).

## [Unreleased]

### Added
* `carta-zarr-bench probe --describe` adds what the image means as the library reads it: stored type and unit, pixel mask, direction, spectral, polarization and time coordinates, the observation, and the beam of every plane.
* `tools/testset`: builds the standard read-path test set -- a pancake (7763 x 4742 x 256) and a cigar (512 x 512 x 30000), each a FITS cube and the Zarr layouts the site's xradio converter makes from it -- and checks that each Zarr holds its FITS cube's pixels.

### Changed
* `tools/zarr-bench/generate.py --synthetic` makes a continuum-subtracted HI cube calibrated against ASKAP's, on ASKAP's frequency axis, instead of continuum point sources inside a circle. `--sources` and `--nan-radius` are replaced by `--line-sources`, `--point-sources`, `--extended-sources`, `--footprint-fill` and `--flagged-channels`; synthetic datasets written before are not reused.

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
