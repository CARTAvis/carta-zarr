# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Before 1.0 a
minor release may change the ABI, and the soname says so ([ADR 0020](docs/adr/0020-before-1-0-the-minor-version-is-the-abi.md)).

## [Unreleased]

### Changed
* `ReadOptions::read_budget_bytes` bounds the memory a read holds beyond the caller's destination, not the bytes it decodes: a chunk is charged about three times what it decodes to, for its compressed bytes and the codec's buffer, plus the buffers the library allocates for it. A budget a caller states buys about a third of the chunks it did. The default aims at two chunks for every decode thread of the context, held between 256 MiB and 2 GiB ([ADR 0021](docs/adr/0021-a-read-budget-bounds-what-a-read-holds.md)).
* `Image::Read` keeps to the budget whether or not it is given a progress callback, and reads a piece too large at one chunk deep in parts, so a plane of an image in large chunks no longer holds all of them at once. A read with no progress callback whose budget affords a chunk a decode thread is cut into layers of its chunks -- a chunk deep along the spectral axis and the whole plane -- rather than to its budget, since cut finer it is gathered in segments and runs slower; a plane is read whole as before, and a read of many planes no longer holds a flag for its whole extent or lets the compressed chunks queue without bound when the decode threads fall behind. A read with a progress callback is cut into pieces of no more than sixteen chunks, unless one chunk along the cut is more, so that a cursor spectrum reports as often on many decode threads as on few.

### Fixed
* A masked read under a budget smaller than one chunk reads a chunk at a time instead of failing with `buffer_too_small`.

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
