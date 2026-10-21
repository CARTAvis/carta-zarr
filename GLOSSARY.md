# carta-zarr

A standalone C++ reader library for XRADIO Zarr images, written as the Zarr reader of
`carta-backend` and kept apart from it so that CARTA and other C++ consumers can read these datasets
without depending on TensorStore. It depends on neither casacore nor CARTA's protobuf.

## Language

### The dataset

**Image dataset**:
One Zarr v3 group holding a complete XRADIO image, marked by the root attribute `type:
"image_dataset"`.
It is the unit a consumer opens by path. Where the subject is the schema version it declares,
"store" is the usual word for it, as ADR 0022 and the changelog use it.
_Avoid_: file, image (on its own), group; "store" for it outside talk of its schema version

**Data variable**:
One named array inside an image dataset. Some are images, others are flags, beam tables, or
normalization factors.
_Avoid_: array, dataset variable, HDU

**Coordinate**:
A named axis of the image dataset (`time`, `frequency`, `polarization`, `l`, `m`, `u`, `v`), stored
once at dataset level and referenced by every data variable that uses it.
_Avoid_: dimension, axis array

**Data group**:
A named set of related data variables produced by one imaging run, mapping a role (`sky`, `flag`,
`primary_beam`, …) to the variable that fills it. Several data groups can reference the same
variable.
_Avoid_: image group, group, dataset group

### What can be opened

**Image**:
A data variable a consumer can open and read pixels from: typed as one of schema 0.0.2's images
(`sky`, `point_spread_function`, `primary_beam` or `mask` on the sky plane; `aperture`, `visibility` or
`uv_sampling` on the aperture plane), carrying every coordinate of that plane, and real-valued.
Identified by its variable name, which is unique within the dataset. An untyped variable is not an
image, whatever its dimensions.
_Avoid_: image array, HDU, plane, layer

**Sky plane**:
The tangent plane to the celestial sphere, carrying the `l` and `m` direction-cosine coordinates.
_Avoid_: image plane, lm plane

**Aperture plane**:
The Fourier-conjugate plane, carrying the `u` and `v` coordinates in wavelengths. In a real-valued
image on it, `u` and `v` play the spatial roles, and its descriptor has an `aperture` coordinate in
place of a direction (ADR 0024). A complex image is not opened, on either plane.
_Avoid_: uv plane, visibility plane, Fourier plane

**Flag**:
A boolean data variable marking which pixels of its image are flagged, recognized by the attribute
`type: "flag"` rather than by its name. True means the pixel is flagged and dropped. Stored as a Zarr
`bool`. Tied to its image by the image's own `flag` attribute or by a root `data_groups` entry whose
`sky` is the image. Otherwise an image that is not a sky image, whose flag XRADIO cannot declare,
takes the flag named `FLAG_<image>`. It becomes the image's pixel mask; it is never itself an image.
_Avoid_: mask, internal mask, blanking

**Deconvolution mask**:
An image marking the region a deconvolution was allowed to work in. It is an ordinary image, opened
and displayed like any other.
_Avoid_: mask, flag

**Beam fit parameters**:
A data variable holding fitted major axis, minor axis, and position angle per plane, named by the
owning image's `beam_fit_params` attribute. Different images may point at different ones.
_Avoid_: beam table, beam array, restoring beam

### Reading

**Schema profile**:
A named, versioned description of how an image dataset is laid out, which the library matches a
store against. The only built-in profile is the XRADIO image profile.
_Avoid_: schema, format, flavor

**Schema version**:
The version of XRADIO's image schema a store names in its root's `schema_version`, `0.0.2` from
XRADIO 1.2.5. A store naming none was written by XRADIO 1.2.4 or earlier. Only 0.0.2 is read
(ADR 0023).
_Avoid_: XRADIO version, format version

**Probe**:
A cheap, read-only inspection that decides whether a path is a supported image dataset, without
opening any image in it. Listing the images is what opening the dataset does.
_Avoid_: detect, sniff, validate

**Qualification**:
Whether this library will open a data variable, and what it has to say about one it will not. It is
one decision rather than one per stage: the probe, the listing a consumer chooses an image from, and
describing an image all ask it, so none of them can answer differently about the same variable.
_Avoid_: validation, classification, screening, eligibility

**Transport**:
Where an image dataset's bytes live, and the only thing the library swaps out to read from
somewhere else. It returns one node's metadata verbatim and lists the nodes; it never decides what a
node means, so every transport is interpreted identically.
_Avoid_: driver, backend, kvstore, store

**Node inventory**:
What a store's hierarchy holds: every node under its one canonical name, whether it is a group, an
array, or something the library does not recognize, and for an array, whether its metadata parses.
Taken once for the life of a store, which is what makes a dataset's listing a snapshot. It decides
what each node is; what to do about one that is not a usable array stays with whoever asked.
_Avoid_: listing, node list, catalog, manifest

**Descriptor**:
The immutable metadata the library reports for a dataset or an image. It never contains pixels.
_Avoid_: info, header, metadata (on its own)

**Diagnostic**:
Something the library has to say about a store it nonetheless accepted: a coded note carried on a
probe result, an image entry, or a descriptor. Its counterpart is an error, and the two divide the
same ground: an error is a refusal, so nothing is returned, while a diagnostic accompanies an answer
the caller still gets. A variable that could not be opened, an axis with no linear description, and
a node skipped because its metadata would not parse are all diagnosed rather than refused, because
the rest of the dataset is still readable.

The rule that follows: a path that does not fill in a value leaves a diagnostic saying so. Without
one, the caller receives a default that looks like an answer, such as a reference pixel of zero or a
chunk shape equal to the image, and cannot tell it from a real one.
_Avoid_: warning, message, note, error (for this)

**Observation**:
What an image's own attributes say about the observing run behind it: the object, the observer, the
telescope and where it stood, and the date the observation carries. Optional throughout — a field
written in a type the library cannot read is skipped, and never closes the image.
_Avoid_: provenance, history, header

**Linear description**:
The reference pixel, reference value, and increment that let a consumer treat a coordinate as an
evenly spaced axis. A coordinate carries one only when its samples support it; without one the
consumer builds a tabular axis from the values instead. A direction axis is linear by construction
and always reports one.
_Avoid_: linear triple, FITS keywords, WCS

**Stored order**:
The order the coordinates appear in a data variable's own dimensions, as written on disk.
_Avoid_: native order, disk order, physical order

**Logical order**:
The coordinate order the library reports and reads in, chosen by the schema profile rather than by
the file.
_Avoid_: canonical order, CARTA order, display order

**Plane selection**:
Which planes of an image a reduction is over: a range along the spectral coordinate, one
polarization, and one time. The three come together on every reduction request, because a reduction
is always over whole planes. An ordinary read gives one range per axis instead, which is why it
takes no plane selection.
_Avoid_: slice, cube selection, channel range, plane range

**Pass**:
One ordered visit to every chunk an image read covers, made once and shared by every reduction that
wants those pixels. It decides how much to decode at a time, when to hand a result over, and how to
split the work across workers. It never knows what a region is and never does the arithmetic the
caller came for.
_Avoid_: walk, traversal, scan, loop

**Occupancy**:
The chunks a set of regions actually occupies, and which of those regions touch each one. Not the
chunks their bounding boxes cover: a thin cut along the diagonal has a box the size of the image and
touches one chunk per row. Worked out once from each region's own mask, and read throughout the pass
that follows.
_Avoid_: bucket, index, incidence, coverage, bounding box

**Footprint**:
One rectangle of chunks a pass reads along the spectrum, a slab at a time: small enough that one
layer of it fits in a read, and counted in the chunks it occupies, which is what progress is measured
in. A whole-plane pass cuts its footprints as bands of the plane; a reduction's are cut by its
occupancy from the chunk runs its regions touch, so a diagonal cut is read as one chunk per row
rather than as its bounding box.
_Avoid_: tile, window, block, rectangle

**Piece**:
One chunk-aligned part of an ordinary read, cut along the slowest-varying axis of which the request
selects more than one element. A read is always made of pieces, as many as its read budget needs;
one that fits is a single piece covering everything. Cutting only there keeps the finished part of
the destination a prefix rather than scattered, so a caller can be told how far along it is, or stop
it. A piece that holds more than the budget allows even at one chunk deep is read in **segments**,
cut along the faster axes, each straight into its place in the destination. A segment does not fill
a contiguous run of the destination, so it is not reported on its own; the piece is, once all of its
segments are in.
_Avoid_: block, slab, chunk, batch

**Read budget**:
How much memory one read may hold at once beyond the caller's own destination: the chunks it is
decoding, each counted at about three times its decoded size, and the buffers the library allocates
for them. Every read, pass and prefetch keeps to it, with or without a progress callback, except a
single plane read with no progress callback, which is read whole once its budget allows a chunk for
each decode thread. It is per read, and a read holds at least one chunk however small the budget.
_Avoid_: ceiling, memory limit, chunk budget

**Cache pool**:
Where a read keeps the chunks it decodes, for a later read to find instead of decoding them again.
A context shares one among every read made through it, which holds what is being looked at. A read
may bring its own instead: one of zero bytes, for a scan that will not come back to a chunk, or one
large enough for a walk that will come back, held only as long as the walk. A chunk is decoded
whole, so a read that needs any of it keeps all of it.
_Avoid_: tile cache, chunk cache policy, bypass

**Provisional histogram**:
What one task of a cube histogram's walk bins its pixels into, before the range is known: fine,
starting around the first pixel it sees and doubling to fit whatever arrives after. Each task keeps
its own, and none is ever merged with another. At the end, or whenever a caller asks for a
snapshot, each is re-aggregated onto the caller's bins over the extremes, which are exact, and the
counts are added. `provisional_bin_count` sets how fine it is.
_Avoid_: partial, accumulator, tile histogram, intermediate histogram

**Spread**:
How far a set of pixels lies from its own mean, given as the sum of their squared deviations from
it, `sum_sq_dev`; a sample variance is that divided by one less than the count. It is kept beside the
sum and the sum of squares rather than computed from them, because from those it is a difference of
two numbers that agree in every digit a double holds when the pixels are far from zero compared with
their spread. Two sets' spreads are combined from their counts and means. See ADR 0018.
_Avoid_: variance (for the sum itself), M2, second moment (on its own)

**Run**:
The chunks one plane of an image decodes (the plane rounded out to whole chunks, and one chunk deep
along every other axis), and so every plane that decodes the same ones. An animation playing through
a run is served from the cache after its first frame; the frame that enters the next one decodes all
of it.
_Avoid_: chunk row, slab, group, block

**Read-ahead**:
Decoding the next run of a playing animation while the frames of this one play from the cache, so
that the frame entering it does not stall. It is told after every frame what was played, whether it
was late, and what comes next. One decode runs at a time, and reading ahead stops once a frame is
late while one runs. A single decode ahead is a prefetch.
_Avoid_: lookahead, preload, prefetching (for the whole policy)
