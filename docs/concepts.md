# Core concepts: datasets, images and reads

## Dataset versus image

An **image dataset** is one Zarr v3 group opened by path. It can contain several
image data variables, shared coordinates, flags and beam fit parameters.
An **image** is one supported data variable opened by its variable name.

```{mermaid}
flowchart LR
    D["Dataset: cube.zarr"] --> S["ImageEntry: SKY"]
    D --> R["ImageEntry: SKY_RESIDUAL"]
    D --> C["Shared coordinates"]
    S -->|"OpenImage(id)"| I["Image handle"]
    I -->|"Read(request, buffer)"| P["float32 pixels"]
```

The names above are illustrative: use the IDs the actual listing reports.
An image role such as sky or residual describes a variable's purpose; it is
not the ID used to open it.
For variables typed `sky`, a non-base data group's name can determine the role,
for example `model`, `residual` or `dirty`. An opened image's `sub_type` instead describes
what it measures, such as intensity or spectral index. See
[schema versions](reference/schema-versions.md) for version-dependent metadata.

| Object | Responsibility |
|---|---|
| `Context` | Shared decoded-chunk cache, I/O and decode resources, reduction workers |
| `Dataset` | Schema, image listing, default image ID and dataset size |
| `ImageEntry` | Information for choosing an image without opening it |
| `Image` | Detailed descriptor, pixel reads, reductions, histograms and beams |
| `CachePool` | An optional separate decoded-chunk cache for selected operations |

`ImageEntry` contains `id`, `openable`, `diagnostics`, `image_role` and `axes`.
The axes match those `OpenImage` would report for an openable image; they are
empty for an unopenable entry. Descriptors contain metadata, not pixel buffers.

## Axes and output order

A `ReadRequest` has one `Range {start, count, stride}` per descriptor axis.
Find an axis with `AxisIndex(axes, AxisRole::spectral)`, for example. The current
XRADIO profile reports l, m, frequency, polarization, time in that logical order,
regardless of the array's stored order. Avoid hard-coding that order.

The output is densely packed with logical axis 0 varying fastest. In the current
XRADIO profile, `l` (x) is axis 0, so it varies fastest in the output. `storage_index`
reports where an axis lives in the stored array, which can differ from its logical position.

## Flag versus region mask

An image's **flag** marks bad pixels: true means excluded. How reads apply it
is one of the [rules below](#rules-that-apply-to-every-operation); set
`apply_pixel_mask = false` to read the stored values without it.

A reduction's **region mask** (`RegionMask`) selects pixels within a bounding
box: nonzero means included. A deconvolution mask, `MASK`, is neither: it is an image of
its own, not a flag. A float-valued one opens like any other image; one
stored as booleans is listed but not openable. These three serve different
purposes.

## Rules that apply to every operation

These hold for all reads, reductions and histograms; the guides link here
rather than repeating them.

- **Check every `Result`** before using its value. Failures carry an error
  code, message and node path.
- **A failed or cancelled read leaves the destination undefined.** Discard it,
  even if progress had already reported a completed prefix.
- **Flags are applied by default.** A flagged pixel is output as NaN and
  excluded from reductions; stored NaNs stay NaN either way.
- **Borrowed views last only for the callback.** A reduction block or histogram
  `Counts` view is valid inside the sink; copy what you need to keep.
- **Handles are shared; destinations belong to you.** Keep your buffers and
  region masks alive until the operation returns.

## Results, diagnostics and lifetimes

A `Result<T>` holds either a value or an error. A diagnostic is different: it
accompanies a successful listing or descriptor and explains an unsupported
candidate or missing optional information.

`Context`, `Dataset`, `Image` and `CachePool` are shared handles created by
factories. Copies share their resources, and moving one preserves the source
handle too. Use `std::optional` when your own object needs a not-yet-created
handle. A reference obtained from `descriptor()` must not outlive the object
that owns the underlying descriptor.

## Threads

- **Handles may be shared across threads.** Several threads may call a
  `Context`, `Dataset`, `Image` or `CachePool` at once, including `Read`,
  `Prefetch`, `OpenImage`, reductions and histograms on one handle.
- **Concurrent reductions take turns.** Reductions running at once share the
  context's worker pool, which computes one read at a time, so they interleave
  read by read rather than run side by side.
- **Do not assign to a handle another thread is using.** A handle is a value
  holding a shared pointer: give each thread its own copy, or synchronise
  assignment to one variable yourself.
- **Call `ReadAhead::Served` from one thread at a time,** the one that serves
  the frames, even through copies of the controller. `Cancel` and `stats` may be
  called from any thread.
- **Callbacks run on known threads.** Progress callbacks, reduction sinks and
  `cancellation_requested` run on the thread that called the operation. A read
  started by `ReadAhead` runs on the controller's own thread, so a
  `cancellation_requested` given to `ReadAhead` must be safe to call from
  another thread, for example by reading a `std::atomic<bool>`.
- **Your data stays yours to guard.** Do not change a destination buffer or
  region mask while an operation uses it, and copy anything you keep from a
  borrowed view before its callback returns.
