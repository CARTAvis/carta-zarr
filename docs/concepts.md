# Datasets, images and reads

## Dataset versus image

An **image dataset** is one Zarr v3 group opened by path. It can contain several
image data variables, shared coordinates, flags and beam fit parameters.
An **image** is one supported data variable opened by its variable name.

```{mermaid}
flowchart LR
    D["Dataset: cube.zarr"] --> S["ImageEntry: SKY"]
    D --> R["ImageEntry: RESIDUAL"]
    D --> C["Shared coordinates"]
    S -->|"OpenImage(id)"| I["Image handle"]
    I -->|"Read(request, buffer)"| P["float32 pixels"]
```

The names above are illustrative: use the IDs the actual listing reports.
An image role such as sky or residual describes a variable's purpose; it is
not the ID used to open it.
For variables typed `sky`, a non-base data group's name can determine the role
(`model`, `residual`, `dirty`). An opened image's `sub_type` instead describes
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

The output is densely packed with logical axis 0 varying fastest. When only
the two spatial axes vary, the x coordinate varies fastest for the current
XRADIO profile. `storage_index` reports where an axis lives in the stored array.

## Flag versus region mask

An image's **flag** marks bad pixels: true means excluded. Reads apply it by
default and output NaN for those pixels. Set `apply_pixel_mask = false` to read
the image's stored values without applying the flag.

A reduction's **region mask** selects pixels within a bounding box: nonzero
means included. A deconvolution mask is an ordinary image, rather than a flag.
These three serve different purposes.

## Results, diagnostics and lifetimes

Public operations return `Result<T>` with an error code, message and node path
on failure. A diagnostic accompanies a successful listing or descriptor and
explains an unsupported candidate or missing optional information.

`Context`, `Dataset`, `Image` and `CachePool` are shared handles created by
factories. Copies share their resources, and moving one preserves the source
handle too. Use `std::optional` when your own object needs a not-yet-created
handle. A reference obtained from `descriptor()` must not outlive the object
that owns the underlying descriptor.

Pixel destinations and region views belong to the caller. Keep them alive
through the operation. Reduction blocks are borrowed only during the sink
callback; copy values you need to retain. On a failed read, discard the
destination even if progress had already reported a completed prefix.
