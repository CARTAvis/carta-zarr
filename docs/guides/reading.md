# List images and read a channel or region

The snippets here assume checked `Context`, `Dataset` and `Image` handles as
created in the [complete example](../quickstart.md), and use namespace `carta::zarr`.
The [runnable read-control example](resources.md#complete-read-control-example)
selects the second channel when available and performs a bounded spatial read
with progress and error handling.

## Inspect the listing

```cpp
for (const auto& entry : dataset.descriptor().images) {
    std::cout << entry.id << " openable=" << entry.openable << '\n';
    for (const auto& axis : entry.axes) {
        std::cout << "  " << axis.name << ": " << axis.length << '\n';
    }
    for (const auto& diagnostic : entry.diagnostics) {
        std::cout << "  " << diagnostic.message << '\n';
    }
}
```

Choose an `openable` entry and pass its `id` to `OpenImage`. A
`default_image_id`, when present, is a suggested choice rather than the only
image. Flags and coordinates do not need to be opened as images to read a sky
variable.

## Select a channel

Start with one element on every axis, then expand the spatial axes as in the
quickstart. Select a spectral channel by role:

```cpp
const auto& axes = image.descriptor().axes;
const auto spectral = AxisIndex(axes, AxisRole::spectral);
if (!spectral || channel >= axes[*spectral].length) {
    // Report an invalid channel to your application's caller.
    return;
}
request.axes[*spectral] = {channel, 1, 1};
```

Select polarization and time the same way. Their default indices are zero;
a spectral channel's index is distinct from its frequency coordinate value.

## Select a spatial rectangle

Set the spatial ranges to `{x_start, width, 1}` and `{y_start, height, 1}`.
Ensure the starts and counts stay within the axes. A stride of two selects
every second sample; `count` is the number of selected samples, not the width
of the interval they span.

The destination needs the product of all selected counts in elements. Check
for multiplication overflow and allocation limits before constructing a buffer.
The quickstart bounds its rectangle and checks this product.

## Flags and non-finite pixels

Flags are applied by default (see the
[core concepts](../concepts.md#rules-that-apply-to-every-operation)). To read
the stored values instead:

```cpp
ReadOptions options;
options.apply_pixel_mask = false;  // Do not apply the image's flag.
auto result = image.Read(request, {pixels.data(), pixels.size()}, options);
if (!result) {
    // Report result.error() and discard pixels.
    return;
}
```

Stored NaNs remain NaNs even with flag application disabled. Inspect pixels
with `std::isfinite` when calculating your own statistics. Reductions exclude
non-finite values and report their count separately.

See the [read API](../api/read.md) for range, buffer and option contracts.
