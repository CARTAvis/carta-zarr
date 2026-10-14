# Reduce regions and compute histograms

The [runnable statistics example](../examples.md) opens an image and works on
up to four channels at the first polarization and time. It reduces a bounded
spatial patch and a single-pixel region together, then computes per-plane and
combined cube histograms.

## Statistics for regions over channels

`PlaneSelection` carries a spectral range, one polarization and one time.
Supply many regions in one request so their shared chunks are visited together.
The example uses two unmasked `RegionMask` rectangles; a nonzero raster value
would include a pixel, while an absent raster includes the whole rectangle.
Keep the region array and any rasters alive until the operation returns.

Consume a `SpectralBlock` only during the callback. The example ignores blocks
where `complete` is false: incomplete blocks are refined and emitted again,
and final values replace earlier estimates. Copy data if you need to retain it
past the callback. Returning false cancels the operation.

`first_channel` is an index into the requested spectral selection. The actual
image channel is `start + (first_channel + i) * stride`, as the program prints.
For a mean, check that `num_pixels > 0`; for sample variance, require more than
one pixel and divide `sum_sq_dev` by `num_pixels - 1`.

A flux density or world-coordinate operation belongs in the application,
which has the scientific context needed to interpret the totals.

## Histograms

`ComputeHistogram` uses caller-specified finite bounds and returns a histogram
for each selected plane. This example uses 64 bins over [-1e6, 1e6], so values
outside that interval and non-finite pixels are excluded. Choose bounds suitable
for your own image. `HistogramBlock::Counts` returns a borrowed view valid only
within the sink; the example consumes it immediately with `std::accumulate`.

`ComputeCubeHistogram` instead finds the range and returns one histogram for
all selected planes. Its provisional bins are re-aggregated onto the final
range, so individual bin counts are approximate and can depend on thread count.
The example checks for finite pixels before printing extrema. Setting
`spatial_sample > 1` samples spatial pixels and sets `sampled`; sampling within
a chunk may save no decompression.

Both histograms cover **whole planes**, not just the region reduction's patch.
All three operations apply the image's flag by default and check their `Result`.

## Complete source

Build and run with the [example commands](../examples.md#build-and-run).
The [shared setup](../examples.md#shared-setup) keeps image selection and shape
checks consistent across the advanced examples.

```{literalinclude} ../../examples/statistics.cc
:language: cpp
```

See the [reduction API](../api/reduce.md) for callbacks, limits and approximation details.
