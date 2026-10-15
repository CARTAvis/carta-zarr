# Reduce regions and compute histograms

This guide follows the statistics example, whose steps are listed under
[what the program does](#what-the-program-does).

## Reduce regions over channels

`PlaneSelection` carries a spectral range, one polarization and one time.
Supply many regions in one request so their shared chunks are visited together.
Each region is a `RegionMask`: a bounding box and an optional `mask` with one
byte per pixel of the box, nonzero to include it. The example's two regions
have no `mask`, so each includes its whole box.

The sink receives `SpectralBlock`s. A block with `complete == false` is an
estimate over the pixels visited so far, and a final block with
`complete == true` replaces it later; the example skips estimates. Returning
false from the sink cancels the operation.

`first_channel` is an index into the requested spectral selection, so the
image channel is `start + (first_channel + i) * stride`, as the program prints.

`sum_sq_dev` is the sum of squared deviations from the mean. For a mean, check
that `num_pixels > 0`; for a sample variance, require more than one pixel and
divide `sum_sq_dev` by `num_pixels - 1`.

A flux density or world-coordinate operation belongs in the application,
which has the scientific context needed to interpret the totals.

## Compute histograms

`ComputeHistogram` uses caller-specified finite bounds and returns a histogram
for each selected plane. This example uses 64 bins over [-1e6, 1e6], so values
outside that interval and non-finite pixels are excluded. Choose bounds suitable
for your own image. `HistogramBlock::Counts` returns a borrowed view valid only
within the sink; the example consumes it immediately with `std::accumulate`.

`ComputeCubeHistogram` instead finds the range and returns one histogram for
all selected planes. Its provisional bins are re-aggregated onto the final
range, so individual bin counts are approximate and can depend on thread count.
The example checks for finite pixels before printing extrema. It does not
sample; `spatial_sample = n` takes every nth pixel along both spatial axes and
sets `sampled`, but saves decoding only when n is large enough to skip whole
chunks.

Both histograms cover **whole planes**, not just the region reduction's patch.
The [rules shared by every operation](../concepts.md#rules-that-apply-to-every-operation)
apply to all three: flags, borrowed views and checking each `Result`.

## What the program does

1. Opens the first openable image and bounds a patch of at most 64 × 64 pixels.
2. Reduces that patch and a single-pixel region together over up to four
   channels, printing each channel's count, mean and variance.
3. Computes a 64-bin histogram for each of those channels over fixed bounds.
4. Computes one cube histogram across all of them and prints its range.
5. Returns zero only if the reduction and the plane histograms delivered a
   final block for every requested channel, and the cube histogram has all
   its bins.

## Complete source

Build and run with the [example commands](../examples.md#build-and-run).
The [shared setup](../examples.md#shared-setup) keeps image selection and shape
checks consistent across the advanced examples.

```{literalinclude} ../../examples/statistics.cc
:language: cpp
```

See the [reduction API](../api/reduce.md) for callbacks, limits and approximation details.
