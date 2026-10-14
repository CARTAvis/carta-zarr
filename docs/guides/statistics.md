# Reduce regions and compute histograms

All operations below act on a `PlaneSelection`: a spectral range plus one
polarization and one time. These snippets assume an open `Image`, valid axis
lengths and the namespace `carta::zarr`.

## Statistics for a region over channels

```cpp
RegionMask region{0, 0, width, height};  // Entire rectangle; no raster mask.
SpectralReduceRequest request;
request.planes.spectral = {0, channel_count, 1};
request.regions = {&region, 1};
request.statistics = Statistic::num_pixels | Statistic::sum | Statistic::sum_sq_dev;

auto result = image.ReduceSpectral(request, [](const SpectralBlock& block) {
    if (!block.complete) {
        return true;  // This consumer only uses final results.
    }
    for (std::uint64_t i = 0; i < block.channel_count; ++i) {
        const auto totals = block.Totals(0, i);
        if (totals.num_pixels > 0) {
            std::cout << block.first_channel + i << ": "
                      << totals.sum / totals.num_pixels << '\n';
        }
    }
    return true;
});
```

Check `result` for failure. Region dimensions must fit within the spatial axes.
Keep the region array and any raster masks alive until the call returns. Supply
many regions in one request so shared chunks are visited together. A nonzero
raster value includes the corresponding pixel; an absent raster includes the
whole rectangle.

`first_channel` is an index into the requested selection. Its image channel is
`start + (first_channel + i) * stride`. Incomplete blocks are refined and emitted
again; replace earlier estimates rather than adding them to final totals.
Pointers returned by `Series` are valid only within the callback.

For a sample variance use `sum_sq_dev / (num_pixels - 1)` only when the count
exceeds one. A flux density or a world-coordinate operation belongs in the
application, which has the required scientific context.

## Per-plane histogram with known bounds

```cpp
HistogramRequest request;
request.planes.spectral = {0, channel_count, 1};
request.bins = 256;
request.lower = -1.0;
request.upper = 1.0;
auto result = image.ComputeHistogram(request, [](const HistogramBlock& block) {
    if (block.complete) {
        for (std::uint64_t channel = 0; channel < block.channel_count; ++channel) {
            const auto* counts = block.Counts(channel);
            // Consume or copy block.bin_count values here.
            (void)counts;
        }
    }
    return true;
});
```

Bounds must be finite, ordered and representable as float. Values outside the
range and non-finite values are excluded. Check the operation's result.

## One cube histogram without known bounds

```cpp
CubeHistogramRequest request;
request.planes.spectral = {0, channel_count, 1};
request.bins = 256;
auto result = image.ComputeCubeHistogram(request);
if (result) {
    const auto& histogram = result.value();
    // counts span [totals.min, totals.max]. Check num_pixels before using extrema.
}
```

This returns one histogram over the selected planes. Its provisional bins are
re-aggregated onto the final range, so individual bin counts are approximate
and can depend on thread count. Extrema and totals are accumulated from the
selected finite pixels. `spatial_sample > 1` samples spatial pixels and sets
`sampled`; sampling within a chunk may save no decompression.

See the [reduction API](../api/reduce.md) for callbacks, limits and approximation details.
