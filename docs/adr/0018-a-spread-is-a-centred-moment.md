# A spread is a centred moment, not a difference of sums

A spectral reduction and a cube histogram report the sum of squared deviations of the pixels from
their own mean, `Statistic::sum_sq_dev`, and a consumer makes its standard deviation from that:
`sqrt(sum_sq_dev / (num_pixels − 1))`. It does not make it from `sum` and `sum_sq`, although it
could, and every consumer of this library did until now.

## Why not from the sums

`sqrt((sum_sq − sum² / n) / (n − 1))` is right in exact arithmetic and wrong in doubles once the
pixels are far from zero against their spread. Both terms of the subtraction are then about
`n · mean²`, they agree in every digit a double holds, and what is left of the subtraction is the
rounding each of them carried. A million float32 pixels, accumulated in doubles:

| pixels | from the sums | centred | the pixels' own |
|---|---|---|---|
| all 1e8 | NaN | 0 | 0 |
| 1e6 + noise of σ 0.5 | 1.918 | 0.4996 | 0.4996 |
| 1e7 + noise of σ 0.5 | 0.932 | 0.5698 | 0.5698 |
| 1e8 + noise of σ 20 | NaN | 20.14 | 20.14 |

(The third row's own answer is 0.57, not 0.5, because float32 holds 1e7 only to the nearest unit.)
NaN is the subtraction coming out below zero. Clamping it to zero, which a consumer did for a day,
turns the first row right and the fourth into a confident zero; and nothing at all can be done about
the second and third, which look plausible. The sums hold no more than that. Only accumulating
something else does.

## How it is accumulated

Updating a running mean per pixel -- Welford's method -- costs a division per pixel and stops the
loop vectorising, and the per-pixel loop is what a reduction's time is (ADR 0005). So it is made in
two steps (`src/reduce/deviations.h`):

- **Within a span**, the loop sums the pixels' distances from a shift, the span's first finite
  pixel, and their squares: one subtraction and one multiply-add more per pixel. A span's pixels
  are near one of their own against their spread, so turning those two sums into the span's
  deviations subtracts nothing large. In the worst case -- the shift an outlier -- a 512-pixel span
  loses about three digits, against the sixteen the difference of sums can lose.
- **Above a span** -- the spans of a row, the rows of a region, a task's partials, a block's reads,
  a cube histogram's accumulators -- two sets are put together by the formula of Chan, Golub and
  LeVeque, from their counts, their means and their deviations, in the order `MergeInOrder` already
  fixes for the sums. It too subtracts nothing large.

The mean those merges need is kept, not had from `sum / num_pixels`, and kept as a base and an
offset from it. A running sum near 1e7 has already lost the digits a mean needs, and so has a mean
held in one double, while the merge multiplies the error in a difference of two means by the size of
the sets. On a 256 × 260 plane of 1e7 with a spread of 0.5, `sum_sq_dev` was 1.3e-10 out with the
mean had from the sum, 5e-11 with it held as one double, and 7e-16 as a base and an offset. The base
is the first span's shift, a pixel, so two bases are near each other and their difference is exact.
The base and the offset are two values per region and channel that a block keeps and never reports.

The loop takes distances only when `sum_sq_dev` was asked for; otherwise it is the loop it was, an
instantiation of its own -- and so is everything under one unit of the walk, the folds included,
since asking the question once per row cost a reduction of one-pixel strips 2% whether or not it
wanted the spread. Asking for it brings `num_pixels` and `sum` with it, since a merge needs both, and
a block reports them. A cube histogram always takes them: its loop is bound by the branch on
finiteness and the histogram it writes, not by arithmetic.

What it costs, measured with `carta_zarr_pass_timing` on the Apple-silicon Mac against the ASKAP cube
(aligned builds, six interleaved rounds, median of seven repeats each, paired by round):

| entry | 1 thread | 4 threads |
|---|---|---|
| ReduceSpectral, two regions, with `sum_sq_dev` | +5.8% | +6.2% |
| ReduceSpectral, the whole plane, with it | +5.5% | +6.8% |
| ReduceSpectral, 64 one-pixel strips, with it | +3.9% | +4.3% |
| ReduceSpectral, a masked ellipse, with it | +6.6% | +7.0% |
| ComputeCubeHistogram, which always counts it | +3.8% | +3.4% |
| every entry not asked for it | within ±0.4%, each round's difference either side of zero | |

Nothing else moves. The sum, the sum of squares, the extrema, the counts and every histogram bin are
the same bits they were, with or without `sum_sq_dev` in the request, which an A/B trace of every
answer over two fixtures and a 1024² × 32 synthetic cube, at one thread and four, checks.

## What it does not reach

The variance is the consumer's to derive, and so is what one pixel's is: the library reports the
sum, as it reports `sum` and not a mean.

Statistics a consumer did not have the library count -- its own over pixels it read some other way,
or a file's precomputed ones, which hold only the sum and the sum of squares -- cannot be helped from
here. CARTA's backend makes sigma from the spread for a Zarr image and leaves every other format's
made from the sums as it always was, so the same pixels stored two ways can report two sigmas when
they are far from zero against their spread.

## Consequences

`SpectralTotals` gained a field, which changes the library's ABI. The version stays 0.1.0: nothing
has been released, and the one consumer is built against each change as it lands. What a version
promises before the first release, and when its soname changes, is decided separately.
