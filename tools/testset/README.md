# Read-path test set

A fixed set of cubes to measure carta-zarr's read paths on, and to compare them with FITS: each cube
is a FITS file, and each of its Zarr layouts is made from that FITS by the site's xradio converter
(`fits_to_zarr_xradio_v1.2.2_v1.py`), as production data is. Every cube is synthetic, so the whole
set rebuilds from nothing on any machine.

```sh
XRADIO_CONVERTER=/path/to/fits_to_zarr_xradio_v1.2.2_v1.py tools/testset/build.sh /data/carta-testset
```

`build.sh` skips what already exists, so a run that stopped is finished by running it again. A FITS
file or a Zarr appears under its own name only once it is complete, and a Zarr only once `verify.py`
has found every one of its pixels, its flag (or the absence of one) and its coordinates to match
the FITS file's (the report is left beside it as `NAME.verify.txt`).
It needs `uv`, the converter, and about 155 GB (190 GB while a cube is being synthesized); on a
28-thread machine with local NVMe it takes about 13 minutes, a third of it verifying.

## Cubes and layouts

Shapes are written l x m x frequency; every other axis is 1.

| cube | shape | FITS |
|---|---|---|
| pancake | 7763 x 4742 x 256 | `pancake.fits` (36 GB) |
| cigar | 512 x 512 x 30000 | `cigar.fits` (30 GB) |

| Zarr | chunk | shard | flag |
|---|---|---|---|
| `pancake_c256x256x16.zarr` | 256 x 256 x 16 (4 MiB) | - | - |
| `pancake_c256x256x16_flag.zarr` | 256 x 256 x 16 (4 MiB) | - | yes |
| `pancake_c512x512x4.zarr` | 512 x 512 x 4 (4 MiB) | - | - |
| `cigar_c128x128x64.zarr` | 128 x 128 x 64 (4 MiB) | - | - |
| `cigar_c64x64x256.zarr` | 64 x 64 x 256 (4 MiB) | - | - |

Each Zarr is about 21 GB.

- A pancake has large planes and few channels; a cigar small planes and many. Both are chunked at
  4 MiB, as sites have chunked their cubes.
- The pancake is chunked two ways at that size. 256 x 256 x 16 is what sites have used; 512 x 512 x 4
  is a shallower chunk, which decodes 4 channels to read one where the other decodes 16, at the cost
  of four times as many chunks along a spectrum.
- The cigar is chunked two ways the other way round. 128 x 128 x 64 is shaped as sites have shaped
  their chunks; 64 x 64 x 256 is a deeper chunk, nearer the depth `docs/storage-tuning.md` points to
  for a cube this long, which reads a spectrum from a quarter as many chunks at the cost of decoding
  256 channels to read one.
- A layout is a chunk, an optional shard (a whole number of chunks) and whether it carries a flag.
  With a flag the converter writes `FLAG_SKY`, true where the pixel is NaN (`--compute_mask`); the
  pixels are the same, so the two differ only in what reading the flag costs. Add a layout by adding
  a line to `layouts` in `build.sh`; its name says all three, e.g. `pancake_c256x256x16_s1024x1024x16_flag.zarr`.
- xradio 1.2.2's FITS reader leaves that flag undeclared: untyped, not named by the image or a data
  group, so carta-zarr does not apply it. The converter has to declare it (type `flag`, the image's
  `flag` attribute, and a `base` data group naming it); `verify.py` refuses a flag carta-zarr would
  not apply, so a converter that does not makes `build.sh` stop at the first flagged layout.
- The converter writes no chunk that is wholly NaN, so a fully flagged run of channels long enough
  to cover a chunk leaves chunks that are not on disk at all.

## What the cubes hold

`tools/zarr-bench/generate.py --synthetic` with its defaults, which imitate a continuum-subtracted
HI cube from an ASKAP mosaic:

- ASKAP's frequency axis: 1295.5 MHz plus 18.5 kHz a channel, rest frequency HI's. The cigar's
  30000 channels span 556 MHz.
- Noise of about 2 mJy/beam whose rms varies by about 15 % from channel to channel and doubles
  towards the footprint's edge, as primary-beam correction makes it, with a heavy tail: 4 % of
  pixels drawn six times as wide.
- Faint line sources, one per 10^8 pixels (94 in the pancake, 79 in the cigar): an elliptical
  Gaussian on the sky with a Gaussian or double-horned profile 50 to 500 km/s wide, peaking at 1 to
  20 times the noise.
- NaN outside an irregular rounded footprint holding 76 % of each plane.
- Runs of flagged channels, flagged whole or on one side of a line across the plane, 2 % of the
  channels; the cigar's happen to include a run covering 64 whole channels. The pancake's runs are
  too short for that, so channels 128 to 192 are flagged whole as well (`--flagged-range`), covering
  whole chunks.

`zarr-to-fits.py` writes the FITS file, with the axes ASKAPsoft writes (RA, Dec, Stokes, frequency)
and the header keywords xradio's FITS reader requires.

## Calibration

`stats.py` measures a dataset the way the synthetic cube was calibrated. The real cube is ASKAP
Hydra (`askap_hydra_extragalactic_256`), converted the same way at 256 x 256 x 16:

| | ASKAP Hydra | synthetic pancake |
|---|---|---|
| compression | 1.39 | 1.25 |
| NaN in a plane | 24.1 % | 24.0 % |
| noise (robust rms) | 2.12 mJy/beam | 2.05 mJy/beam |
| beyond +5 rms | 0.863 % | 0.858 % |
| beyond -5 rms | 0.862 % | 0.855 % |
| channels wholly NaN | 0 | 66 |

The synthetic cube compresses about 10 % worse. Its noise is independent from pixel to pixel, where
a restored image's is correlated over the beam (ASKAP Hydra: 0.90 between neighbouring pixels, 0.66
two apart), and neighbours that are alike shuffle and compress better. That is left as it is: the
test set is for comparing layouts and formats on one set of pixels, not for passing as real data.

The synthetic cube flags more than ASKAP Hydra does, on purpose: a run that covers whole chunks
leaves them NaN, or not on disk at all, and that is a read path real cubes take too (channels 128 to
192 of the pancake, and one run of the cigar, 16 chunks). The cigar measures the same per plane, and
compresses 1.51 times; there is no real cigar to set it against.

The noise is the robust rms (1.4826 x MAD) of a plane's finite pixels; the tails are the fraction
of those pixels beyond +/-5 of it, the NaN outside the footprint not counted. Compression is blosc zstd 5 with shuffle, decoded bytes over bytes on disk of
the chunks on disk, so that chunks a writer left out for being all NaN do not count.
