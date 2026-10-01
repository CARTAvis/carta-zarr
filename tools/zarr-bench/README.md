# zarr-bench

Tools for finding which Zarr layout, and which carta-zarr settings, read fastest on a given storage
system -- Lustre and BeeGFS in particular, where the answer differs from a local disk.

## generate.py

Writes one XRADIO image dataset in a chosen layout. It is a self-contained
[`uv`](https://docs.astral.sh/uv/) script, so it needs nothing installed but `uv`.

Rewrite a real cube, cropped to something that can be rewritten many times over:

```sh
./generate.py --source /data/cube.zarr --crop frequency=0:2000 \
    --chunk l=512,m=512,frequency=16 --shard frequency=256 --codec blosc:zstd:5:shuffle \
    --output-root /lustre/scratch/zarr-bench
```

Or synthesize one, when there is no real cube to hand:

```sh
./generate.py --synthetic --shape frequency=4000,polarization=1,l=4096,m=4096 \
    --chunk l=256,m=256,frequency=64 --codec zstd:3 --stripe lustre:count=4,size=4M \
    --output-root /lustre/scratch/zarr-bench
```

- **Axes** are named as XRADIO names them -- `time`, `frequency`, `polarization`, `l`, `m` -- and an
  axis left out of `--chunk` is 1. An axis left out of `--shard` is its chunk.
- **Codecs:** `none`, `zstd[:level]`, `gzip[:level]`, `blosc[:cname[:level[:shuffle]]]`.
  `--keep-bits N` rounds pixels to N mantissa bits before compressing them.
- **Real pixels are the better source.** What a layout costs depends on how well the pixels
  compress, and synthetic noise does not compress like a real image. The synthetic cube is
  deterministic: two layouts of the same `--seed` and `--shape` hold the same pixels.
- **`--stripe`** sets Lustre (`lfs setstripe`) or BeeGFS (`beegfs-ctl --setpattern`, or
  `beegfs entry set` on BeeGFS 8) striping on the output directory before anything is written,
  since striping applies only to files created after it. BeeGFS normally reserves this for root;
  the generator stops rather than writing an unstriped dataset.
- **`--output-root`** names the dataset after a hash of everything that decides its bytes, and
  reuses one that is already complete. `--output` names it yourself.

Each dataset carries `bench-manifest.json`: its source, layout, striping as the filesystem reports
it, file count and compression ratio. The manifest is written last, so a dataset without one is
unfinished. The last line on stdout is the dataset's path.

`ctest` checks the generator against the library when configured with
`-DCARTA_ZARR_BUILD_BENCH=ON` and `uv` is on the path.
