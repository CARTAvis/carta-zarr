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

- **Memory:** each worker holds about three times what it writes at once -- `--block-mib`, or one
  whole shard when shards are larger than that. Lower `--workers` for large shards.
- **Python:** the script runs on uv's own Python build, never the system's. Ubuntu 24.04's 3.12.3
  segfaults as soon as zarr starts its event-loop thread.

Each dataset carries `bench-manifest.json`: its source, layout, striping as the filesystem reports
it, file count and compression ratio. The manifest is written last, so a dataset without one is
unfinished. The last line on stdout is the dataset's path.

`ctest` checks the generator against the library when configured with
`-DCARTA_ZARR_BUILD_BENCH=ON` and `uv` is on the path.

## carta-zarr-bench

Measures how fast the library reads one dataset the way CARTA reads it. Built with
`-DCARTA_ZARR_BUILD_BENCH=ON`, as `bench/carta-zarr-bench` in the build tree; measure with a Release
build, since a Debug one is five to seven times slower and says so on every run.

```sh
carta-zarr-bench probe /lustre/scratch/zarr-bench/<dataset>
carta-zarr-bench run /lustre/scratch/zarr-bench/<dataset> --processes 8 \
    --io-threads 16 --decode-threads 8 --cache-bytes 2G --csv results.csv --resume
```

`probe` prints what the library sees of a dataset as one line of JSON, and fails when it would not
open. `run` writes one CSV row per operation; `--help` lists its options.

- **Modes.** `plane` reads a whole plane at a random channel, `spectrum` every channel at a random
  pixel, `region` reduces every statistic over a box covering 5% of the plane, `cube-histogram` bins
  the cube, and `open` times `Context::Create`, `Dataset::Open` and `OpenImage` together.
- **Positions** come from `--seed` and the trial number and the cube's shape, never its layout, so
  every layout of a cube is read at the same places. They are not aligned to chunks.
- **Processes.** `--processes N` stands for N users, since carta-controller starts one backend per
  user: N processes, each with its own context, released together. They read different positions;
  a cube histogram splits the channels between them. When there are not enough positions to go
  round, the rows that share one say `overlap`.
- **Trials.** Each trial empties the caches, then forks fresh processes, so no trial inherits a cache
  or a thread pool from the one before. Within a trial a process keeps its context, as a backend
  would: later operations may find chunks an earlier one decoded, and an `open` after the first
  finds the metadata in the page cache.
- **Cold reads.** `--cold auto` uses `--drop-cache-cmd` when given, then `drop_caches` as root, then
  `posix_fadvise` on every file of the dataset. Only the command can reach a parallel filesystem's
  servers; the other two empty this host's page cache alone, and `cold_method` says which was used.
  On macOS only the command is available.
- **Deadline.** `--trial-timeout` becomes each read's `ReadControl::deadline`. An operation it stops,
  and any not yet started, is a `timeout` row; a process still running 30 s past it is killed.
- **Resuming.** Rows are written a trial at a time, so an interrupted run loses at most the trial it
  was in. `--resume` skips the trials the CSV already holds for the same settings -- the `run_key`
  column -- that finished without an error.
- **`checksum`** fingerprints what each operation returned, taken after the clock stops: the pixels
  for a read, and the pixel counts and extremes for a reduction or a histogram, which every layout of
  the same pixels must agree on exactly. Sums and histogram counts are left out, because their
  rounding follows the order the chunks were visited in.
- **`storage_read_bytes`** is what `/proc/self/io` says the operation fetched from storage, so it is
  empty off Linux.
