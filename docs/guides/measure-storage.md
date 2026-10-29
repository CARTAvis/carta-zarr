# Measure your own storage

[Storage tuning](../storage-tuning.md) reports what layouts and reader settings
did on the machines measured so far. This page takes a site administrator from a
clone of the repository to a `summary.md` that names the layout and the
carta-backend settings for their own file system. Words the tools and the report
use are defined under [Terms](#terms); every option is in the
[benchmark tools reference](../reference/zarr-bench.md).

Three tools do the work, and you run only the last:

- `generate.py` writes an image dataset in a given layout, from a real cube or a
  synthetic one.
- `carta-zarr-bench` reads a dataset the way CARTA does and times every
  operation.
- `sweep.py` runs the other two over every layout and setting you list, and
  writes `summary.md`.

## Before you start

You need:

- **A C++ toolchain** to build carta-zarr, as the
  [installation guide](../installation.md) describes.
- **[`uv`](https://docs.astral.sh/uv/)** on the path. `generate.py` runs under
  it, on uv's own Python, so nothing else has to be installed.
- **Space on the file system you are measuring** for one dataset uncompressed,
  and a tenth more. The sweep writes one layout at a time and deletes it once it
  is measured. It checks the file system's free space, not your quota.
- **A way to empty caches** before each trial, or every read after the first
  is served from memory and the sweep measures RAM rather than storage. Which way
  depends on where you run it:

  | Where | What empties the caches | Root? |
  |---|---|---|
  | Linux, local disk | `cold = "auto"` writes to `drop_caches` | yes |
  | Linux, local disk | `cold = "auto"` falls back to `posix_fadvise` on every file | no |
  | Lustre | `drop-lustre-cache.py`, client and servers | no |
  | BeeGFS | `posix_fadvise` on this client only; the servers keep theirs | no |
  | macOS | `/usr/sbin/purge`, run as a command | yes |

  Without any of these the sweep stops before it starts, unless every dataset is
  larger than RAM; see [Troubleshooting](#troubleshooting).

## Build the bench

Build in Release with the bench enabled. A Debug build reads five to seven times
slower, and the report says so in its first lines.

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DCARTA_ZARR_BUILD_BENCH=ON
cmake --build build-release --target carta_zarr_bench -j 8
```

The bench is `build-release/bench/carta-zarr-bench`.

## Write a configuration

A first sweep should finish in minutes, so that a mistake costs little. This one
writes a synthetic cube of 256 channels of 1024 × 1024 pixels in two layouts, 1
GiB each uncompressed:

```toml
[paths]
bench = "/path/to/carta-zarr/build-release/bench/carta-zarr-bench"
work = "/lustre/scratch/zarr-bench"   # on the file system you are measuring

[source]
synthetic = "frequency=256,polarization=1,l=1024,m=1024"

[users]
target = 2    # people using CARTA at once at the busiest time
typical = 1   # and on an ordinary day

[measure]
trials = 3
cold = "auto" # or "command", with drop_cache_cmd

[stage1]
chunk = ["l=512,m=512,frequency=4", "l=128,m=128,frequency=64"]
shard = [""]
codec = ["zstd:3"]

[stage2]
top = 1
trials = 2
io_threads = [2, 8]
decode_threads = ["cores"]
cache_mib = [1024]
read_budget = ["0"]
```

Everything left out takes the default that
[`example-sweep.toml`](https://github.com/CARTAvis/carta-zarr/blob/dev/tools/zarr-bench/example-sweep.toml)
shows and explains. To empty the caches with a command, as on Lustre without
root or on macOS, set it under `[measure]`:

```toml
cold = "command"
drop_cache_cmd = "/path/to/carta-zarr/tools/zarr-bench/drop-lustre-cache.py /lustre/scratch/zarr-bench"
```

On macOS, `drop_cache_cmd = "sudo -n /usr/sbin/purge"` works when `sudo -v` has
been run first and is kept fresh, for example with
`while true; do sudo -n true; sleep 50; done &` for the length of the sweep.

For the real sweep, start from `example-sweep.toml` instead. Point
`source.path` at a cube the site holds and `crop` it to 50 to 100 GB, since real
pixels compress differently from synthetic ones and what a layout costs depends
on it. List the chunk shapes the site uses now among the layouts, and set the
users to your own numbers.

## Check it

Every command runs from `tools/zarr-bench`:

```sh
./sweep.py my-sweep.toml --dry-run
```

```text
cold reads by: command
stage1: 2 layouts x 2 bench runs = 4 runs
  L01: chunk l=512,m=512,frequency=4, zstd:3
  L02: chunk l=128,m=128,frequency=64, zstd:3
stage2: up to 1 layouts x 2 settings = 2 runs
each dataset: up to 1.0 GiB uncompressed; one exists at a time
```

It says how caches will be emptied, names each layout, and ends with a
`problem:` line for anything that would stop the sweep: a bench that is not
there, caches that cannot be emptied, a stripe that cannot be set. Fix those
before running. Once stage 1 has run, `--dry-run` also estimates how long
stage 2 will take, which is the long one in a real sweep.

## Run it

```sh
./sweep.py my-sweep.toml
```

The sweep above took six minutes on a desktop with a local disk. A real one,
with a crop of 50 GB and a dozen layouts, takes hours to a day, most of it in
stage 2. It runs in five stages, each into one `results.csv`:

1. **stage1** writes each layout, measures it at the baseline settings with one
   user and with `users.target`, and deletes it.
2. **stage2** tries every reader setting on the best layouts of stage 1.
3. **confirm** measures the recommended settings, their nearest rivals and the
   baseline with one user and with `users.typical`.
4. **validate** writes the recommended layout again from a copy larger than
   RAM, and checks that what the smaller copy said still holds. It runs only
   when `validate_crop` or `validate_synthetic` is given.
5. **warm** measures builds with other tuning constants, only when `[warm]`
   lists them. It is reported and never recommended.

If it stops, run the same command again: finished runs are skipped, and a run
cut short resumes from what it wrote. `--report-only` rewrites `summary.md`
from what there is, mid-sweep or after.

## Read summary.md

The output directory, `zarr-bench-results` beside the configuration unless
`--output` names another, holds `summary.md`. This is the start of the one the
sweep above wrote, on a desktop with cold reads emptied by `purge`.

### The warnings

```text
## Conclusion

- **Too few first touches for plane on L02**: with 2 users, 0 of 64 operations read no chunk
  another one read, so plane is ranked on every operation there, cache hits included, and so is
  every layout compared with it. The crop is small for chunks of 128×128×64. 32 planes a trial
  fall in 4 runs of chunks 64 channels deep: crop at least 2048 channels, or set
  measure.ops plane = 2.
```

Warnings come first, and each says what to change. Two more like this one
followed, for region and spectrum on L01. This one says the cube is too shallow
for chunks 64 channels deep: two users reading 16 planes each always land in a
run of chunks someone else read, so none of their reads measured the disk. A
warning of this kind means its mode compared cache hits rather than reads, and
the sweep should be run again on a larger cube before its answer is trusted. For
a synthetic source, "crop" means a larger `synthetic` shape.

Other warnings say that caches were emptied on this host only, that a layout
failed to write or read, that the recommendation gives up a mode, or that the
recommendation did not survive validation.

### The recommendation

````text
**Once the data is in the recommended layout:** io 8, decode 12, cache 1024 MiB, budget library

```
--zarr_io 8 --zarr_decode 12 --zarr_cache 1
```

or in settings.json: `{"zarr_io": 8, "zarr_decode": 12, "zarr_cache": 1}`

**Recommended layout:** L02 (chunk l=128,m=128,frequency=64, zstd:3).
````

These are the flags to give carta-backend, and the layout to convert cubes to.
A sweep of a real cube also gives settings for the data as it is, the `current`
layout, which apply until the cubes are converted. A synthetic source has no
current layout, so this one gives only the second.

### How each mode fared

```text
| mode | best layout | median | current is |
|---|---|---|---|
| plane | L01 | 5.7 ms | – |
| animation | L02 | 198.9 ms | – |
| spectrum | L02 | 4.2 ms | – |
| region | L02 | 22.9 ms | – |
| cube-histogram | L01 | 252.7 ms | – |
| open | L01 | 7.6 ms | – |
```

No layout is best at everything. Chunks 4 channels deep read a plane fastest;
chunks 64 deep read a spectrum seven times faster and a plane 4.6 times slower.
The recommendation is the layout whose slowdowns, mode by mode against each
mode's best, have the smallest weighted geometric mean, the score:

```text
| layout | score | worst mode |
|---|---|---|
| L02 | 1.35× | 4.60× |
| L01 | 1.65× | 6.91× |
```

`worst mode` is the largest of those slowdowns, which is what users of that
mode will notice. If planes matter more than spectra at your site, raise
`weights.plane`; the report's "Plane against spectrum" section says when
doubling a weight would change the answer:

```text
Doubling weights.plane makes L01 the best layout instead of L02: the choice turns on how much
plane matters against spectrum.
```

### The settings

```text
| settings | score | worst | plane | animation | spectrum | region | cube-histogram | open |
|---|---|---|---|---|---|---|---|---|
| io 8, decode 12, cache 1024 MiB, budget library | 1.03× | 1.19× | 27.7 ms | 198.9 ms | 5.6 ms | 26.8 ms | 285.5 ms | 7.3 ms |
| io 2, decode 12, cache 1024 MiB, budget library | 1.07× | 1.26× | 28.1 ms | 199.0 ms | 5.7 ms | 22.6 ms | 331.6 ms | 9.2 ms |
```

Stage 2 ranks the reader settings on the best layouts the same way. Here they
differ by a few per cent, which on a local disk is near the noise; on a parallel
file system the number of reads in flight, `io`, usually matters more. The
section after it, "Sensitivity to the number of users", says what the
settings chosen for the peak cost one user alone and the typical number.

The rest of the report is the tables behind these: every mode of every layout,
with medians, 90th percentiles and first-touch counts, and the configuration as
it ran, every default filled in.

## Give carta-backend the settings

Pass the flags on carta-backend's command line, or put the `settings.json` line
in its settings file. carta-controller starts one backend per user, so each user
gets them, and a server with sixteen users at once holds sixteen caches of
`--zarr_cache`. Write cubes in the recommended layout when they are next
converted to Zarr. [Storage tuning](../storage-tuning.md) says what each setting
does.

## Terms

- **Layout:** How a cube is stored: its chunk shape, shards, codec, and on Lustre
  or BeeGFS its striping. The sweep names them `L01`, `L02`, … in the order of
  `[stage1]`'s grid.
- **Modes:** what the bench times, each the way CARTA does it.
  - `plane` reads a whole plane at a random channel.
  - `animation` plays 32 consecutive planes at 5 frames a second, as CARTA's
    animator does. Its time is per frame: the 200 ms a frame's turn takes at
    that rate, plus however late the frame was. About 200 ms means no frame
    stalled; it is not how long a read takes.
  - `spectrum` reads every channel at a random pixel.
  - `region` computes every statistic over a box covering 5% of the plane, for
    every channel.
  - `cube-histogram` bins the whole cube, in two passes as carta-backend does.
  - `open` opens the dataset and its image.
- **current:** The layout of the source cube as it is. Every other layout is
  compared with it. A synthetic source has none.
- **User:** One `carta-zarr-bench` process, standing for one carta-backend:
  carta-controller starts one per user. `users.target` is how many run at once.
- **Cold read:** A read with the caches emptied first, so that it reaches
  storage. `cold reads` in the report says how they were emptied; "on this host
  only" means the servers of a parallel file system kept theirs.
- **First touch:** An operation that reads no chunk another operation of the same
  trial read before it or alongside it, so its time is the storage's rather than
  a cache's. `plane`, `spectrum` and `region` are ranked on first touches.
- **Shares chunks:** The opposite of a first touch. The bench marks such an
  operation in its `shares_chunks` column.
- **Run of chunks:** The chunks one plane decodes, a chunk deep along the
  spectrum, and so every plane that decodes the same ones. An animation is
  served from the cache within a run and stalls on the frame that enters the
  next.
- **Baseline:** The settings carta-backend runs with now, `[baseline]`, which
  every recommendation is compared with.
- **Score:** The geometric mean of each mode's slowdown against its own best,
  weighted by `[weights]`. 1.00× is best at everything; lower is better.
- **Front:** The layouts no other layout beats in every mode. Which of them to
  take is what the weights decide.
- **Makespan:** The slowest user's whole trial, from start to finish.
- **Read budget:** How much memory one read may hold, `budget library` being the
  library's own choice. carta-backend has no flag for it, so it is measured and
  never recommended.

## Troubleshooting

- **`problem: caches cannot be emptied here`:** Nothing can empty the page cache
  with the `cold` method chosen, and the dataset fits in RAM, so every trial
  after the first would be read from memory. Give `cold = "command"` and a
  `drop_cache_cmd` that works here, or measure a cube larger than RAM.
  `[checks] ram = false` runs anyway, and the report then says reads were not cold; its
  numbers are memory's.
- **`carta-zarr-bench was not a Release build`:** The bench came from a Debug
  tree. Build it as in [Build the bench](#build-the-bench) and run again into a
  new output directory.
- **Too few first touches:** The cube is shallow, or its plane small, for the
  chunks and the number of users: operations kept landing in chunks another one
  read. Do what the warning says, usually a deeper crop or fewer operations of
  that mode a trial.
- **No room for a crop:** The dry run reports how large each dataset is
  uncompressed. Crop fewer channels, or set `paths.work` on a file system with
  more room. Only one dataset exists at a time unless `keep_datasets` is true.
- **BeeGFS without root:** `beegfs-ctl --setpattern` belongs to root unless the
  file system lets users set it, so `stage1.stripe` may not be sweepable;
  `--dry-run` says so. Nothing an ordinary user can run empties a storage
  server's cache, so the report says reads were cold on this host only, and a
  layout read twice in a row reads faster the second time.
- **Lustre without root:** Use `drop-lustre-cache.py` as the command, as above.
  Without its second step, which asks the servers to drop their cache, a sweep
  measures the servers' memory: on one Lustre 2.15 system a whole-cube read took
  1.2 s warm, 2.0 s with only the client's cache emptied, and 5.5 s with both.
- **`holds runs of another workload`:** The output directory holds a sweep of a
  different source, bench build or measurement setting. Write this one to
  another directory with `--output`.
