# Tuning Zarr storage for CARTA on Lustre and BeeGFS

For administrators of a CARTA server whose image cubes are XRADIO Zarr datasets on a parallel file
system. It says what can be tuned, what measurements on two file systems found, and how to measure
your own. The tools themselves are described in
[`tools/zarr-bench/README.md`](../tools/zarr-bench/README.md).

## In short

- **Choose the chunk depth by the size of the plane.** The
  [rule of thumb](#the-best-chunk-depth-depends-on-the-size-of-the-plane) gives a starting point:
  about 1 to 4 channels deep for a 7763 x 4742 plane, 16 to 32 for 2048 x 2048 with 2048 channels,
  and over 100 for 512 x 512 with thousands of channels.
- **Start from 4 MiB chunks on Lustre; on a local disk 1 MiB are as good.** Depth matters more than
  size: [chunk size matters less than depth](#chunk-size-matters-less-than-depth).
- **Set `--zarr_file_io_threads 8` on Lustre**, rather than the default of 2:
  [more file-reading threads](#use-more-file-reading-threads-on-lustre).
- **For animations, let the cache hold two runs of chunks**, a run being the plane's area times the
  chunk depth times four bytes. With less, the backend reads nothing ahead and stalls where each run
  begins: [animations](#animations-stall-where-a-run-of-chunks-begins).
- **On a busy parallel file system, gather small chunks into shards** to cut the number of files:
  [striping and shards](#striping-and-shards).
- **Measure your own system** before rewriting an archive: these numbers describe the machines they
  came from. [Measuring your own system](#measuring-your-own-system) says how.

## What can be tuned

Two things, by two different people at two different times:

- **The layout of the data**, chosen when a cube is converted to Zarr: the chunk shape, whether
  chunks are grouped into shards, the compressor, and on Lustre the striping of the directory it is
  written to. Changing it means rewriting the cube.
- **carta-backend's reader settings**, chosen when the server is configured, as flags or in
  `settings.json`. They belong to carta-backend's Zarr support, which is not yet in any of its
  releases:

  | flag | default | what it sets |
  |---|---|---|
  | `--zarr_file_io_threads` | 2 | how many reads of chunk files run at once |
  | `--zarr_data_copy_threads` | every logical core | how many threads decode chunks |
  | `--zarr_cache_size` | 1024 (MiB) | the decoded-chunk cache each backend process keeps |
  | `--zarr_histogram_method` | `exact` | how a cube histogram is computed |

  carta-controller starts one backend process per user, so every user gets these, and a server with
  sixteen users at once has sixteen caches.

The constants inside carta-zarr's reductions are not settings; [ADR 0014](adr/0014-tuning-constants-are-measured-not-configured.md)
explains why, and how they are measured.

## What the measurements found

On a server (two Xeon Gold 6134, 32 logical cores, 1 TB RAM) against its Lustre 2.15 and BeeGFS file
systems, with eight users reading at once and every cache emptied before each trial, and on a desktop
with 28 logical cores against a local NVMe disk for comparison. The numbers below are medians per
operation; a plane read is what jumping to another channel costs, since carta-backend reads the whole
plane of a Zarr cube into memory and downsamples it there.

### The best chunk depth depends on the size of the plane

A chunk deep in channels holds a long run of each pixel's spectrum, so a spectrum reads few chunks;
but a plane then reads every channel of its chunks, not only its own. Which matters more depends on
how large the plane is.

A 512 x 512 x 4096 crop of an ASKAP cube on Lustre:

| chunk (l x m x channels) | plane | spectrum | region (5% of the plane) |
|---|---|---|---|
| 512 x 512 x 1 (the cube's own) | 6 ms | 2.40 s | 6.02 s |
| 512 x 512 x 16 | 37 ms | 2.36 s | 3.54 s |
| 256 x 256 x 32 | 41 ms | 529 ms | 1.38 s |
| 128 x 128 x 128 | 146 ms | 222 ms | 382 ms |

A synthetic 2048 x 2048 x 2048 cube on the same Lustre:

| chunk (l x m x channels) | plane | spectrum | region |
|---|---|---|---|
| 2048 x 2048 x 1 | 205 ms | 14.8 s | 28.1 s |
| 512 x 512 x 16 | 478 ms | 1.55 s | 2.93 s |
| 256 x 256 x 32 | 1.08 s | 838 ms | 4.12 s |
| 128 x 128 x 128 | 2.81 s | 246 ms | 4.04 s |

On the 512-square plane eight users' spectra and regions could not all read chunks of their own --
there are too few chunks across it -- so those columns include some reads the page cache answered,
which if anything flatters the layouts with large chunks across the sky.

A real 7763 x 4742 x 128 ASKAP cube on the same Lustre, in the 512 x 512 x 4 chunks it was written in
and in others:

| chunk (l x m x channels) | size | plane | spectrum | region | score |
|---|---|---|---|---|---|
| 1024 x 1024 x 4 | 16 MiB | 578 ms | 451 ms | 2.50 s | 1.83 |
| 512 x 512 x 16 | 16 MiB | 1.47 s | 119 ms | 1.97 s | 1.89 |
| 512 x 512 x 4 (the cube's own) | 4 MiB | 1.00 s | 232 ms | 2.91 s | 1.94 |
| 256 x 256 x 16 | 4 MiB | 2.39 s | 67 ms | 2.69 s | 2.03 |
| 1024 x 1024 x 1 | 4 MiB | 279 ms | 792 ms | 3.84 s | 2.03 |
| 256 x 256 x 64 | 16 MiB | 4.56 s | 40 ms | 1.62 s | 2.62 |

The score is the weighted geometric mean of each mode's slowdown against its best layout, lower being
better, over every mode the sweep measures. The best layout was only 4% better overall than the
cube's own, and doubling the weight of spectra or of planes moved first place to a different layout
each way: for a plane this size, 4 MiB chunks four channels deep are about as good as any.

With every mode weighted alike, 128 x 128 x 128 served the 512-square cube best, and 512 x 512 x 16
the 2048-square one, where the same 128-deep chunks made a plane take nearly three seconds. A chunk
one plane deep -- what a cube converted plane by plane gets -- gives the fastest channel changes and
the slowest spectra by a factor of ten or more. Which to prefer is a judgement about how a site's
users work; the sweep's `[weights]` makes it explicit, and its report says when doubling the weight of
planes or of spectra would change the answer. For the 512-square cube it would: weighting planes
twice as heavily brings the cube's own layout back to first place.

**A rule of thumb.** Each operation took time roughly in proportion to the pixels it decoded: a plane
decodes the plane's area times the chunk depth, a spectrum the chunk's area times the number of
channels. For chunks of about the same size the two balance at a depth of about

    depth ≈ sqrt(chunk pixels × channels / plane pixels)

which gives about 180 for the 512-square cube (128 measured best) and 30 to 45 for the 2048-square one
(16 to 32 measured best), with chunks of 2 to 4 million pixels. It is a starting point for the grid a
sweep measures, not an answer: regions, animations and cube histograms pull in other directions, and
it was fitted to two cubes.

For the 7763 x 4742 x 128 cube the rule gives about 4 for chunks of 4 million pixels and about 2 for
chunks of a million; 1024 x 1024 x 4 measured best, and the million-pixel chunks one and four deep
came close behind. By the rule a million-pixel chunk -- 4 MiB, the size CARTA's own tests have mostly
used, as 512 x 512 x 4 or 256 x 256 x 16 -- suits a 2048-square, 2048-channel cube at about 22 channels
deep, close to 256 x 256 x 16, and a 512-square, 4096-channel one at about 126, far deeper than
512 x 512 x 4. Put the chunk shapes a site already uses into its sweep's grid, so that the report
compares every candidate against them.

### Chunk size matters less than depth

A later set of measurements, [chunk size and depth](measurements/chunk-size.md), read synthetic cubes
in chunks of 1 to 16 MiB at several depths, on a Mac's SSD, a local NVMe disk and this Lustre. A cold
plane took time in proportion to the pixels it decoded -- the plane's area times the chunk's depth --
plus a small cost for each chunk; the report gives each machine's rate and cost. So depth decides what
a channel change costs; between 1 and 16 MiB at the same depth, size changed it by less than a factor
of two.

- **On local disks 1 MiB chunks are as good as 4 MiB.** Their planes were as fast on the SSD and about
  20% slower on NVMe, and their spectra and regions faster on both.
- **On Lustre 4 MiB chunks were the best with several users,** for planes and regions alike: with
  eight users a plane of the 7763 x 4742 pancake took 815 ms in 512 x 512 x 4 chunks, 1.34 s in
  256 x 256 x 4 and 1.51 s in 1024 x 1024 x 4. 16 MiB chunks gave the slowest regions.
- **1 MiB chunks in 4 MiB shards read planes on Lustre nearly as fast as 4 MiB chunks** and kept the
  1 MiB chunks' spectra; 16 MiB shards gave most of that back with eight users.

So 4 MiB remains the size to start a sweep from on Lustre, and on a local disk 1 MiB is as good.

### Animations stall where a run of chunks begins

An animation reads plane after plane through the backend's cache. Within a run of chunks along the
spectrum every frame after the first finds its chunks decoded and costs a few milliseconds; the frame
that enters the next run reads and decodes all of it, which costs what jumping to a channel does. At
CARTA's 5 frames a second, a 7763 x 4742 cube in 512 x 512 x 4 chunks stalled for up to 135 ms every
fourth frame for one user, and in 512 x 512 x 16 chunks for up to 400 ms every sixteenth; 512-square
cubes barely stalled at all.

carta-backend's Zarr support reads ahead, and so does the bench with `--animation-prefetch`:

- after the first frame of each run, it decodes the next run into the cache, one run at a time;
- once a frame is late while a prefetch is under way, it stops for the rest of the animation;
- a frame that catches a prefetch waits for it rather than decoding the same chunks again
  ([ADR 0015](adr/0015-a-cached-chunk-is-not-checked-again.md)).

For one user reading ahead hid every stall, at 5 and at 10 frames a second, in every case measured,
including the 512 x 512 x 4 chunks with the default 1 GiB cache, which holds only 1.8 runs: a run is
the plane's area times the chunk depth times four bytes, 589 MB for the 512 x 512 x 4 chunks and 2.4
GB for 512 x 512 x 16. An earlier version that read the next run's first plane, 147 MB to allocate and
fill, still stalled there at 10 frames a second.

With eight users animating the same large cube on one machine, frames late out of 496 per setting,
or 1008 for the 16-deep chunks:

| Chunks, cache | fps | Without | With prefetch |
|---|---|---|---|
| 512 x 512 x 4, 1 GiB | 5 | 111 | 45 |
| 512 x 512 x 4, 4 GiB | 5 | 112 | 35 |
| 512 x 512 x 16, 8 GiB | 5 | 48 | 59 |
| 512 x 512 x 4, 1 GiB | 10 | 141 | 159 |
| 512 x 512 x 4, 4 GiB | 10 | 136 | 155 |
| 512 x 512 x 16, 8 GiB | 10 | 203 | 343 |

At CARTA's 5 frames a second it cut the 4-deep chunks' late frames by more than half. In the 16-deep
chunks it shortened the longest stall from 1.3 s to 1.2 s but made 11 more frames late. At 10 frames
a second it made every case worse, because the one prefetch each animation made before
stopping still decoded a whole run while the machine had no time over. Stopping on the second
prefetch the animation caught up with, an earlier rule, was worse still: 541 late in the last row.
So a backend that prefetches should do it only while the cache holds two runs for each cube being
animated, and the sweep's trade-off table gives the run size of each layout.

The [chunk-size measurements](measurements/chunk-size.md) found the same over 15 layouts of the 7763 x
4742 x 256 pancake, 32 frames at a time at 5 frames a second. With a 5 GiB cache, which holds two runs
of every layout up to sixteen channels deep, no layout had more than 6% of its frames late for one
user, on the SSD, NVMe or Lustre. With the default 1 GiB, which on that plane holds two runs only of
chunks one or two channels deep, the frames late on NVMe were 22% in chunks four deep, 9% eight deep,
3% sixteen deep and 97% thirty-two deep. With four or eight users animating at once and the 5 GiB
cache, chunks one channel deep had 71 to 94% of their frames late on NVMe and 14 to 80% on Lustre, and
chunks sixteen deep 0 to 8% on either, after a first frame of several seconds.

### Through carta-backend, the rule holds for a pancake

The read-path test set (`tools/testset`) holds a synthetic pancake, 7763 x 4742 x 256, made from one
FITS file in 4 MiB chunks of two shapes: 256 x 256 x 16, which sites have used, and 512 x 512 x 4,
which the rule of thumb above points to (about 2.7 channels deep). Measured through carta-backend's
own read paths on a machine with 28 logical cores and local NVMe and a warm page cache, the shallow
chunk was the better one for planes -- a channel shown, changed or animated -- and the deep one for
spectra. A moment's peak is how far it raised carta-backend's resident memory above what the session
held before it:

| operation | 256 x 256 x 16 | 512 x 512 x 4 |
|---|---|---|
| animation, no cache (mean ms a frame) | 190 ms | 60 ms (the FITS file: 60 ms) |
| animation, 1 GiB cache (mean ms a frame) | 104 ms | 26 ms |
| eight channel changes, 1 GiB cache | 0.75 s | 0.18 s |
| cursor spectrum, no cache | 8 ms | 30 ms |
| cube histogram, exact, 1 GiB cache | 17.4 s | 12.7 s |
| moment, 1 GiB cache | 27 s, peak +446 MiB | 24 s, peak +1355 MiB |

A plane in 16-deep chunks decodes 16 channels to show one; in 4-deep chunks, 4. A spectrum pays the
other way, reading four times as many chunks, and still takes tens of milliseconds. A moment's peak
grows with the chunk's area rather than its depth: carta-backend reads it in slabs a whole column of
chunks deep, here 512 x 512 x 256 against 256 x 256 x 256. Neither layout let the backend read ahead
with the default 1 GiB cache, since two runs of chunks come to 4.6 GiB and 1.25 GiB; the shallow
chunk did not need it.

### And for a cigar

The test set's cigar, 512 x 512 x 30000, is chunked 128 x 128 x 64 as sites have shaped their chunks,
and 64 x 64 x 256, nearer the depth the rule points to for it (about 350). Measured the same way, with
a 1 GiB cache unless noted, the deep chunk was the better one for what a cigar is opened for (a moment's
peak as above, above what the session held before it):

| operation | 128 x 128 x 64 | 64 x 64 x 256 |
|---|---|---|
| cursor spectrum | 300 ms | 70 ms |
| region spectrum | 1.20 s | 0.51 s |
| moving the cursor nearby (mean of seven moves after the first) | 271 ms | 1.3 ms |
| moving the cursor across the image (mean of seven moves after the first) | 213 ms | 62 ms |
| animation, no cache (mean ms a frame) | 12.9 ms | 29.4 ms |
| animation, 1 GiB cache (mean ms a frame) | 1.4 ms | 3.0 ms |
| channel change, no cache | 7 ms | 24 ms |
| moment | 42 s, peak +4280 MiB | 41 s, peak +1146 MiB |

A spectrum in 64-deep chunks reads 469 of them, 1.8 GiB decoded, which is more than half the cache,
so carta-backend reads it past the cache and moving the cursor reads every chunk again. In 256-deep
chunks it is 118, 472 MiB, which the cache keeps, so moving the cursor within the same chunks reads
nothing. The cost is on the plane: one channel decodes 64 chunks of 256 channels, 256 MiB, to show
1 MiB. That is still under 30 ms a frame on average without a cache, 40 ms at most, which leaves most of
an animation's 200 ms a frame for the rest of a frame's work -- these times are the plane's read alone,
over the 32 frames of one chunk's depth -- and a 1 GiB cache holds the run of chunks it reads ahead.

The exact cube histogram took 27 s against 14 s. Its plane histograms woke the worker pool once per
channel of each band a read was cut into, which narrow chunks made small (#6); carta-zarr since #7
splits a read of many channels by channel instead, and those plane histograms take 5.7 s rather than
10.7 s, as long as at 128 x 128 x 64 (called directly, without a cache). Through carta-backend, with the
1 GiB cache, the whole exact cube histogram then takes 16.8 s against 13.6 s at 128 x 128 x 64.

### Use more file-reading threads on Lustre

`--zarr_file_io_threads 8` was better than the default of 2 on Lustre in every sweep: a plane of the
128-deep layout went from 146 ms to 104 ms and a spectrum from 222 ms to 157 ms; a region of the
cube's own layout from 6.35 s to 3.33 s; a plane of the 2048-square cube's 512 x 512 x 16 layout from
478 ms to 358 ms; and in the 7763 x 4742 cube's own 512 x 512 x 4 chunks a plane from 995 ms to
838 ms and a spectrum from 240 ms to 167 ms. On BeeGFS the difference was within the noise.
`--zarr_data_copy_threads` at half the logical cores was as good as all of them.

### Striping and shards

Striping a dataset's directory across one OST or across eight made no difference beyond 5% for chunks
of a few megabytes: a chunk file smaller than the stripe size lives on one OST either way. What
striping cannot change is the number of files, which every open and every stat sends to the metadata
server that all clients share: a 32 GiB cube in 128-deep chunks is 4,000 files, and a terabyte cube in
the same chunks would be over 100,000. Shards -- chunks gathered into one file, each still read on its
own -- are the remedy, and are worth including in a sweep on a busy file system. The sweep accepts
them relative to the chunk, `frequency*8` or `l*4,m*4`. The sweep ranks only first touches, operations
that read no chunk another operation of the same trial read before or alongside them. A sharded
layout's first touches are counted by chunk, as an unsharded one's are, so a read that finds its
shard's index already in the page cache still counts; if the file system reads ahead within a shard's
file, a sharded layout's planes and spectra may look faster than a first touch would be. The
[chunk-size measurements](measurements/chunk-size.md#the-pancake-in-shards) compare 1 MiB chunks in 4
and 16 MiB shards with the unsharded layouts.

## Measuring your own system

The answer on another file system, another cube size or another number of users can differ, and the
tools exist so that it can be measured rather than assumed.

1. **Build** carta-zarr in Release with `-DCARTA_ZARR_BUILD_BENCH=ON`; a Debug build is several times
   slower and the report says so.
2. **Write a configuration** from [`example-sweep.toml`](../tools/zarr-bench/example-sweep.toml):
   the bench, a work directory on the file system to measure, a real cube to rewrite (or a synthetic
   shape), the number of users at the busiest time and on an ordinary day, the layouts to try, and the
   reader settings to try on the best of them. If the site holds cubes of very different sizes, list
   each as a `[[source.shape]]`.
3. **Check it:** `./sweep.py my-sweep.toml --dry-run` says what would run, how much space it needs,
   how caches will be emptied, and whether striping can be set, and writes nothing.
4. **Run it:** `./sweep.py my-sweep.toml`. It writes, measures and deletes one layout at a time, and
   the same command resumes it after an interruption. A sweep of four layouts of a 32 GiB cube took
   about two hours on that server's Lustre.
5. **Read `summary.md`**: the layout and the flags to use, as they would be pasted into the backend's
   command line or `settings.json`, then every warning, then the tables behind them.

### Empty the right caches

A read that finds its data in memory measures memory. Between a CARTA server and the disks there are
at least three caches: carta-zarr's own (each read in the bench starts without it), the client's page
cache, and the storage servers'. On that server's Lustre a whole-cube read took 1.24 s with everything
cached, 2.0 s with the client's page cache emptied, and 5.5 s with the servers' emptied too -- but
11.9 s the first time after the cube was written. Something below what an ordinary user can empty
keeps data that was read recently, and no later read was as slow as the first.

The sweep therefore runs the bench once, untimed, every time it writes a dataset, before measuring
it. Every layout is then measured in the same state, which is what a comparison needs, and its times
are those of data read recently rather than data never read: there, a third to half faster for a
region of the cube's own layout. A first view of a cube nobody has opened in a while will be slower
than the report says.

- With root, `drop_caches` empties the client's page cache; give the sweep a `drop_cache_cmd` that
  empties the servers' too if the storage administrators provide one.
- Without root on Lustre,
  [`drop-lustre-cache.py`](../tools/zarr-bench/drop-lustre-cache.py) empties both: the client's by
  `posix_fadvise`, and the servers' by `lfs ladvise -a dontneed`, which Lustre 2.9 and later accept
  from a file's owner.
- Without root on BeeGFS nothing reaches the storage servers' caches, and the report says reads were
  cold on this host only. Layouts measured second can then look faster than they are.

### Things that mislead a measurement

- **Too shallow a crop.** A plane or a spectrum is timed as the first read of its chunks. A crop of
  1,024 channels holds only eight chunks 128 deep, and sixteen users each reading sixteen planes
  share them, so most reads find their chunks already in the page cache. The report counts the reads
  that did not, ranks on those, and warns when there are too few, with what to change.
- **The order of runs on Lustre.** With the caches emptied the same way before each, the first run
  on a dataset after it was written, or after other datasets were read, was slower than the ones that
  followed it: a region of a cube in one-channel chunks took 1.2 s to 1.7 s the first time and 0.8 s
  every time after. Opening every file first, reading a byte of each, reading all of them, or waiting
  two minutes did not change it; only a run of the bench did. A comparison of two builds once showed
  a 50% to 70% gain that vanished when the order was alternated. The sweep's untimed run puts every
  measured run after one; compare builds or settings by hand the same way.
- **One user is not eight.** Settings that help one user can do nothing for eight, whose processes
  already keep every core busy; the report gives the target number of users, and how the choice does
  with one and with the typical number.
- **Data that fits in memory.** A crop smaller than RAM measures caching unless caches are emptied
  before every trial; the sweep refuses one when they cannot be. Validation re-measures the
  recommendation on a copy larger than RAM, which is impossible on a server with as much memory as
  that one's.
- **Synthetic pixels.** A layout's cost depends on how well its pixels compress, and noise compresses
  worse than sky. Rewrite a real cube when there is one.

## Very large planes

Cubes from the SKA or the ngVLA may be tens of thousands of pixels on a side and hundreds of channels
deep. Two things change.

**Every channel change reads the whole plane.** carta-backend reads the plane of a Zarr cube in full
and downsamples it for display, so a 30,000 x 30,000 plane is 3.6 GB of float32 per channel change,
however shallow the chunks. No layout makes that interactive; what would is a downsampled copy beside
the full one, which carta-backend does not read yet. Until it does, the rule of thumb above points to
chunks one channel deep for such a plane, and a spectrum then reads a chunk per channel.

**A layout cannot be written at full depth for every try.** A sweep shape may keep the full plane and
fewer channels, with `channels` saying how deep the cubes it stands for are:

```toml
[[source.shape]]
name = "ska"
synthetic = "frequency=64,polarization=1,l=32768,m=32768"
channels = 512
```

Planes and animations are then measured at their full size, and the summary scales the modes that
read every channel to the full depth, as an estimate. A chunk deeper than the shape's channels is
skipped, since it would not be the chunk a full cube has.

## Where the numbers came from

The sweeps behind this page ran on 2026-10-01 and 2026-10-02 with carta-zarr at the commits that
added each tool, eight users, two or three trials, and every cache emptied as described above. Their
configurations and full reports were kept on the machines they ran on. The measurements since keep
their medians in `docs/measurements/`, beside a report of the method and every number: the
chunk-size measurements ran on 2026-10-14 and 2026-10-15, and
[their report](measurements/chunk-size.md#reproducing) says how to run them again. Either way they
describe those machines; measure yours rather than reuse them.

The comparison through carta-backend ran on 2026-10-08 and 2026-10-09, on the read-path test set built
by `tools/testset/build.sh`, with carta-zarr at 21f18a2 and carta-backend's `test_zarr` branch at
6e406c80 (its `MeasureReadPaths` test, one operation per process), on a single machine with local NVMe
rather than a parallel file system. The pancake's 256 x 256 x 16 figures are the earlier runs of the
same test at `test_zarr` 18d82c15, and its cube histogram and moment at 8a118389. Those differ from
6e406c80 only in reading a spectrum around the session cache when its chunks outgrow half of it,
which the pancake's spectrum, 16 chunks or 64 MiB, does not. The cigar's plane histograms before and
after #7 were timed by calling carta-zarr directly, at 21f18a2 and at #7's first version (59789f8), on
the same machine, and the exact cube histogram again through carta-backend on 2026-10-10 with that
version swapped in, the mean of two runs. #7's final version (1c2ac7c) splits by channel only where its
busiest worker bins no more rows that way; on 28 threads it timed the same as the first version, 5.7 s
and 16.3 s without a cache.
