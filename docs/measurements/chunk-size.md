# Chunk size and depth, measured on SSD, NVMe and Lustre

How the size of a chunk and its depth in channels decide what a CARTA user waits for, on three kinds
of storage, measured with `carta-zarr-bench` from 2026-10-14 to 2026-10-15
([#17](https://github.com/CARTAvis/carta-zarr/issues/17)). [Storage tuning](../storage-tuning.md)
keeps the advice drawn from it; this report keeps the method and every number, and
[`chunk-size/`](chunk-size) the medians and the script that makes them.

## What the measurements say

- **Depth decides how long a channel takes to show; size matters much less.** 512 × 512 chunks
  sixteen channels deep made a plane of the 7763 × 4742 pancake 6.5 (Lustre) to 14 (Mac) times slower
  than chunks one deep; going from 1 to 16 MiB at the same depth changed it by less than a factor of
  two either way. A [model of a cold plane](#a-model-of-a-cold-plane) fits every machine.

- **What a plane gains from shallow chunks, a spectrum and a region lose.** The best spectrum is the
  chunk smallest across the sky and deepest in channels; the best region, one small across the sky.
  No layout is best at everything: on acdc the pancake's fastest plane, 1024 × 1024 × 1 (76 ms), took
  495 ms for a spectrum, and its fastest spectrum, 128 × 128 × 16 (15 ms), 1.38 s for a plane. The
  rule of thumb in [storage tuning](../storage-tuning.md), depth ≈ √(chunk pixels × channels ÷ plane
  pixels), points at the right end on all three cubes: 1 to 3 for the pancake in 1 to 4 MiB chunks,
  90 to 180 for the 512 × 512 × 8192 cigar, below 1 for the 16384 × 16384 × 32 plane.

- **On local storage 1 MiB chunks are as good as 4 MiB; on Lustre 4 MiB are better.** On the Mac's SSD
  the cost of a chunk is too small to see: 1 MiB chunks read planes within 5% of 4 MiB ones at the same
  depth, and spectra and regions faster. On acdc's NVMe their planes were about 20% slower, their
  spectra and regions faster. On Lustre a chunk costs more than twice what it does on NVMe, and with
  several users the 4 MiB layouts were the best for planes and regions: with eight, 512 × 512 × 4 read
  a plane in 815 ms against 1.34 s for 256 × 256 × 4 (1 MiB) and 1.51 s for 1024 × 1024 × 4
  (16 MiB), and 256 × 256 × 16 a region in 4.22 s against 6.12 s for 128 × 128 × 16 and 8.54 s for
  512 × 512 × 16. With several users 16 MiB chunks were mostly the slowest for regions on acdc too.

- **On Lustre, shards of 4 MiB give 1 MiB chunks the planes of 4 MiB ones.** 256 × 256 × 4 chunks in
  512 × 512 × 4 shards read a plane for eight users in 646 ms against 1.34 s unsharded, near
  512 × 512 × 4's 815 ms, and kept the 1 MiB chunks' spectra. Shards of 16 MiB lost much of that gain
  with eight users. On acdc's NVMe shards gained a little on planes and lost some regions.

- **An animation keeps up for one user if its cache holds two runs of chunks.** With a 5 GiB cache,
  which holds two runs of every pancake layout up to sixteen channels deep, no layout of the pancake
  had more than 6% of its frames late for one user, on any machine. carta-backend's default 1 GiB
  holds two runs only for chunks one or two channels deep on a plane this size -- a run is the plane's
  area times the depth times four bytes -- so in deeper chunks it reads nothing ahead and stalls where
  each run begins. On acdc that made 22% of frames late in chunks four deep, 9% eight deep, 3% sixteen
  deep and 97% thirty-two deep, and Lustre was much the same; the Mac decoded fast enough to stay on
  time in chunks up to eight deep.

- **With several users at once, shallow chunks make animations late.** In chunks one channel deep each
  frame costs every user a cold plane. With a 5 GiB cache, on acdc 71 to 94% of frames were late in
  chunks one channel deep with four or eight users, and 2 to 3% in chunks sixteen deep, which decode a
  run once and serve sixteen frames from the cache -- after a first frame of 4 to 5 seconds. On Lustre
  14 to 80% were late one deep and 0 to 8% sixteen deep.

- **A 16384 × 16384 plane cannot animate at 5 frames a second on these machines.** A cold plane took
  0.19 s at best on the Mac and 0.52 s on acdc, and deeper chunks multiplied that by roughly their
  depth. One user had 10% of frames late at best on the Mac and 22% on acdc, and several users 89 to
  97%. Its regions took about the same time in every layout.

- **The synthetic cube reads like a real one,** within 3 to 4% for cold planes on acdc
  ([check](#synthetic-and-real-pixels)), and **`purge` empties the Mac's caches**
  ([check](#did-purge-empty-the-macs-caches)).

## A model of a cold plane

A cold plane read took the pixels it decoded -- the plane's area times the chunk's depth -- at a rate
set by the machine, plus a small cost for each chunk. Fitted to every 1-user plane, the model is off
by 5 to 10% at the median:

| machine and storage | cold rate | cost per chunk | median error of the fit |
|---|---|---|---|
| Mac mini, Apple SSD | 2.5 Gpx/s | 0.016 ms | 10%, 37 layouts |
| acdc, local NVMe | 0.59 Gpx/s | 0.13 ms | 5%, 37 layouts |
| almat3, Lustre | 0.90 Gpx/s | 0.35 ms | 10%, 15 layouts |

Gpx/s is a billion float32 pixels decoded a second, for one user. Lustre's fit has only the pancake.
The error is largest on the cigar's planes, which take a few milliseconds.

## How they were measured

### Cubes and layouts

Every cube is synthetic (`tools/zarr-bench/generate.py --synthetic --seed 1`), float32, compressed
with Blosc zstd 5 with byte shuffle, and the same on every machine:

| cube | l × m × channels | measured on | layouts |
|---|---|---|---|
| pancake | 7763 × 4742 × 256 | Mac, acdc, Lustre | 1, 2, 4, 8 and 16 MiB, three depths each |
| pancake in shards | 7763 × 4742 × 256 | acdc, Lustre | 1 MiB chunks in 4 and 16 MiB shards |
| cigar | 512 × 512 × 8192 | Mac, acdc | 1, 4 and 16 MiB, depths 1 to 1024 |
| wide | 16384 × 16384 × 32 | Mac, acdc | 1, 4 and 16 MiB, depths 1, 4 and 16 |

Chunks are square across the sky with sides a power of two, so that their edges fall on CARTA's tile
edges at every zoom. The wide cube's spectra are 32 channels long, so they do not stand for a
spectrum of a deep cube. The layouts, as l × m × channels of a chunk:

| cube | 1 MiB | 2 MiB | 4 MiB | 8 MiB | 16 MiB |
|---|---|---|---|---|---|
| pancake | 512²×1, 256²×4, 128²×16 | 512²×2, 256²×8, 128²×32 | 1024²×1, 512²×4, 256²×16 | 1024²×2, 512²×8, 256²×32 | 2048²×1, 1024²×4, 512²×16 |
| cigar | 512²×1, 256²×4, 128²×16, 64²×64, 32²×256 | | 512²×4, 256²×16, 128²×64, 64²×256 | | 512²×16, 256²×64, 128²×256, 64²×1024 |
| wide | 512²×1, 256²×4, 128²×16 | | 1024²×1, 512²×4, 256²×16 | | 2048²×1, 1024²×4, 512²×16 |

The pancake in shards keeps its 1 MiB chunks 512²×1, 256²×4 and 128²×16, and gathers each into
shards of four chunks (1024²×1, 512²×4, 256²×16) and of sixteen (2048²×1, 1024²×4, 512²×16).

### What was timed

For each layout, a script wrote the cube, read one plane of it untimed
(a cube just written reads slower the first time on Lustre), measured it at each user count, and
deleted it. Each user is a process of its own with its own cache, and each trial empties every cache
first, then starts fresh processes. Three trials. In each trial every user does:

- **plane**: 16 single-channel planes at random channels, as a cold channel switch reads them;
- **spectrum**: 32 cursor spectra through every channel, at random pixels;
- **region**: one region of 5% of the plane, reduced over every channel;
- **animation**: two runs of 32 frames at 5 frames a second, reading ahead (`--animation-prefetch`).

Every mode is measured with a 5 GiB cache, which holds two runs of chunks -- the plane's area times
the depth, twice -- of every pancake layout up to sixteen channels deep, and the animation again with
carta-backend's default 1 GiB. A frame is late when it is not ready by the time it is due.

A plane, spectrum or region counts only if it was a *first touch*: it read no chunk that another
operation of the same trial read before or alongside it, as `tools/zarr-bench/sweep.py` ranks them;
otherwise its time could be the page cache's. Where fewer than three operations a cell, one a trial,
were first touches -- deep chunks with several users, every plane of the wide cube with several --
the tables show "–".

### Machines

| | Mac mini | acdc | almat3 |
|---|---|---|---|
| CPU | Apple M4 Pro, 12 cores | Intel i7-14700F, 28 logical cores | 2 × Xeon Gold 6134, 32 logical cores |
| memory | 64 GiB | 125 GiB | 1 TiB |
| storage | internal SSD (APFS) | local NVMe (btrfs) | Lustre 2.15, 8 stripes of 4 MiB |
| caches emptied by | `purge` | `posix_fadvise(DONTNEED)` | client page cache and the servers' caches (`lfs ladvise`) |
| I/O, decode threads a user | 4, 12 | 4, 28 | 8, 16 |
| users | 1 | 1, 4, 8 | 1, 4, 8 |

`carta-zarr-bench` was built in Release from d846091 on the Mac and acdc, and from 81bfbed on almat3,
which differs from it only in the README. [Reproducing](#reproducing) gives the commands.

## The pancake, 7763 × 4742 × 256

Each cell is the median time of one operation over three trials, first touches only; "–" marks a
cell with fewer than three. Columns name the mode and the machine. Planes slow with depth on every
machine while spectra speed up, and on Lustre with several users the 4 MiB layouts read planes and
regions fastest.

### One user

| chunk | MiB | plane, mac | spectrum, mac | region, mac | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 18 ms | 42 ms | 616 ms | 89 ms | 201 ms | 2.25 s | 129 ms | 258 ms | 3.50 s |
| 256 × 256 × 4 | 1 | 67 ms | 11 ms | 556 ms | 342 ms | 49 ms | 1.94 s | 419 ms | 68 ms | 3.00 s |
| 128 × 128 × 16 | 1 | 260 ms | 3.8 ms | 432 ms | 1.38 s | 15 ms | 1.71 s | 1.49 s | 26 ms | 2.62 s |
| 512 × 512 × 2 | 2 | 34 ms | 43 ms | 620 ms | 159 ms | 179 ms | 2.01 s | 176 ms | 212 ms | 3.45 s |
| 256 × 256 × 8 | 2 | 126 ms | 11 ms | 548 ms | 594 ms | 47 ms | 1.70 s | 637 ms | 60 ms | 3.07 s |
| 128 × 128 × 32 | 2 | 507 ms | 4.1 ms | 448 ms | 2.39 s | 16 ms | 1.53 s | 2.11 s | 24 ms | 2.18 s |
| 1024 × 1024 × 1 | 4 | 20 ms | 136 ms | 969 ms | 76 ms | 495 ms | 3.13 s | 98 ms | 631 ms | 3.52 s |
| 512 × 512 × 4 | 4 | 66 ms | 45 ms | 589 ms | 285 ms | 165 ms | 1.81 s | 240 ms | 166 ms | 2.23 s |
| 256 × 256 × 16 | 4 | 248 ms | 12 ms | 522 ms | 1.11 s | 44 ms | 1.56 s | 861 ms | 47 ms | 2.22 s |
| 1024 × 1024 × 2 | 8 | 38 ms | 135 ms | 952 ms | 150 ms | 494 ms | 3.02 s | 145 ms | 484 ms | 3.02 s |
| 512 × 512 × 8 | 8 | 128 ms | 47 ms | 574 ms | 564 ms | 162 ms | 1.80 s | 422 ms | 139 ms | 1.84 s |
| 256 × 256 × 32 | 8 | 486 ms | 15 ms | 519 ms | 2.20 s | 48 ms | 1.61 s | 1.58 s | 50 ms | 1.81 s |
| 2048 × 2048 × 1 | 16 | 33 ms | 469 ms | 1.31 s | 82 ms | 1.53 s | 4.09 s | 124 ms | 1.45 s | 3.92 s |
| 1024 × 1024 × 4 | 16 | 75 ms | 135 ms | 945 ms | 282 ms | 464 ms | 2.96 s | 242 ms | 420 ms | 2.68 s |
| 512 × 512 × 16 | 16 | 258 ms | 47 ms | 561 ms | 1.08 s | 164 ms | 1.81 s | 834 ms | 147 ms | 1.66 s |

### Four users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 320 ms | 586 ms | 5.71 s | 191 ms | 414 ms | 4.32 s |
| 256 × 256 × 4 | 1 | 1.21 s | 170 ms | 5.53 s | 744 ms | 111 ms | 4.13 s |
| 128 × 128 × 16 | 1 | 3.51 s | 45 ms | 5.21 s | 2.47 s | 31 ms | 3.91 s |
| 512 × 512 × 2 | 2 | 519 ms | 526 ms | 5.14 s | 294 ms | 335 ms | 3.76 s |
| 256 × 256 × 8 | 2 | 2.10 s | 155 ms | 4.98 s | 1.27 s | 101 ms | 4.29 s |
| 128 × 128 × 32 | 2 | – | 40 ms | 4.63 s | – | 28 ms | 3.28 s |
| 1024 × 1024 × 1 | 4 | 243 ms | 1.26 s | 7.33 s | 124 ms | 928 ms | 3.74 s |
| 512 × 512 × 4 | 4 | 1.01 s | 470 ms | 4.81 s | 444 ms | 290 ms | 2.93 s |
| 256 × 256 × 16 | 4 | 3.03 s | 140 ms | 4.67 s | 1.52 s | 81 ms | 3.21 s |
| 1024 × 1024 × 2 | 8 | 464 ms | 876 ms | 6.85 s | 372 ms | 699 ms | 5.01 s |
| 512 × 512 × 8 | 8 | 1.91 s | 469 ms | 4.93 s | 1.47 s | 383 ms | 4.13 s |
| 256 × 256 × 32 | 8 | – | 127 ms | 4.53 s | – | 109 ms | 4.29 s |
| 2048 × 2048 × 1 | 16 | 168 ms | – | – | 170 ms | – | – |
| 1024 × 1024 × 4 | 16 | 964 ms | 1.03 s | 6.64 s | 883 ms | 954 ms | 6.18 s |
| 512 × 512 × 16 | 16 | 2.85 s | 406 ms | 4.30 s | 2.57 s | 425 ms | 4.75 s |

### Eight users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 603 ms | 1.06 s | 9.46 s | 370 ms | 684 ms | 7.30 s |
| 256 × 256 × 4 | 1 | 2.34 s | 335 ms | 8.82 s | 1.34 s | 206 ms | 6.17 s |
| 128 × 128 × 16 | 1 | – | 88 ms | 9.07 s | – | 47 ms | 6.12 s |
| 512 × 512 × 2 | 2 | 1.07 s | 945 ms | 8.29 s | 566 ms | 554 ms | 6.14 s |
| 256 × 256 × 8 | 2 | 3.34 s | 302 ms | 7.75 s | 1.72 s | 191 ms | 5.73 s |
| 128 × 128 × 32 | 2 | – | 77 ms | 8.59 s | – | 38 ms | 5.07 s |
| 1024 × 1024 × 1 | 4 | 504 ms | – | – | 210 ms | – | – |
| 512 × 512 × 4 | 4 | 2.00 s | 917 ms | 8.27 s | 815 ms | 450 ms | 4.39 s |
| 256 × 256 × 16 | 4 | – | 271 ms | 6.79 s | – | 134 ms | 4.22 s |
| 1024 × 1024 × 2 | 8 | 956 ms | – | – | 630 ms | – | – |
| 512 × 512 × 8 | 8 | 2.50 s | 679 ms | 8.07 s | 1.89 s | 535 ms | 7.66 s |
| 256 × 256 × 32 | 8 | – | 248 ms | 6.96 s | – | 214 ms | 7.13 s |
| 2048 × 2048 × 1 | 16 | 430 ms | – | – | 413 ms | – | – |
| 1024 × 1024 × 4 | 16 | 1.81 s | – | – | 1.51 s | – | – |
| 512 × 512 × 16 | 16 | – | 728 ms | 7.91 s | – | 669 ms | 8.54 s |

### Animation

Each cell is the share of frames late with a cache that holds two runs (5 GiB), the same with the
default 1 GiB, and the median first frame with 5 GiB.

| chunk | MiB | mac, 1 user | acdc, 1 user | acdc, 4 users | acdc, 8 users | lustre, 1 user | lustre, 4 users | lustre, 8 users |
|---|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 0% / 0%, 25 ms | 0% / 1%, 89 ms | 94% / 94%, 355 ms | 93% / 93%, 326 ms | 0% / 0%, 128 ms | 49% / 49%, 250 ms | 60% / 54%, 229 ms |
| 256 × 256 × 4 | 1 | 0% / 0%, 74 ms | 0% / 22%, 351 ms | 21% / 22%, 1.38 s | 21% / 22%, 1.18 s | 0% / 22%, 411 ms | 19% / 22%, 924 ms | 18% / 22%, 899 ms |
| 128 × 128 × 16 | 1 | 0% / 3%, 284 ms | 0% / 3%, 1.41 s | 3% / 91%, 5.50 s | 3% / 96%, 5.02 s | 0% / 3%, 1.40 s | 3% / 51%, 3.66 s | 8% / 88%, 3.58 s |
| 512 × 512 × 2 | 2 | 0% / 0%, 40 ms | 0% / 0%, 162 ms | 45% / 45%, 615 ms | 45% / 45%, 581 ms | 0% / 0%, 180 ms | 6% / 8%, 315 ms | 29% / 23%, 347 ms |
| 256 × 256 × 8 | 2 | 0% / 0%, 129 ms | 0% / 9%, 602 ms | 9% / 9%, 2.43 s | 9% / 9%, 2.36 s | 0% / 9%, 630 ms | 4% / 9%, 1.60 s | 6% / 9%, 1.58 s |
| 128 × 128 × 32 | 2 | 0% / 95%, 507 ms | 0% / 97%, 2.60 s | 0% / 97%, 9.77 s | 0% / 97%, 9.06 s | 0% / 97%, 2.09 s | 0% / 97%, 5.40 s | 1% / 97%, 6.68 s |
| 1024 × 1024 × 1 | 4 | 0% / 0%, 30 ms | 0% / 0%, 78 ms | 94% / 94%, 292 ms | 91% / 90%, 282 ms | 2% / 2%, 90 ms | 14% / 28%, 172 ms | 63% / 35%, 244 ms |
| 512 × 512 × 4 | 4 | 0% / 0%, 67 ms | 0% / 22%, 317 ms | 20% / 22%, 1.17 s | 21% / 22%, 1.19 s | 0% / 12%, 240 ms | 0% / 22%, 604 ms | 12% / 22%, 562 ms |
| 256 × 256 × 16 | 4 | 0% / 3%, 243 ms | 0% / 3%, 1.17 s | 3% / 87%, 4.68 s | 3% / 94%, 4.68 s | 0% / 3%, 822 ms | 0% / 64%, 2.39 s | 5% / 92%, 2.22 s |
| 1024 × 1024 × 2 | 8 | 0% / 0%, 46 ms | 0% / 0%, 145 ms | 44% / 44%, 548 ms | 42% / 43%, 551 ms | 6% / 4%, 162 ms | 45% / 45%, 509 ms | 48% / 33%, 507 ms |
| 512 × 512 × 8 | 8 | 0% / 0%, 128 ms | 0% / 9%, 565 ms | 9% / 9%, 2.24 s | 7% / 9%, 2.17 s | 0% / 9%, 430 ms | 8% / 9%, 1.88 s | 8% / 9%, 1.94 s |
| 256 × 256 × 32 | 8 | 0% / 97%, 483 ms | 0% / 97%, 2.23 s | 0% / 97%, 8.92 s | 0% / 97%, 8.69 s | 0% / 96%, 1.63 s | 0% / 97%, 8.35 s | 0% / 97%, 8.33 s |
| 2048 × 2048 × 1 | 16 | 0% / 0%, 40 ms | 0% / 0%, 86 ms | 71% / 72%, 267 ms | 86% / 88%, 289 ms | 4% / 4%, 140 ms | 69% / 68%, 296 ms | 80% / 74%, 339 ms |
| 1024 × 1024 × 4 | 16 | 0% / 0%, 87 ms | 0% / 22%, 280 ms | 20% / 22%, 1.10 s | 20% / 22%, 1.12 s | 0% / 45%, 231 ms | 40% / 49%, 1.09 s | 42% / 47%, 1.10 s |
| 512 × 512 × 16 | 16 | 0% / 3%, 258 ms | 0% / 3%, 1.07 s | 3% / 84%, 4.52 s | 2% / 96%, 4.26 s | 0% / 3%, 861 ms | 3% / 79%, 4.29 s | 3% / 95%, 4.38 s |

## The pancake in shards

1 MiB chunks, each gathered into shards of 4 MiB (the first row of each pair) and of 16 MiB. On
Lustre, 4 MiB shards bring the planes of 1 MiB chunks near those of 4 MiB chunks; 16 MiB shards lose
much of that with eight users.

### One user

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 in 1024 × 1024 × 1 | 1 | 76 ms | 205 ms | 2.76 s | 94 ms | 326 ms | 3.55 s |
| 512 × 512 × 1 in 2048 × 2048 × 1 | 1 | 79 ms | 198 ms | 2.68 s | 85 ms | 249 ms | 3.61 s |
| 256 × 256 × 4 in 512 × 512 × 4 | 1 | 285 ms | 52 ms | 1.66 s | 234 ms | 82 ms | 2.91 s |
| 256 × 256 × 4 in 1024 × 1024 × 4 | 1 | 281 ms | 52 ms | 2.20 s | 227 ms | 58 ms | 3.82 s |
| 128 × 128 × 16 in 256 × 256 × 16 | 1 | 1.10 s | 17 ms | 1.58 s | 816 ms | 32 ms | 2.14 s |
| 128 × 128 × 16 in 512 × 512 × 16 | 1 | 1.09 s | 15 ms | 1.94 s | 850 ms | 20 ms | 2.23 s |

### Four users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 in 1024 × 1024 × 1 | 1 | 198 ms | 574 ms | 6.73 s | 130 ms | 443 ms | 4.06 s |
| 512 × 512 × 1 in 2048 × 2048 × 1 | 1 | 256 ms | 641 ms | 6.20 s | 169 ms | 359 ms | 6.30 s |
| 256 × 256 × 4 in 512 × 512 × 4 | 1 | 1.01 s | 200 ms | 5.21 s | 427 ms | 145 ms | 3.72 s |
| 256 × 256 × 4 in 1024 × 1024 × 4 | 1 | 983 ms | 167 ms | 5.72 s | 857 ms | 94 ms | 4.46 s |
| 128 × 128 × 16 in 256 × 256 × 16 | 1 | 2.83 s | 49 ms | 4.60 s | 1.36 s | 37 ms | 2.92 s |
| 128 × 128 × 16 in 512 × 512 × 16 | 1 | 2.95 s | 44 ms | 4.71 s | 2.54 s | 26 ms | 4.51 s |

### Eight users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc | plane, lustre | spectrum, lustre | region, lustre |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 in 1024 × 1024 × 1 | 1 | 478 ms | 705 ms | 12.1 s | 174 ms | 727 ms | 5.80 s |
| 512 × 512 × 1 in 2048 × 2048 × 1 | 1 | 513 ms | 886 ms | 10.8 s | 396 ms | 563 ms | 11.7 s |
| 256 × 256 × 4 in 512 × 512 × 4 | 1 | 2.07 s | 356 ms | 9.33 s | 646 ms | 254 ms | 4.89 s |
| 256 × 256 × 4 in 1024 × 1024 × 4 | 1 | 1.96 s | 359 ms | 9.61 s | 1.70 s | 163 ms | 8.31 s |
| 128 × 128 × 16 in 256 × 256 × 16 | 1 | – | 102 ms | 8.48 s | – | 62 ms | 4.13 s |
| 128 × 128 × 16 in 512 × 512 × 16 | 1 | – | 84 ms | 8.61 s | – | 43 ms | 7.83 s |

### Animation

Cells as in the pancake's animation table: late with 5 GiB / late with 1 GiB, first frame.

| chunk | MiB | acdc, 1 user | acdc, 4 users | acdc, 8 users | lustre, 1 user | lustre, 4 users | lustre, 8 users |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 in 1024 × 1024 × 1 | 1 | 1% / 0%, 77 ms | 91% / 91%, 293 ms | 91% / 90%, 282 ms | 0% / 0%, 93 ms | 2% / 1%, 181 ms | 42% / 24%, 230 ms |
| 512 × 512 × 1 in 2048 × 2048 × 1 | 1 | 0% / 0%, 83 ms | 92% / 91%, 287 ms | 93% / 91%, 300 ms | 0% / 0%, 87 ms | 59% / 62%, 219 ms | 66% / 65%, 249 ms |
| 256 × 256 × 4 in 512 × 512 × 4 | 1 | 0% / 22%, 284 ms | 21% / 22%, 1.08 s | 20% / 22%, 1.14 s | 0% / 17%, 248 ms | 1% / 22%, 618 ms | 11% / 21%, 591 ms |
| 256 × 256 × 4 in 1024 × 1024 × 4 | 1 | 0% / 22%, 296 ms | 19% / 22%, 1.14 s | 21% / 22%, 1.15 s | 0% / 17%, 223 ms | 19% / 21%, 1.08 s | 17% / 22%, 1.14 s |
| 128 × 128 × 16 in 256 × 256 × 16 | 1 | 0% / 3%, 1.16 s | 3% / 88%, 4.46 s | 3% / 96%, 4.76 s | 0% / 3%, 813 ms | 0% / 57%, 2.19 s | 6% / 89%, 2.33 s |
| 128 × 128 × 16 in 512 × 512 × 16 | 1 | 0% / 3%, 1.10 s | 3% / 81%, 4.68 s | 3% / 96%, 4.11 s | 0% / 3%, 809 ms | 3% / 57%, 4.36 s | 5% / 90%, 4.40 s |

## The cigar, 512 × 512 × 8192

No layout had a late frame in its animation, for any number of users, so it is not tabulated; its
first frames are in `medians.csv`. A chunk as wide as the plane leaves spectra and regions with
several users no first touches.

### One user

| chunk | MiB | plane, mac | spectrum, mac | region, mac | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 1.6 ms | 1.16 s | 1.25 s | 5.6 ms | 4.32 s | 4.39 s |
| 256 × 256 × 4 | 1 | 2.0 ms | 289 ms | 327 ms | 6.2 ms | 1.18 s | 1.36 s |
| 128 × 128 × 16 | 1 | 3.3 ms | 73 ms | 108 ms | 9.8 ms | 309 ms | 375 ms |
| 64 × 64 × 64 | 1 | 10 ms | 21 ms | 117 ms | 33 ms | 85 ms | 321 ms |
| 32 × 32 × 256 | 1 | 32 ms | 5.9 ms | 115 ms | 141 ms | 26 ms | 382 ms |
| 512 × 512 × 4 | 4 | 3.7 ms | 1.13 s | 1.16 s | 12 ms | 3.57 s | 3.60 s |
| 256 × 256 × 16 | 4 | 4.8 ms | 285 ms | 322 ms | 16 ms | 938 ms | 1.01 s |
| 128 × 128 × 64 | 4 | 10 ms | 73 ms | 105 ms | 31 ms | 233 ms | 313 ms |
| 64 × 64 × 256 | 4 | 31 ms | 24 ms | 113 ms | 114 ms | 76 ms | 327 ms |
| 512 × 512 × 16 | 16 | 12 ms | 1.09 s | 1.11 s | 14 ms | 3.40 s | 3.46 s |
| 256 × 256 × 64 | 16 | 17 ms | 278 ms | 307 ms | 38 ms | 944 ms | 1.00 s |
| 128 × 128 × 256 | 16 | 38 ms | 75 ms | 117 ms | 115 ms | 219 ms | 356 ms |
| 64 × 64 × 1024 | 16 | 120 ms | 30 ms | 121 ms | 434 ms | 83 ms | 388 ms |

### Four users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 2.7 ms | – | – |
| 256 × 256 × 4 | 1 | 7.9 ms | – | 3.83 s |
| 128 × 128 × 16 | 1 | 37 ms | – | 1.32 s |
| 64 × 64 × 64 | 1 | 113 ms | 200 ms | 1.01 s |
| 32 × 32 × 256 | 1 | 549 ms | 74 ms | 1.15 s |
| 512 × 512 × 4 | 4 | 8.6 ms | – | – |
| 256 × 256 × 16 | 4 | 28 ms | – | 2.92 s |
| 128 × 128 × 64 | 4 | 98 ms | – | 928 ms |
| 64 × 64 × 256 | 4 | 453 ms | 200 ms | 1.19 s |
| 512 × 512 × 16 | 16 | 25 ms | – | – |
| 256 × 256 × 64 | 16 | 96 ms | – | 2.93 s |
| 128 × 128 × 256 | 16 | 326 ms | – | 958 ms |
| 64 × 64 × 1024 | 16 | – | 165 ms | 1.09 s |

### Eight users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 3.1 ms | – | – |
| 256 × 256 × 4 | 1 | 16 ms | – | 5.04 s |
| 128 × 128 × 16 | 1 | 73 ms | – | 2.38 s |
| 64 × 64 × 64 | 1 | 211 ms | 295 ms | 1.88 s |
| 32 × 32 × 256 | 1 | 1.07 s | 142 ms | 2.20 s |
| 512 × 512 × 4 | 4 | 12 ms | – | – |
| 256 × 256 × 16 | 4 | 56 ms | – | 3.98 s |
| 128 × 128 × 64 | 4 | 193 ms | – | 1.73 s |
| 64 × 64 × 256 | 4 | 694 ms | 272 ms | 2.03 s |
| 512 × 512 × 16 | 16 | 48 ms | – | – |
| 256 × 256 × 64 | 16 | 189 ms | – | 4.01 s |
| 128 × 128 × 256 | 16 | 648 ms | – | 1.76 s |
| 64 × 64 × 1024 | 16 | – | 192 ms | 1.98 s |

## The wide plane, 16384 × 16384 × 32

No layout keeps up with an animation at 5 frames a second, and regions take about the same time in
every layout. With several users every plane read chunks another user read, so planes are shown for
one user only. With eight users and the default cache the 16-deep layouts' animations ran past the
bench's 600-second limit; 69 of their operations, the only failures in the measurements, are left out
of the animation table, whose other operations were 97% late.

### One user

| chunk | MiB | plane, mac | spectrum, mac | region, mac | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 193 ms | 5.8 ms | 399 ms | 637 ms | 27 ms | 1.53 s |
| 256 × 256 × 4 | 1 | 469 ms | 2.4 ms | 390 ms | 2.51 s | 9.5 ms | 1.43 s |
| 128 × 128 × 16 | 1 | 1.99 s | 1.4 ms | 391 ms | 9.42 s | 4.8 ms | 1.34 s |
| 1024 × 1024 × 1 | 4 | 216 ms | 24 ms | 463 ms | 520 ms | 86 ms | 1.22 s |
| 512 × 512 × 4 | 4 | 447 ms | 7.9 ms | 428 ms | 2.02 s | 26 ms | 1.26 s |
| 256 × 256 × 16 | 4 | 1.70 s | 4.3 ms | 409 ms | 7.75 s | 12 ms | 1.18 s |
| 2048 × 2048 × 1 | 16 | 251 ms | 93 ms | 385 ms | 525 ms | 309 ms | 1.25 s |
| 1024 × 1024 × 4 | 16 | 473 ms | 29 ms | 443 ms | 1.97 s | 91 ms | 1.22 s |
| 512 × 512 × 16 | 16 | 1.75 s | 16 ms | 423 ms | 7.44 s | 31 ms | 1.20 s |

### Four users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|
| 512 × 512 × 1 | 1 | – | 87 ms | 4.66 s |
| 256 × 256 × 4 | 1 | – | 22 ms | 4.26 s |
| 128 × 128 × 16 | 1 | – | 5.8 ms | 4.13 s |
| 1024 × 1024 × 1 | 4 | – | 250 ms | 3.83 s |
| 512 × 512 × 4 | 4 | – | 66 ms | 3.94 s |
| 256 × 256 × 16 | 4 | – | 18 ms | 3.60 s |
| 2048 × 2048 × 1 | 16 | – | 752 ms | 3.82 s |
| 1024 × 1024 × 4 | 16 | – | 238 ms | 3.89 s |
| 512 × 512 × 16 | 16 | – | 71 ms | 3.89 s |

### Eight users

| chunk | MiB | plane, acdc | spectrum, acdc | region, acdc |
|---|---|---|---|---|
| 512 × 512 × 1 | 1 | – | 171 ms | 9.21 s |
| 256 × 256 × 4 | 1 | – | 42 ms | 8.30 s |
| 128 × 128 × 16 | 1 | – | 11 ms | 7.78 s |
| 1024 × 1024 × 1 | 4 | – | 499 ms | 8.02 s |
| 512 × 512 × 4 | 4 | – | 130 ms | 7.89 s |
| 256 × 256 × 16 | 4 | – | 35 ms | 7.25 s |
| 2048 × 2048 × 1 | 16 | – | 1.03 s | 7.64 s |
| 1024 × 1024 × 4 | 16 | – | 427 ms | 7.87 s |
| 512 × 512 × 16 | 16 | – | 132 ms | 7.75 s |

### Animation

Cells as in the pancake's animation table: late with 5 GiB / late with 1 GiB, first frame.

| chunk | MiB | mac, 1 user | acdc, 1 user | acdc, 4 users | acdc, 8 users |
|---|---|---|---|---|---|
| 512 × 512 × 1 | 1 | 10% / 18%, 209 ms | 47% / 47%, 398 ms | 93% / 93%, 757 ms | 95% / 95%, 1.66 s |
| 256 × 256 × 4 | 1 | 22% / 97%, 468 ms | 22% / 93%, 1.42 s | 92% / 97%, 2.57 s | 96% / 97%, 3.32 s |
| 128 × 128 × 16 | 1 | 97% / 97%, 1.58 s | 97% / 97%, 5.53 s | 97% / 97%, 7.66 s | 97% / 97%, 10.3 s |
| 1024 × 1024 × 1 | 4 | 94% / 94%, 230 ms | 47% / 55%, 357 ms | 94% / 94%, 908 ms | 95% / 95%, 1.83 s |
| 512 × 512 × 4 | 4 | 22% / 97%, 433 ms | 22% / 97%, 1.24 s | 95% / 97%, 2.02 s | 95% / 97%, 3.35 s |
| 256 × 256 × 16 | 4 | 97% / 97%, 1.49 s | 97% / 97%, 4.74 s | 97% / 97%, 6.83 s | 97% / 97%, 11.1 s |
| 2048 × 2048 × 1 | 16 | 94% / 94%, 259 ms | 89% / 92%, 372 ms | 94% / 94%, 1.04 s | 95% / 94%, 2.56 s |
| 1024 × 1024 × 4 | 16 | 22% / 97%, 459 ms | 22% / 97%, 1.21 s | 89% / 97%, 1.88 s | 96% / 97%, 3.52 s |
| 512 × 512 × 16 | 16 | 97% / 97%, 1.57 s | 97% / 97%, 4.61 s | 97% / 97%, 7.02 s | 97% / 97%, 11.7 s |

## Checks

### Synthetic and real pixels

The first round of these measurements, on acdc with the same settings but for its sweep's own
operation counts, read a real 7763 × 4742 × 256 ASKAP cube (compression 1.62) rewritten into three
layouts, and the synthetic pancake (1.53) in the same ones. One user, medians of every operation:

| chunk | plane, real | plane, synthetic | spectrum, real | spectrum, synthetic | region, real | region, synthetic |
|---|---|---|---|---|---|---|
| 1024 × 1024 × 1 | 74.5 ms | 77.4 ms | 435 ms | 485 ms | 3.20 s | 3.03 s |
| 512 × 512 × 4 | 277 ms | 288 ms | 166 ms | 174 ms | 1.91 s | 1.88 s |
| 256 × 256 × 16 | 1.08 s | 1.11 s | 43.9 ms | 45.3 ms | 1.64 s | 1.63 s |

The synthetic cube compresses a little less, so it reads a little more and decodes about as much.

### Did `purge` empty the Mac's caches?

After its measurements, each layout's planes were read again with nothing emptied: fresh processes,
so no decoded chunk survived, but the page cache untouched. First touches only:

The pancake:

| chunk | MiB | plane, mac | plane warm, mac |
|---|---|---|---|
| 512 × 512 × 1 | 1 | 18 ms | 18 ms |
| 256 × 256 × 4 | 1 | 67 ms | 66 ms |
| 128 × 128 × 16 | 1 | 260 ms | 153 ms |
| 512 × 512 × 2 | 2 | 34 ms | 33 ms |
| 256 × 256 × 8 | 2 | 126 ms | 85 ms |
| 128 × 128 × 32 | 2 | 507 ms | 315 ms |
| 1024 × 1024 × 1 | 4 | 20 ms | 20 ms |
| 512 × 512 × 4 | 4 | 66 ms | 64 ms |
| 256 × 256 × 16 | 4 | 248 ms | 164 ms |
| 1024 × 1024 × 2 | 8 | 38 ms | 37 ms |
| 512 × 512 × 8 | 8 | 128 ms | 102 ms |
| 256 × 256 × 32 | 8 | 486 ms | 329 ms |
| 2048 × 2048 × 1 | 16 | 33 ms | 32 ms |
| 1024 × 1024 × 4 | 16 | 75 ms | 73 ms |
| 512 × 512 × 16 | 16 | 258 ms | 172 ms |

The cigar:

| chunk | MiB | plane, mac | plane warm, mac |
|---|---|---|---|
| 512 × 512 × 1 | 1 | 1.6 ms | 1.5 ms |
| 256 × 256 × 4 | 1 | 2.0 ms | 1.7 ms |
| 128 × 128 × 16 | 1 | 3.3 ms | 3.1 ms |
| 64 × 64 × 64 | 1 | 10 ms | 8.6 ms |
| 32 × 32 × 256 | 1 | 32 ms | 31 ms |
| 512 × 512 × 4 | 4 | 3.7 ms | 3.5 ms |
| 256 × 256 × 16 | 4 | 4.8 ms | 15 ms |
| 128 × 128 × 64 | 4 | 10 ms | 9.8 ms |
| 64 × 64 × 256 | 4 | 31 ms | 30 ms |
| 512 × 512 × 16 | 16 | 12 ms | 12 ms |
| 256 × 256 × 64 | 16 | 17 ms | 17 ms |
| 128 × 128 × 256 | 16 | 38 ms | 38 ms |
| 64 × 64 × 1024 | 16 | 120 ms | 82 ms |

The wide plane:

| chunk | MiB | plane, mac | plane warm, mac |
|---|---|---|---|
| 512 × 512 × 1 | 1 | 193 ms | 192 ms |
| 256 × 256 × 4 | 1 | 469 ms | 380 ms |
| 128 × 128 × 16 | 1 | 1.99 s | 1.30 s |
| 1024 × 1024 × 1 | 4 | 216 ms | 217 ms |
| 512 × 512 × 4 | 4 | 447 ms | 410 ms |
| 256 × 256 × 16 | 4 | 1.70 s | 1.23 s |
| 2048 × 2048 × 1 | 16 | 251 ms | 245 ms |
| 1024 × 1024 × 4 | 16 | 473 ms | 453 ms |
| 512 × 512 × 16 | 16 | 1.75 s | 1.25 s |

In the pancake's and the wide plane's chunks eight or more channels deep, whose planes read the most,
the warm read was 20 to 41% faster, which it could not have been had `purge` left the page cache full.
Elsewhere the two mostly agree within noise: the SSD delivers chunks about as fast as twelve cores
decode them.

## Reproducing

Write each layout with `generate.py`, for example:

```sh
tools/zarr-bench/generate.py --synthetic --shape frequency=256,polarization=1,l=7763,m=4742 --seed 1 \
    --chunk l=512,m=512,frequency=4 --codec blosc:zstd:5:shuffle --output cube.zarr
```

and measure it with the machine's thread counts and cache-emptying method (`--cold fadvise` on acdc;
`--cold command` with `--drop-cache-cmd /usr/sbin/purge` on the Mac, or with
`tools/zarr-bench/drop-lustre-cache.py` on Lustre), once with each cache, writing every cube's
results to its own CSV:

```sh
carta-zarr-bench run cube.zarr --mode plane,spectrum,region,animation --cache-bytes 5G \
    --processes 4 --io-threads 4 --decode-threads 28 --trials 3 --animation-prefetch --cold fadvise \
    --label "LAYOUT|users=4|cache=5G" --csv pancake.csv
carta-zarr-bench run cube.zarr --mode animation --cache-bytes 1G \
    --processes 4 --io-threads 4 --decode-threads 28 --trials 3 --animation-prefetch --cold fadvise \
    --label "LAYOUT|users=4|cache=1G" --csv pancake.csv
```

On the Mac, planes were also read again with `--cold off`, labelled `LAYOUT|users=1|cache=5G|warm`.
With each machine's CSVs in one directory:

```sh
docs/measurements/chunk-size/summarize.py mac=RESULTS_MAC acdc=RESULTS_ACDC lustre=RESULTS_LUSTRE > medians.csv
```

[`medians.csv`](chunk-size/medians.csv) is that file, one row per machine, cube, layout, user count,
cache and mode. The raw results, about 90 MB, were kept on the machines.

The real-pixel check ran `tools/zarr-bench/sweep.py` on acdc over the pancake's three 4 MiB layouts,
once on a 7763 × 4742 × 256 ASKAP cube that is not in the repository and once on the synthetic
pancake, with the same settings as above: 1 and 8 users, three trials, `--cold fadvise`, 4 I/O
threads, a decode thread per core and a 5 GiB cache.
