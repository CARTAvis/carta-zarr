# A read budget bounds what a read holds beyond its destination

`ReadOptions::read_budget_bytes` says how much memory one read may hold at once, not counting the
caller's own destination, and every read keeps to it: `Image::Read`, `Image::Prefetch` and every pass.
A chunk is charged what it holds while it is read, which is about three times what it decodes to, and
a read holds at least one chunk however small its budget.

## What went wrong without it

On an image of large chunks, carta-backend computing a cube histogram of a Zarr store exceeded 20 GB
of RAM where the FITS file it was made from stayed under 400 MiB. Measured on a 7763 x 4742 x 128
float32 cube in 512 x 512 x 64 chunks -- 64 MiB decoded each, blosc zstd, 11 GB on disk -- on a
28-thread Linux machine:

| What | Peak resident |
|---|---|
| `Image::Read` of one plane, 147 MiB of pixels | 2.6-4.7 GB |
| the same with a 64 MiB budget | 1.4 GB |
| a cube histogram at the default budget (256 MiB) | 0.9-1.0 GB |
| carta-backend, open and cube histogram | 2.5-2.7 GB |

Three things made those numbers, each on its own enough to defeat a budget.

1. `Image::Read` was cut only when it was given a progress callback or a budget. carta-backend's
   plane reads did neither, so one plane sent all 160 of its chunks, 10 GiB decoded, to TensorStore
   at once.
2. It was cut along one axis, the slowest it selects, so its least piece was a whole row of chunks:
   sixteen here, about 1 GiB, whatever the budget said.
3. The budget counted the bytes a chunk decodes to. A chunk being decoded also holds its compressed
   bytes and the codec's own buffer, and the measured peak fit
   `base + destination + chunks x (2 x decoded + compressed)` -- about 3.5 times the budget.

## What a read holds is decided by its chunks, not the machine

The same cube histogram at one, two and four chunks a read peaked at 282 MiB, 457 MiB and 894 MiB.
With four and twenty-eight decode threads the peaks were within 3% of each other; one thread held
20-25% less, because fewer chunks were in flight at once. More threads never held more. So a bound
stated in chunks a read holds on any machine, and the thread count only lowers what it reaches.

Compressed bytes are at most about the decoded ones, so a chunk holds at most three times what it
decodes to. `kHeldBytesPerDecodedByte` is that three, on every image. A figure per codec would be
closer for a codec that holds less, and would go quietly wrong when TensorStore's codecs change;
three is a statement of what the budget means, and nothing overrides it the way ADR 0014's tunings
can be.

## The rules

- **A chunk is charged what it holds**: three times what it decodes to, the flag it brings included
  in the flag's own chunks, plus the buffers the library allocates for its share of the read -- the
  folded-in flag at a byte an element, and the pixels at four bytes an element wherever the library
  rather than the caller holds them. `ReadCost` answers both ways (`PixelsHeld`).
- **The caller's destination is not charged.** It exists before the call, at a size the caller chose,
  and no cutting makes it smaller; charging it would buy fewer chunks and save nothing.
- **The default budget aims at two chunks for every decode thread of the image's context, held
  between 256 MiB and 2 GiB.** Two a thread keeps every thread a chunk ahead of the one it is
  decoding. A count fixed for every machine does not: eight a read, the first aim, was measured on a
  five-core machine, and on twenty-eight threads it left most of them idle and plane reads of the
  64 MiB chunks above ran 2.4 times slower. 256 MiB keeps small chunks to at least sixty-four a read,
  as before; 2 GiB holds the 64 MiB chunks above to about ten on twenty-eight threads, which trades
  some of their speed for the bound.
- **Every read is held to its budget**, the caller's or the library's, with a progress callback or
  without, in pieces of at most what the budget affords -- with one exception, below. A read is not
  a fixed number of chunks in flight however large it is: TensorStore decodes no more at once than
  the context has decode threads, but it reads every chunk the request spans as fast as its I/O
  threads allow, and when the decoders fall behind the compressed chunks queue for them. On the
  28-thread machine, one read with no progress callback of 128 whole planes of the 7763 x 4742 test
  cube in 256 x 256 x 16 chunks, issued whole with two I/O threads and no cache, held beside its
  destination 292 MiB on twenty-eight decode threads, 5.1 GiB on four and 9.6 GiB on one -- dev
  alike -- and 4.7 to 9.2 GiB with a flag, whose folded-in mask is a byte an element of the whole
  read. With the I/O threads left to TensorStore's default the queue grows on twenty-eight threads
  too (3.0 GiB).
- **A single plane read with no progress callback is read whole**, once its budget affords a chunk
  for every decode thread. That is the exception. A plane of the test cube is 589 chunks, more than
  any default budget affords, and cut to the budget it is read in segments one after another, which
  ran a quarter slower -- 0.33 s against 0.26 s on twenty-eight threads -- while holding no less:
  what a plane holds is bounded by the plane itself, its decode threads' chunks and the queue ahead
  of them and a flag of one channel, about 130 MiB on twenty-eight threads and 370 MiB on four. A
  plane is what a channel change and an animation wait on. Any read of more than one plane is cut to
  its budget: it was cut into layers of its chunks, a chunk deep along the spectral, polarization
  and time axes and the whole plane, because a segment used to be gathered in a buffer of the
  library's and copied into place, which ran the 128-plane read 25-65 % slower; read straight into
  place, it runs within 4-18 % of the layers, and holds 90 MiB rather than 950 MiB on four threads,
  90 MiB rather than 2.3 GiB with the flag, and 230 MiB rather than 830 MiB with the flag on
  twenty-eight. The layers' memory grew with the I/O threads, too: with TensorStore's default,
  1.4 GiB on twenty-eight threads where the budget holds 420 MiB. A read with a progress callback is
  always cut to its budget, so that it has pieces to report, and to sixteen chunks a piece unless
  one chunk along the cut is more, so that they come often: a budget aimed at the decode threads
  grows with them, and a cursor spectrum of the cigar test cube cut to it alone came in nine pieces
  on 28 threads instead of thirty, 160 ms apart on a cold local disk rather than 55, and further
  apart on slower storage.
- **A piece too large at one chunk deep is read in segments**, cut along the axes below the piece's
  down to a single chunk, each read straight into its place in the destination through the strides
  of the piece it is part of. Progress is still reported a piece at a time, so the finished part
  stays a prefix, which carta-backend's cursor spectrum relies on. A segment's pixels are the
  caller's, so the budget pays only for the chunks it decodes and the flag beside them. It was
  gathered in a buffer of the library's and copied into place, and that copy, transposing each
  chunk into a small dense buffer, nearly doubled the CPU the decode threads spent on four threads.
- **Under a budget smaller than one chunk, a read decodes one chunk at a time** and holds what that
  chunk holds. It used to refuse a masked read with `buffer_too_small` when the folded-in flag of its
  least piece would exceed the budget, while the unmasked read beside it went on to decode a gigabyte;
  the flag is part of what a chunk holds now, and a caller whose data is chunked that way still gets
  its pixels.

## A reduction's reads are not cut as a watched read's are

*Added for 0.3.0 ([#73](https://github.com/CARTAvis/carta-zarr/issues/73)). The rules above are
unchanged.*

A reduction calls its sink and checks for a stop between reads, so its reads decide how often a
caller hears from it. They are cut only to the budget, never to sixteen chunks as a watched read is,
and that was measured before being kept. One user, the region of a spectral profile over the cigar
(512 x 512 x 30000 in 128 x 128 x 64 chunks) and pancake (7763 x 4742 x 256 in 256 x 256 x 16) test
cubes, cold, against an `emit_every_channels` that gives reads of about sixteen chunks, the size such
a cut would:

| | Default: longest wait between sink calls | Sixteen chunks: longest wait | Sixteen chunks: total time |
|---|---|---|---|
| Local NVMe, 4 I/O and 28 decode threads | 44-65 ms | about 20 ms where a read was more than sixteen chunks | up to +23% |
| Local NVMe, 2 and 4 | 25-40 ms | unchanged, but 22 ms for a small pancake region | unchanged, but -18% for that region |
| Lustre, 4 and 32, small regions | 0.15-0.30 s | 0.04-0.20 s | -26% to +9% |
| Lustre, 4 and 32, cigar, a quarter of the plane or more | 0.29-0.34 s | 0.12-0.15 s | +23-30% |
| Lustre, 4 and 32, the whole pancake plane | 0.31-0.34 s | unchanged | unchanged |
| Lustre, 2 and 4 | 0.13-0.23 s | unchanged, but 0.07 s for a small pancake region | unchanged, but -36% for that region |

A budget aimed at the decode threads makes a read larger as decoding gets faster, so from the local
disk a read took about the same 25-65 ms on four threads and on twenty-eight of one machine, and with
few threads its floor already makes reads of about twenty chunks. A machine whose threads decode more
slowly, or a costlier codec, would make each read longer; of the configurations measured, only
storage slower than decoding stretched the wait, to about 0.3 s on Lustre, and there the cut to
sixteen cost a large region a quarter of its time, because fewer chunks were in flight. A stop took
effect within 0.1 s everywhere, cut or not.

A region whose chunk layer is more than a read finishes a block only once a layer: for the whole
pancake plane, every 0.8-1.4 s on NVMe and every 2.6-5 s on Lustre, whatever the cut, with its
unfinished block handed over once a read. That is the plane's size, not the read's.

A caller that wants results sooner still can ask for them with a smaller `emit_every_channels`, which
can shorten the reads with the blocks, unless they were already that short. What that does depends
on the region and the storage, as above: later totals for large regions, and sooner results and
sometimes a shorter total for small ones.

## What it does not do

The budget is per read. Reads running at once each hold their own, and a process-wide bound across
them is a separate question with its own costs -- which read waits, and whether a scan may starve the
tiles a user is looking at.

It does not make a large chunk cheaper to read through. A plane of 512 x 512 x 64 chunks still decodes
sixty-four planes' worth to return one; what changes is that it no longer holds them all at once.

ADR 0005 keeps `Image::Read` outside the pass; this is what the two share.
