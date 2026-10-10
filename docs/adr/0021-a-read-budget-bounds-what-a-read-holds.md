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

1. `Image::Read` was cut only when somebody watched it or stated a budget. carta-backend's plane
   reads did neither, so one plane sent all 160 of its chunks, 10 GiB decoded, to TensorStore at
   once.
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
- **Every read is held to its budget**, the caller's or the library's, whether or not it is watched,
  in pieces of at most what the budget affords -- with one exception, below. A read is not a fixed
  number of chunks in flight however large it is: TensorStore decodes no more at once than the
  context has decode threads, but it reads every chunk the request spans as fast as its I/O threads
  allow, and when the decoders fall behind the compressed chunks queue for them. On the 28-thread
  machine, one unwatched read of 128 whole planes of the 7763 x 4742 test cube in 256 x 256 x 16
  chunks, issued whole with two I/O threads and no cache, held beside its destination 292 MiB on
  twenty-eight decode threads, 5.1 GiB on four and 9.6 GiB on one -- dev alike -- and 4.7 to 9.2 GiB
  with a flag, whose folded-in mask is a byte an element of the whole read. With the I/O threads left
  to TensorStore's default the queue grows on twenty-eight threads too (3.0 GiB).
- **A read nobody watches is cut into layers of its chunks**, a chunk deep along the spectral,
  polarization and time axes and the whole of its sky plane, or into what its budget affords when
  that is more, once its budget affords a chunk for every decode thread. That is the exception: a
  layer of the test cube is 589 chunks, more than any default budget affords, and cut finer each piece
  spans the plane's chunks only
  in part and is gathered in segments (below), which ran the 128-plane read 25-65 % slower (9.0 s
  against 7.0 s on twenty-eight threads, 12.4 s against 7.5 s on four). In layers it ran as fast as
  issued whole and held 277 MiB on twenty-eight threads and 960 MiB on four -- 0.8 and 2.0 GiB with
  the flag -- bounded by the plane and the chunk depth rather than by the read's extent. The cut goes
  on the slowest axis the read selects more than one of, so this holds for a read of one Stokes
  parameter, as carta-backend's are; one across several Stokes parameters and many channels is cut
  along polarization and, a polarization chunk deep being more than a layer, falls to segments as it
  did before. A plane is
  one layer and is read whole, as before, at the same speed. A watched read is always cut to its
  budget, so that it has pieces to report. Reading segments straight into the destination rather
  than gathering them would let a layer be cut to the budget at no cost, and is left to a later
  change.
- **A piece too large at one chunk deep is read in segments**, cut along the axes below the piece's
  down to a single chunk, each gathered in a buffer of the library's and put in place. Progress is
  still reported a piece at a time, so the finished part stays a prefix, which carta-backend's cursor
  spectrum relies on. A gathered segment is a copy of what its chunks decode to and a fraction of what
  they hold, so it costs nothing measurable beside the decode.
- **Under a budget smaller than one chunk, a read decodes one chunk at a time** and holds what that
  chunk holds. It used to refuse a masked read with `buffer_too_small` when the folded-in flag of its
  least piece would exceed the budget, while the unmasked read beside it went on to decode a gigabyte;
  the flag is part of what a chunk holds now, and a caller whose data is chunked that way still gets
  its pixels.

## What it does not do

The budget is per read. Reads running at once each hold their own, and a process-wide bound across
them is a separate question with its own costs -- which read waits, and whether a scan may starve the
tiles a user is looking at.

It does not make a large chunk cheaper to read through. A plane of 512 x 512 x 64 chunks still decodes
sixty-four planes' worth to return one; what changes is that it no longer holds them all at once.

ADR 0005 keeps `Image::Read` outside the pass; this is what the two share.
