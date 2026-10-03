# An array is held to its own document once, by the store, and an image when it opens

A store that consolidated its metadata describes every array twice: once in the root's copy, which
is what the store reads, and once in the array's own `zarr.json`, which is what the values are
decoded from -- TensorStore opens an array from its own document, and so does the label decoder. The
two are meant to be the same document and are not always: a copy left stale by a rewrite, or a store
written by hand. `Store::VerifyArray` holds the one to the other, once per node, and every value read
asks it before it asks where the array lives. Describing an image asks it of the image and its flag
before any value is read.

## What it replaced

The question was asked by each reader on its own. The numeric reader opened the array through
TensorStore and compared TensorStore's view of it with the store's metadata; the pixel reader did
the same on every piece of every read; the label reader, which does not use TensorStore, asked
nothing. Two fixes in a row (c809e95, 1462ee5) were the same check carried to one more reader, and
the third reader was still open: a copy saying one label a chunk over chunks of three read the first
label and then the fill value, and a label array whose directory was gone read as nothing but fill.

Nothing read an image's own array or its flag until a pixel read, so an image whose own document
disagreed with the copy opened and failed every read -- and one whose flag disagreed failed only
every masked one. The same disagreement surfaced at four different moments depending on which array
it was in.

## What is compared, and against what

The array's own document, read through the transport and parsed as any other, against the
`ArrayMetadata` the store holds -- not TensorStore's view of the opened array. That is what lets one
check serve the label decoder as well as TensorStore, and what lets it run in the profile tests,
which link no TensorStore: the in-memory transport serves a child document that disagrees with the
root's copy as easily as one that agrees. Comparing against TensorStore's view also needed an
exemption for a 1-D array naming no dimensions, because our parser takes `dimension_names` from the
attributes where TensorStore does not; comparing two of our own parses needs none.

Held to each other on what the store takes from the copy: the extent, the names and order of the
dimensions, the data type, and the chunks and shards a read is planned from. Not on how the chunks
are encoded -- codecs, chunk keys, fill value, a string's length -- because the copy is never used to
decode one. What `VerifyArray` hands back is the own document, and that is what a reader decodes
with. A stricter rule, the two documents equal but for attributes, was considered and not taken: a
child rewritten by a later zarr-python that spells a default out is the same array, and refusing it
would refuse a store every other reader opens.

A copy naming an array whose own document is missing is `invalid_metadata`, naming the array: ADR
0004 already calls such a store malformed, and said its error should name the missing document.

## Why an image is held to it when it opens

ADR 0009 is about qualification, and a stale copy is outside what qualification can see -- it reads
the copy. What an image that opens promises is that it can be read, so `DescribeImage` asks of the
image right after qualifying it, before reading any coordinate, and asks of the flag as soon as the
flag is chosen. The coordinates and the polarization labels are held to theirs as describing reads
them; a beam table and its labels, when `ReadBeams` reads them.

Probing and listing never ask. They are metadata, and consolidated metadata exists to save exactly
the reads checking would cost; a store reached over a network pays one round trip for each. A
listed image whose arrays disagree with the copy is therefore refused when it is opened, not left
out of the listing.

## What it costs

One read of the array's own document per node, once for the life of the store, and only for a node
the copy accounted for: any other node's metadata was read from its own document already, and there
is nothing to hold it to. Opening an image of a consolidated store reads the image's, its flag's and
its coordinates' own documents, five to seven small files, most of which TensorStore reads again a
moment later when it opens the coordinates.

Measured on this machine, release build, local disk with the page cache warm, `Context::Create`,
`Dataset::Open` and `OpenImage("SKY")` with nothing kept between them, tenth percentile of 400:
the conformance fixture went from 0.61-0.62 ms to 0.77-0.80 ms and the minimal fixture from 0.64-0.65
ms to 0.83-0.86 ms, over three interleaved rounds. The pixel fixture, which is not consolidated, did
not move: 0.97-1.02 ms before and after. The copy still saves
more than checking spends, on the one image opened; the listing keeps all of it.

A pixel read asks once per piece, and asking costs one lookup keyed by the name as asked. Keyed
through `ReadArrayMetadata` it normalized the name once more per piece, and a whole-cube read came out
2-3% slower; looked up directly, the paired runs straddle zero.

## Not changed

A backslash in a variable's name. The store's name rule accepts one and the filesystem transport
refuses it when asked where the array lives, which looked like the second opinion a transport is not
allowed. It is TensorStore's: its file kvstore reads a backslash as a separator. Such a variable is
listed and opens, and every read of it is refused saying why. Listing it as unopenable would need the
store to know what TensorStore can locate, for a name no XRADIO writer produces.
