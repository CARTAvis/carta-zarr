# A dataset size says which question it answers, and nothing about how the two compare

`Dataset::Size` answers one of two different questions. Within its timeout it walks the store and
returns the sum of the sizes of the files it holds. When that walk fails, it reads the metadata
instead and returns the uncompressed size every array declares.

Until now the second was labelled `is_upper_bound = true`, and carta-backend-2 passed that through
to `FileInfo::size_is_upper_bound`, which the file browser renders as a `≤` in front of the number.
That label is a claim about how the two answers compare, and this library is not in a position to
make it. It now reports `SizeBasis::declared` or `SizeBasis::measured` -- which question was
answered -- and leaves the comparison to whoever knows the store.

This records why, because the shape invites two suggestions that both look obvious and are both
wrong.

## The comparison cannot be made, by construction

The declared answer is reached *because* the store could not be measured. There is no measurement
to compare it against; that is the precondition of the path.

Nor can the relationship be derived. Two effects push in opposite directions:

- Compression pushes what the store occupies **below** what its arrays declare.
- Per-node `zarr.json` documents and sharding indices are in the store but not in the declared
  total, pushing it **above**.

Which one wins is a property of the store that this path could not read. In this repository's own
fixtures the second wins by a wide margin -- `tests/data/images/zarr/xradio/minimal` declares 4,140
bytes and occupies 34,459, of which 33,123 is `zarr.json` -- so the declared answer is nearly an
order of magnitude *below* the store.

### Making it a true bound does not work either

The obvious repair is to add what is missing: the metadata bytes (the fallback already reads every
`zarr.json`, so they are nearly in hand), an allowance for compression expansion, and the shard
index, which is derivable from `grid_shape`.

That still fails, and not by a margin that can be closed. A store is a directory, and a directory
can hold files that are not Zarr's: an editor's leftovers, a sibling README, a half-written chunk.
Bounding what a directory contains requires enumerating it, and this path runs precisely when
enumeration failed. An upper bound on a container you cannot enumerate does not exist.

Narrowing the claim -- "an upper bound on this dataset's own arrays and metadata" -- would be true,
but it would no longer be a bound on the same quantity the measured answer reports, so the two could
not share a field or a column.

## The failure is not common, and that is not why this changed

The second suggestion is the opposite one: the declared answer is fine because a real image is a
compressed cube where the metadata is noise.

That is true, and stronger than it looks. `zarr.json` is written once per node, not once per chunk,
so the metadata is roughly a constant times the number of data variables -- tens of kilobytes. An
image small enough to be outweighed by it walks in microseconds, so it never reaches the fallback's
50 ms timeout in the first place. **The two conditions are close to mutually exclusive**: the stores
where the bound fails are the stores that get measured.

Reaching the bad label in production takes a small store *and* a walk that failed for a reason other
than time -- a directory that refused to be read.

The numbers above were obtained by forcing the fallback, with a zero timeout or an in-memory
transport. They show what the declared answer is worth, not how often anyone sees it. This is
recorded because it was got wrong once already, by measuring the fixtures and concluding the defect
was common without checking whether production reaches that path.

So the reason for the change is not frequency. It is that `DatasetSizeBytes` cannot tell which of
its four failures it hit -- an expired deadline says the store is large and the declared size is
almost certainly above it; an unreadable directory says nothing at all -- and a library that cannot
tell should not be the one asserting.

## What this does not decide

`FileInfo::size_is_upper_bound` and the frontend's `≤` were left as they were, for the protocol and
the frontend to decide. They have since followed this: the field is `size_is_declared`, the backend
passes the basis through without inferring a comparison, and the file browser marks a declared size
with `~` rather than `≤`.
