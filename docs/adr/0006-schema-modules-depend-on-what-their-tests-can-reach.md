# A schema module's dependency is chosen by what its tests can reach

Three rules that read metadata and nothing else -- the observation, the flag selection, and the beam
table assembly -- were file-local inside `schema/xradio/image.cc`, whose only entry point,
`DescribeImage`, reads coordinate values. Each now has its own module. The two that touch a store
take different dependencies, and the difference is deliberate rather than an inconsistency to tidy
up later.

`DetermineFlag` takes the whole `Store`. `DescribeBeams` takes the arrays already read.

## Considered options

Giving `DetermineFlag` a narrower interface -- the two methods it uses, `ReadArrayMetadata` and
`ListNodes` -- was the obvious alternative and was rejected. What varies underneath is the
transport, which has two adapters already: the filesystem one serves production and the in-memory
one in `tests/support/` serves the profile tests. A seam above that would have exactly one adapter,
`Store` itself, and a seam built for a caller that does not exist is guessed rather than designed.
The substitutability this module needs is a layer down and already paid for: `carta_zarr_profile_tests`
compiles `store.cc` and `transport.cc`, stubs the value and pixel readers, and links no TensorStore,
so a rule that asks a `Store` only for metadata is exercised there without a directory tree.

Giving `DescribeBeams` a `Store` for symmetry with it was rejected on the opposite evidence.
Reading a beam table is `ReadStringArray1D` for the labels and `ReadNumericArray` for the values,
and both of those are *value* reads: they resolve an array directory and go to the filesystem. They
are the two functions the profile-test build replaces with stubs that can only fail. So a
`DescribeBeams` that read its own arrays would be unreachable in the one build that has no
TensorStore, and its four tests would stay what they were -- a store written to a temporary
directory, opened, and read back, to check index arithmetic over a flat buffer.

## Consequences

The rule that decides the shape is: **a module takes the narrowest thing its tests can stand up.**
For flag selection that is a `Store`, because a `Store` over an in-memory transport is cheap. For
the beam table it is two vectors, because the arrays behind them are not.

`ReadBeams` keeps its name and its place in `image.h`. What it does now is three store reads and one
call; everything it used to decide -- which dimension is which, that the labels name a major axis, a
minor axis and a position angle, that a table missing a dimension is an error rather than an empty
list, and that time varies slowest -- moved into `DescribeBeams` and is checked from two vectors.
`DescribeBeams` takes the table's node name as well, only because `ArrayMetadata` does not carry the
name of the node it came from and the error should say which variable it is about.

`RequireUsableFlag` is inline in `flag.h` rather than beside `DetermineFlag` in the translation
unit. This is the part that looks like a style choice and is not: `flag.cc` references `Store`, so a
test target compiling it has to link `store.cc`, `transport.cc` and the stubs. Inline, the three
ways a mask is refused -- not marked as a flag, not boolean, not the image's own dimensions in the
image's own order -- are checked by a target that links nothing at all. Those three were checked by
opening three stores on disk.

`IsFlag` moved with the flag module, because it asks whether a variable is a flag. The rest of the
taxonomy in `image.cc` -- `HasAllAxes`, `KnownImageRank`, and the three axis-name lists that
disagree about order -- did not move, and is a separate question.

What the move is worth, in the terms it was chosen by: `schema_probe_test.cc` loses six cases and
about a quarter of its length, and what is left there is what belongs there -- a store opens, it
reads, and a failure below reaches the caller. `TestBeamTableWithUnreadableLabels` stays for exactly
that reason: it checks that a label array which cannot be read is reported rather than swallowed,
which is a store failure and not arithmetic.

One thing this does not change: `DescribeObservation` keeps its signature, returning a value rather
than a `Result` and reporting no diagnostic. Optional metadata that cannot be read is skipped and
the image still opens. Making it speak up is a change to the contract, and this was a move.
