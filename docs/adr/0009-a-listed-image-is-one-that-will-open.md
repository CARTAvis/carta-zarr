# A listed image is one that will open, and one module decides it

`ImageEntry::openable` is the only thing a consumer has to go on when it decides which images to
offer. carta-backend-2 reads exactly that field to build a file's HDU list, so it is the flag that
puts a variable in front of a user.

It was decided by one rule and contradicted by another. Discovery published a variable as readable
on its dimension names, its data type and its `type` attribute; `DescribeImage` then refused the
same variable because its extent disagreed with the coordinate it named. The field was a summary of
what a variable looked like, and it was read as a promise that the variable would open.

It is that promise now. One module decides whether this profile will open a variable, everything
that needs the answer asks it, and the field was renamed from `readable` so that the word matches
the question -- `RequireOpenableDataset` already spelled it that way, and `readable` was the one
place saying the same thing with a different word.

## The rule existed twice, and the copies differed

Coordinate agreement -- an image's length along an axis must be the length of the dataset coordinate
of that name -- was implemented twice, because it was written for two different questions:

| | `ProbeReport::RequireCoordinateOf` | `RequireMatchingCoordinates` |
|---|---|---|
| asked about | the default image only | the image being described |
| walked | the five sky axes | the image's own `dimension_names` |
| also checked | rank, dimension name, data type kind | nothing else |
| coordinate that cannot be read | failed the probe | skipped |
| outcome | probe reports `invalid` | `invalid_metadata`, at open |

Classification existed twice as well. Discovery required the whole sky axis set, a real data type
and a `type` attribute that is not `flag`; `DescribeImage` opened by requiring `l`, `m`, a real data
type and not a flag. Those answered the same question and did not agree. Nothing reached the weaker
one, because `RequireOpenable` ran first, which is to say the divergence was invisible rather than
absent.

## What moved, and what stayed

The split is between what is true of one variable and what is true of the dataset.

**Per image, in `schema/xradio/qualification`**: classification, agreement with each coordinate the
variable names, and the reason for every refusal. These decide `openable`.

**Per dataset, in `ProbeReport`**: that a required coordinate array exists, that it is
one-dimensional, that it names itself, that it holds the right kind of data, and that
`coordinate_system_info` describes a direction coordinate. These decide whether the store matches
the profile at all.

`RequireCoordinateOf` split along that line exactly: its length comparison went, the rest stayed. It
still walks the default image's axes rather than the five unconditionally, because an axis the image
does not carry is deliberately not required -- that is why a continuum image with no frequency
coordinate opens. Nothing is lost by giving up the length there: the default image is the first one
that qualified, so it agrees with every coordinate it can read before the probe sees it.

Ordering and default selection stayed where they were. Ranking `SKY` first is presentation, and the
default is "the first qualified image", which is a consequence of qualification rather than part of
it. `DescribeAxes`'s logical axis order is untouched: it is a third list of axis names, but it
describes the order this library reports, which is a public contract and a separate question.

## Both directions, not one

A listing that offers only what will open is half of it. The other half is that nothing outside the
listing opens, and it is the half that is easy to lose: `Store::ReadArrayMetadata` reaches a node by
name whether or not the listing mentions it, so a qualification that read the name straight from the
store would have opened a variable written after the dataset was opened, which the listing a
consumer is holding does not have.

So qualification decides against `Store::ListNodes`, which is the snapshot the dataset's images were
enumerated from. `TestMetadataCache` is what says so, and it is not incidental to this: a Store is a
read-only view, and an image that appears in a view of a store that has already been listed is a
different thing from an image in it.

## A refusal says which kind of refusal it is

Qualification distinguishes two things the code did not:

- **This library does not do that.** A complex sky-plane variable, an aperture-plane variable. The
  store is well formed and the answer would not change if it were rewritten.
- **This store is malformed.** A variable whose extent disagrees with its coordinate, a node whose
  metadata will not parse as an array.

Three things follow from the distinction, and the first is why it is needed at all.

**A store with nothing openable in it is `invalid` when something in it was malformed, and
`no_match` otherwise.** A store whose only image disagreed with a coordinate was reported `invalid`
before, by the probe's own copy of the rule. Without the distinction it would have become `no_match`
-- every image disqualified, no default, nothing matched -- and a caller would be told that no
built-in profile recognised a store that a profile recognised perfectly well and found broken. A
store whose only image will not parse now answers the same way, which is a change: it used to be a
non-match with its reason dropped.

**A refusal reports the code its own reason implies.** `RequireOpenable` wrote
`unsupported_data_type` for every refusal; a coordinate disagreement reports `invalid_metadata` and
a complex variable reports `unsupported_data_type`, because the thing that knows the reason is the
thing that answers.

**The reason is put first.** `RejectionMessage` takes the first diagnostic, and node enumeration
order would otherwise decide which that is -- alphabetically, so a store holding both an
aperture-plane variable and a malformed one named the aperture plane. Diagnostic order carries no
contract, so this makes the data fit an assumption the reporting side already made.

The kind does not reach `ImageEntry`, and it does not reach `ImageDiscovery` either: it is a field
on the qualification of one node, and `InspectImages` reads it to decide one match kind. Promoting
it to shared vocabulary would be shaping a concept for a second profile that does not exist -- the
thing ADR 0006 refused to do for `DetermineFlag`.

### It needed a valid group to stop looking malformed

`Store::ListNodes` returns every node, and `ParseArrayMetadata` refuses a group with "Zarr node is
not an array". Discovery turned that into an `unreadable_array` diagnostic, so **every store with a
nested group emitted one**. A promotion rule that counted it would report a perfectly well-formed
store as invalid.

Qualification now checks `node_type` before treating a parse failure as a malformation: a group is
passed over in silence, and only a node claiming to be an array and failing to parse is malformed.
That is a defect on its own -- a valid group was diagnosed as a broken array whatever was decided
here -- so it landed as its own change, before the rule that depends on it. Naming a group where an
image is expected answers `not_found` now, the way naming a flag or a coordinate already did.

## There is one place that refuses

`RequireOpenable` existed because two callers had to agree: the facade asked its kept listing, and
the profile asked a store it had just inspected. Once `DescribeImage` asks qualification about the
one variable it was handed, there are no longer two callers, and the shared function is the
duplication it was written to prevent. It is gone, and with it the `Describe` / `DescribeVerified`
pair -- "Verified" named a precondition the caller no longer establishes, so one `Describe` is left.

One thing had to be carried over rather than dropped. `RequireOpenable` answered for a variable that
was never listed with `not_found` and "Image variable was not found". Reaching the store instead
yields the transport's own `not_found` -- "Zarr node is missing zarr.json", against a filesystem path
-- which leaks storage into an answer about images. Qualification translates it.

ADR 0007 recorded the profile's surviving members as `For`, `Probe`, `Discover`, `Describe`,
`DescribeVerified` and `ReadBeams`, and placed `RequireOpenable` beside `RequireOpenableDataset` as
the profile's question rather than the facade's. Both of those readings were right and neither is
disturbed here: the question stays the profile's, and it is now answered one level further in, by
the module that also decides what an image is. `DescribeVerified` leaves that list;
`RequireOpenableDataset` and `RejectionMessage` do not move.

## Why the module takes a `Store`

By ADR 0006's rule: a module takes the narrowest thing its tests can stand up. Qualification reads
array metadata and nothing else, and a `Store` over the in-memory transport is cheap --
`carta_zarr_profile_tests` links no TensorStore. The coordinate-agreement rule, which could only be
exercised against a directory tree, is a map entry now.

Cost is not a reason to hesitate. `Store::ReadArrayMetadata` is memoized, a dataset carries at most
five coordinates, and discovery already reads every node's metadata. Qualifying every image rather
than one costs comparisons, not reads -- and `DescribeImage` asking again about its own variable
costs a table lookup.

## Considered options

**Handing `DescribeImage` a token instead.** Qualification could return a value only a qualified
variable can produce, making the precondition structural rather than checked. It was rejected for
its reach: `SchemaProfile`'s dispatch is function pointers taking a name, and `Dataset::Impl` keeps
the public descriptor rather than the internal discovery, so a token means changing the profile's
dispatch and the facade's cached state. Asking again, against one memoized lookup, gets the same
guarantee without either.

**Leaving the refusal in `DescribeImage` and only sharing the rule.** This keeps every error code
and every test as it is, and keeps a listing that offers an image the library will refuse. It is the
contradiction, minus the duplication.

**Promoting every malformation to `invalid`.** A `MODEL` that disagrees with a coordinate would
close a dataset whose `SKY` is fine. CONTEXT.md already decides this the other way: a variable that
cannot be opened is diagnosed rather than refused, because the rest of the dataset is still
readable.

## Consequences

`openable` is a promise. A consumer that lists only openable images -- which is what a file browser
does -- stops offering one that will fail. carta-backend-2's `FileInfoLoader` changed with the
rename, in the same change, because the package-consumer test reinstalls the prefix it builds
against.

`DescribeImage` lost its opening classification block. It is reached for a variable qualification
has accepted, so re-deciding would be either duplication or disagreement.

`TestOneRuleDecidesWhatIsOpenable` lost its subject and kept its point as
`TestARefusalCarriesTheVariablesOwnReason`: the property worth holding is that a refusal carries the
variable's own diagnostic rather than a generic message, and it is asserted against `Describe`.
`TestImageDisagreeingWithACoordinate` became two tests with two subjects -- the rule, in memory, and
the three steps between it and a consumer, on disk. `TestFirstFaultIsTheOnlyDiagnostic` kept its
subject and changed its fixture: it is about the report latching after the first unmet requirement,
and one of the two faults it used no longer reaches the report.

ADR 0006 closed with a note that `HasAllAxes`, `KnownImageRank` and the three axis-name lists that
disagree about order had not moved and were a separate question. This is that question, and it
answers it in two parts: `HasAllAxes` was qualification all along, and the other two were never
about qualification at all.

## What this does not decide

**The other readings of a parse failure.** `DetermineFlag` ignores one and `TotalArraySizeBytes`
inspects `node_type` itself, so callers of `ListNodes` still reconstruct what a node is in more than
one way. Qualification stopped being one of them; the listing still hands out names.

**Whether a consumer should see why an image was refused.** The refusal's kind stays inside the
profile. Putting it on `ImageEntry` is a change to the public descriptor and should arrive with a
consumer that wants it.

**Diagnostics dropped on `no_match`.** `ProbeStore` fills neither `diagnostics` nor `schema_id` on
the branch where nothing matched, so a store of nothing but complex variables is still reported as
unrecognised with no reason attached, although `InspectImages` says in a comment that the discovery
diagnostics explain it. The promotion rule above rescues the malformed half of that case and leaves
the capability half exactly as it is.

*Since decided:* `ProbeStore` now keeps what the profile said when it did not match, so that store is
refused with the profile's own reason ("Complex sky-plane variables are not openable") instead of
"No built-in schema profile matched the Zarr store". Only diagnostics the profile produced are kept
-- a store with nothing in it the profile recognised still carries none, and falls back to the
generic message -- and `schema_id` stays empty, as `ProbeResult` promises for a kind that is not a
match. The kind itself does not change: the store is still one this library is not for. With a
second profile, the first to say anything would be the one heard, the same rule the invalid branch
already follows.
