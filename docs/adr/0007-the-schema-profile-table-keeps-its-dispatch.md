# The schema profile table keeps its dispatch with one entry in it

`SchemaProfile::BuiltIn()` holds exactly one profile, the XRADIO image one. Everything above it is
shaped for a table: `For` looks an identifier up and reports `unsupported_schema` when it misses,
`ProbeStore` asks every entry and refuses a store that more than one of them claims, and each
question a profile answers goes through a function pointer.

Read on its own, that is a seam with one adapter, and this repository has a rule about those --
ADR 0006 rejected a narrower interface for `DetermineFlag` because "a seam built for a caller that
does not exist is guessed rather than designed." Applying that rule here says to collapse the
dispatch: let the facade call `xradio::` directly and keep `RequireOpenable`.

This records why that reading is wrong, so that it is not re-derived every time someone reads
`profile.cc` and counts the entries.

## The caller exists, and it is the public interface

ADR 0006's rule is about *substitutability* seams: an interface introduced so that something can be
swapped, with nothing to swap it for. The thing it protects against is inventing a shape for a
consumer nobody has.

The schema profile is not that. A consumer already names a schema:

- `ProbeSchema(location, schema_id)` takes one from the caller.
- `DatasetDescriptor::schema_id` and `SchemaProbeResult::schema_id` report one.
- `kXradioImageSchema` is a public constant, which only means something if it is one of several.

So the identifier is shipped, and with it the promise that asking about a schema this library does
not have is a defined answer rather than a compile error. `schema_profile_test` pins it:
`SchemaProfile::For("future.schema")` reports `unsupported_schema`. Collapsing the table would leave
that promise standing with a string comparison behind it, which is the same dispatch with the name
taken off.

CONTEXT.md also names the concept -- "Schema profile: a named, versioned description of how an image
dataset is laid out, which the library matches a store against" -- so the module is honest about
intent rather than speculative. A second profile is expected; this is not a seam waiting for a
caller, it is a table waiting for a row.

## What was wrong with it, and what changed

The complaint underneath the count was real: the interface was nearly as large as the body, and part
of what a profile is for was being done by its callers.

Two schema rules lived in the facade. `RejectionMessage` decided that a probe's own diagnostic is a
better thing to show a consumer than "not a supported dataset", and the facade decided, per call
site, what a probe's refusal meant: nothing matched, something matched and was malformed, or a
profile matched a store with no images in it. Those are three different errors about a schema, and
they were being distinguished by the thing that had asked the question rather than by the thing that
had answered it.

`RequireOpenableDataset` is those three, beside `RequireOpenable`, which is the same question one
level down about an image within a dataset. `RejectionMessage` moved with them, since it is the rule
for turning any profile's refusal into something a consumer can read.

Three members went the other way. `SchemaProfile::id()` had no callers at all. `Inspect` is private:
it is the one enumeration `Probe` and `Discover` are each half of, and the only caller that wants
both halves is `ProbeStore`, which is a friend and reaches the entry directly. What is left is
`For`, `Probe`, `Discover`, `Describe`, `DescribeVerified` and `ReadBeams`, and the pass-throughs
among those are the dispatch itself rather than a layer over it.

## Consequences

A second profile supplies three functions -- inspect, describe, read beams -- and one identifier, and
nothing above it changes. That is the whole reason the table is kept.

It also turns on two things that cannot be exercised today. `ProbeStore`'s ambiguity branch is
unreachable while there is one entry, so nothing tests it; its diagnostic now takes its code name
from `ErrorCode::ambiguous_schema` rather than from a literal beside it, so at least the two cannot
drift apart in the meantime. And `ErrorCode::ambiguous_schema` is still emitted nowhere, because
`ProbeStore` answers with a `ProbeResult` rather than an `Error` and `RequireOpenableDataset` maps
every `invalid_dataset` to `invalid_metadata`. Distinguishing an ambiguity there would mean matching
on a diagnostic's code name, which is worse than leaving it; the right time to fix it is when a
second profile makes the branch reachable, and it should arrive with the test that reaches it.

`ErrorCode::unsupported_schema_version` is unemitted for the same kind of reason: `schema_version` is
reported but nothing refuses a store for carrying the wrong one. That is a gap in what profiles do,
not in how they are dispatched, and this ADR does not close it.

What would reopen this decision is the opposite of what opened it: if the second profile does not
arrive, or arrives as a version of the XRADIO one that `InspectImages` handles internally rather than
as a row of its own, then the table has one entry for a reason that has expired, and collapsing it
becomes the honest move.
