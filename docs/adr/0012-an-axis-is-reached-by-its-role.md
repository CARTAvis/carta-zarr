# An axis is reached by its role, and where it sits is the profile's

An image's axes are reported in a logical order, and a read's ranges follow it. Which order that is
was said two ways at once. The README promised `l`, `m`, `frequency`, `polarization`, `time` as part
of the API, and `kXradioImageAxisOrder` exported it for a consumer to index by. `descriptor.h`, next
to that same constant, said the order was the schema profile's choice and that a consumer which
must not assume should read `descriptor().axes` instead. [ADR 0009](0009-a-listed-image-is-one-that-will-open.md)
called it a public contract in passing.

The order is the profile's. `kXradioImageAxisOrder` left the public header, the XRADIO profile
keeps its order as `kLogicalAxes`, and `AxisIndex(axes, role)` is how a consumer finds an axis. The
index it returns is also where that axis's `Range` goes in a `ReadRequest`, so nothing about how a
request is written changed except where its indices come from.

## Nothing that reads axes was using the promise

The library's own walks never indexed by position. `AxisMap` is shared by every walk that reads
planes because they "all have to find x, y and the spectrum without assuming they are the first
three", and it finds them by role. carta-backend-2 maps a Zarr image onto CARTA's four axes in one
place, `CartaZarrAxes`, and does it by role too; it never named the constant. The only readers of
`kXradioImageAxisOrder` were the profile that produced the order and a conformance test that
pinned it.

So the promise constrained nothing a consumer did, and it would have constrained the library.

## What the promise would have cost

**A second profile.** [ADR 0007](0007-the-schema-profile-table-keeps-its-dispatch.md) keeps the
profile table's dispatch because the public interface already names a schema, and asking about one
this library does not have is a defined answer. A second profile would
either report XRADIO's order whatever its data looked like or break a contract every consumer had
been told to rely on. The constant's own comment named it for the profile "because the order is the
profile's choice", which is the same conclusion reached from the other end.

**An image with fewer axes.** `DescribeAxes` skips a dimension the variable does not carry, so
position and role coincide only because qualification currently requires all five sky axes. That is
another module's rule, decided for another reason in ADR 0009. A promise about positions would have
tied the two together without either saying so.

## Considered options

**Keep the constant and promise the order.** It is the smaller change, and it is what the README
said. It makes the order part of every profile's contract, for no reader that depends on it.

**Replace positional ranges with ranges named by role.** A `ReadRequest` keyed by role would make
the question disappear from requests altogether. It is a larger change to the one call every
consumer makes, for the same answer `AxisIndex` gives: the index is still needed to address the
destination, which is dense in the same logical order.

## Consequences

Removing `kXradioImageAxisOrder` breaks a consumer that used it. carta-backend-2 did not.

The README says axes are found by role, names the XRADIO order as the profile's rather than as a
promise, and its example sets its two spatial ranges through `AxisIndex` instead of writing five
ranges by position. `conformance_test` pins the profile's order directly, as a statement about what
the profile does, and checks that `AxisIndex` finds each role where it is and finds nothing for
`other`.

`kLogicalAxes` pairs each name with its role. The two lists it replaced matched by position -- the
names in `DescribeAxes` and the roles in the constant -- which was one more place where an index
meant something only because two arrays happened to agree.

## What this does not decide

**Whether the XRADIO order should change.** It does not, and nothing here gives a reason to: it is
still what every consumer of this profile has seen, and the conformance test would say so if it
moved.

**Whether a role can be played twice.** `AxisIndex` finds the first axis playing a role. Every role
but `other` is played by at most one axis of an image this library describes, and no image it
describes has an `other` axis today. If one ever has two, finding "the" axis for that role is a
question for whoever introduces them.
