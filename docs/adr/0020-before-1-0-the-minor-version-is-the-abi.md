# Before 1.0 the minor version is the ABI

A release of this library says which other releases can stand in for it in three places: the soname a
binary records when it links, the version node every exported symbol carries on ELF, and the package
version file `find_package(CartaZarr <version>)` asks. Until 1.0 all three name the major and the
minor -- `libcarta-zarr.so.0.1`, `CARTA_ZARR_0_1`, `SameMinorVersion` -- and from 1.0 on they name the
major alone. `CARTA_ZARR_ABI_VERSION` in CMakeLists.txt is the one place that decides it, and the
`carta-zarr-version-compatibility` test holds the version file and the soname to it.

## Why the minor

Before this, all three named the major, and the major was 0. So 0.2 kept the soname 0.1 had and
satisfied `find_package(CartaZarr 0.1)`, although 0.x has already changed its ABI without warning:
`SpectralTotals` gained a field in ADR 0018 and the version did not move. A backend built against 0.1
and handed 0.2 by its package manager would then load it, read a struct of the wrong size, and say
nothing.

A library before 1.0 is expected to change its interface while its first consumers find out what it
should have been, and this one will: its consumer's changes go upstream as separate pull requests,
and their review is where the interface is most likely to move. That leaves two honest promises.
One is 1.0 now, with every such change a new major. The other is 0.x with the minor as the
compatibility boundary, the way libraries in this position usually say it. The second spends
nothing on version numbers that mean "the review asked for a rename".

Within one minor a patch release keeps the ABI, as it will within one major after 1.0.

## What it costs

A consumer asking for `0` alone asks for 0.0 and is refused 0.1, because before 1.0 those are
different libraries. A consumer writes the minor it was built against, as carta-backend's
`find_package(CartaZarr 0.1 ...)` already does.

Raising the minor is now a decision with consequences downstream: every binary linked against the old
one has to be rebuilt. That is the point of it, and it is the reason a change to a public struct or a
public signature raises the minor in the same commit that makes it.

## Considered

**1.0 now.** Rejected for the reason above. Nothing has been released, the protocol and schema
changes the consumer needs are not upstream yet, and their review is likely to ask for changes to
this interface.

**Keep the major-only soname and write the rule in the README.** That is a promise the package
version file contradicts: `find_package` would still hand 0.2 to a consumer that asked for 0.1, and
nothing would check the README.
