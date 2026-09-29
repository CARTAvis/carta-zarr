# A handle always refers to something, and moving one copies it

`Context`, `Dataset` and `Image` are shared handles: each holds a `shared_ptr` to an implementation,
and copying one shares it. The only way to get one is from the factory that made what it refers to
-- `Context::Create`, `Dataset::Open`, `Dataset::OpenImage` -- and none of them hands out an empty
one.

Moving one did. The move operations were defaulted, and a defaulted move of a `shared_ptr` leaves
null behind, so every handle had an empty state reachable by `std::move` and by nothing else. No
consumer could construct one, and the public interface gave no way to ask whether a handle was in
it. The library still had to answer for it: every `Image` entry point opened by refusing an empty
handle with `invalid_argument`, `Dataset::Open` refused an empty context, and `descriptor()` and
`chunk_geometry()` returned a static empty object rather than dereference null.

Each handle now moves by delegating to its copy, so the handle left behind still refers to what it
did, and every one of those checks is gone.

## The empty state was the library's, not the consumer's

A state a consumer can reach only by moving from a handle and then using it anyway is a state the
library created and then defended against. The defence was not free, and it was not uniform: before
the checks were gathered into `WithReadableImage`, two entry points refused an empty handle and
three also refused a null store that `Image::Impl`'s constructor cannot produce, so which of the two
a reader was looking at took working out.

A consumer that genuinely needs a handle it fills in later already has the standard answer.
carta-backend-2 held its image in a `std::optional<carta::zarr::Image>` before this change and holds
it the same way after; the empty case is then the optional's, where the compiler and every reader
can see it.

## Considered options

**A default-constructible handle with `valid()` or `operator bool`.** This makes the empty state
first class: a consumer can construct it and test for it. It also keeps it, so every entry point
keeps refusing it and every consumer has one more thing to check before each call -- the cost the
checks were already paying, now paid on both sides.

**Leave the moves defaulted and keep the checks.** Nothing is wrong with the checks as code. They
guard a state that exists only because of how the move was declared.

**Move by copying.** Chosen. The cost is that a move is one atomic reference-count increment rather
than a pointer taken, beside calls that read metadata or decompress chunks. The moves stay
`noexcept`, which is what lets a `std::vector` of handles relocate by moving; a test asserts it.

## Consequences

Calling through a handle that has been moved from used to report `invalid_argument` and now works.
carta-backend-2 never did it.

`schema_probe_test` moves each of the three handles and uses the one left behind, down to reading a
pixel. With the defaulted moves restored the test binary crashes, since nothing checks for null any
more; that is the regression it is there to catch.

carta-backend-2's `GetZarrContext` returned a `shared_ptr<const Context>`, two layers of sharing
around one TensorStore context. It returns the `Context` by value now, and keeps the process-wide
one in a `std::optional` that is empty only before the first call.

The README states the rule among what the API promises, and `carta_zarr.h` says it above the three
classes.

## What this does not decide

**Concurrent assignment to one handle.** A handle is a value. Reading through one handle from
several threads is safe, as the README says; assigning to the same handle object from one thread
while another reads it is a race, exactly as it is for the `shared_ptr` inside it.
