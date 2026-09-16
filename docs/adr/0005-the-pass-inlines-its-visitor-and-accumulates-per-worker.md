# The pass inlines its visitor, and accumulates once per worker

A **pass** is one ordered visit to every chunk an image read covers, shared by every reduction that
wants those pixels. Two things about its shape look like accidents of implementation and are not,
so they are recorded here: the visitor is a template parameter rather than a virtual interface or a
`std::function`, and the accumulators it hands out are keyed by worker rather than by task.

Both are the kind of thing a later reading would tidy up. A pass that took a `Visitor&` interface
would be easier to describe, and an accumulator keyed by task would be easier to reason about than
one keyed by a worker index the pool happens to expose. Each of those changes costs a measured
amount of time, and neither cost is visible in a Debug build or in any test this repository runs.

## The visitor is a template parameter

The per-pixel loop inlines through the visitor. `AccumulateRow` is instantiated four ways on
`(unit stride, masked)`, chosen once per row, and the histogram's `Add` and its row lambdas reach
the same state by being local types behind a template parameter. Type-erasing the visitor puts an
indirect call at the slab boundary and stops the enclosing loop nest from inlining.

`54731c1`, which made the accumulation loop vectorise, is worth about 25% of a whole-region
reduction: a 240 x 240 x 250 profile settles at 24-28 ms against 33 ms, with the reads unchanged at
10.7 ms. Its own message records what makes this dangerous: "That figure is from a Release build.
The same comparison in the default Debug build dir shows no difference at all, which is how this
change was nearly discarded as useless."

`WorkPool::Run` does take a `std::function`, and that is not in tension with this: it is paid once
per task, outside the pixel loop, not once per slab inside it.

## Accumulators are keyed by worker, not by task

A provisional histogram is half a megabyte of scattered writes. Keyed by task it is dragged from one
core's cache to another's once per plane. Measured on a 512x512x7776 ASKAP cube, warm, against 8.8 s
for not splitting at all: four workers 4.9 s, eight 8.7 s, twenty-eight 16.5 s. Keying by task
measured slower than not splitting at all, which is why the cap is the memory budget divided by what
one accumulator costs, and why `Accumulator` is `alignas(64)` — two accumulators sharing a cache
line trade it between cores once per pixel.

This is why the pass hands the visitor a worker index rather than a task index, and why
`WorkPool::Run` passes both.

## Consequences

The pass cannot be given a non-template entry point for convenience, and a second overload taking a
`std::function` would be a trap rather than a shortcut: it would compile, pass every test, and cost
a quarter of the reduction. If one is ever wanted for a cold path, it belongs behind a name that
says so.

The numbers above cannot be reproduced by anything in this repository. Every fixture is far too
small to show either effect, no test asserts a duration, and the default `build/` tree is Debug, so
the loop change measures as nothing there. A timing harness run by hand against `build-release/` is
what stands in for a regression test, and it is the only thing that does.

Neither decision reaches the public interface. `src/reduce/` is entirely
`carta::zarr::internal`, and only `include/carta-zarr/` is installed.
