# Control memory, cancellation and read-ahead

## Size the cache and temporary memory

`ContextOptions::cache_bytes` sizes the context's cache of decoded chunks,
which every operation shares unless it names a cache pool of its own. An absent
value keeps TensorStore's default; zero retains no decoded chunks. Size it from
the chunks your application will revisit, not from its output pixels: a chunk
is cached as it decodes, in the image's stored type, and the image's flag
takes cache space of its own. `Image::DecodedChunkBytes` gives a chunk's size.

`ReadOptions::read_budget_bytes` bounds the memory one read holds beyond the
caller's destination while it works, separately from what the cache retains.
Zero uses the library's default. It is not a process-memory cap: a read holds
at least one chunk however small the budget, reads running at once each have
their own, and a single plane read with no progress callback may hold the
whole plane (see [`ReadOptions`](../api/read.md)).

A scan that will not revisit its chunks can use a cache pool of its own that
holds nothing:

```cpp
auto pool = context.NewCachePool(0);
if (!pool) {
    // Report pool.error().
    return;
}
ReadOptions options;
options.control.cache_pool = pool.value();
```

This keeps the scan from evicting the chunks an interactive session reuses.

## Report progress and cancel

```cpp
ReadOptions options;
options.control.deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
options.control.cancellation_requested = [&stop] { return stop.load(); };
auto result = image.Read(request, {pixels.data(), pixels.size()}, options,
    [](std::size_t written, std::size_t total) {
        std::cout << written << "/" << total << '\n';
        return true;  // Return false to cancel.
    });
```

Here `stop` is a live `std::atomic<bool>`. Cancellation and deadlines are checked
at storage-operation boundaries; an in-flight decode is not forcibly interrupted.
On failure, including cancellation, discard the destination. During progress
callbacks, the reported prefix is complete.

`ReduceSpectral` and `ComputeHistogram` deliver results to a callback, their
sink, which cancels the operation by returning false. A cube histogram reports
progress through a separate callback that can also provide a snapshot of what
has been read so far. Producing that snapshot costs work; request one only when
needed. A small cube histogram can finish without emitting progress callbacks.

## Run operations at once

`Image::Read` and `Image::Prefetch` can run on one handle from several threads
at once. Concurrent reductions share the context's worker pool and interleave
their computation read by read.

## Read ahead during animation

`Image::Prefetch` decodes selected chunks into the cache without producing
pixels. `ReadAhead` uses it to decode the frames an animation will show next:

1. Create it with `ReadAhead::For({{image, options}})` and check the result.
2. After each displayed frame, call `Served` with the frame's start time,
   whether it was late, the shown planes and upcoming frames, nearest first.
3. Keep the controller alive while playing and call `Cancel` when stopping.

A **run** is the chunks one plane's read decodes: the plane's area times the
chunk depth. The controller decodes the next run on its own thread while the
frames of the current one are shown. Each cache an animated image reads through
must hold two runs of every image that reads through it, the one playing and
the one decoded ahead; otherwise creation reports `buffer_too_small`. A late
frame while a prefetch is running stops further read-ahead. Cancellation stops
scheduling promptly but still waits for storage work already under way. See the
[read-ahead API](../api/read-ahead.md) for precise controller lifetimes and policy.

## Complete read-control example

The scan uses a cache pool that holds nothing; the prefetch and the read after
it use the context's 64 MiB cache. The cancelled read starts with its stop flag
already set, so the demonstration does not depend on timing.

1. Creates a context with a 64 MiB cache and reads a bounded patch at channel 1
   (0 for a single-channel image) with progress, a deadline and the empty pool.
2. Sets the stop flag and reads again, checking for the `cancelled` error.
3. Prefetches the patch into the shared cache, then reads it.

Build and run it using the [example commands](../examples.md#build-and-run).
Its helpers are documented in [shared setup](../examples.md#shared-setup).

```{literalinclude} ../../examples/read_control.cc
:language: cpp
```

## Complete animation example

This program follows the steps [above](#read-ahead-during-animation) over a
cropped patch of at most four channels, then reports how many prefetches ran
and whether the controller stopped. It does not measure a display deadline: a
real viewer supplies its own frame schedule and lateness.

Its 64 MiB cache may not hold two runs of a large dataset; creation then fails
with the capacity it needs rather than silently reading nothing ahead. The
number of prefetches in this fast CLI is not a performance measurement.

```{literalinclude} ../../examples/animation.cc
:language: cpp
```
