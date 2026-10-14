# Memory, cancellation and animation

## Cache and temporary memory

`ContextOptions::cache_bytes` controls the shared decoded-chunk cache. An absent
value keeps TensorStore's default; zero retains no decoded chunks. Choose a
cache from the chunks your application will revisit, not just its output pixels:
chunks are cached in their stored type, with flag chunks counted separately.

`ReadOptions::read_budget_bytes` controls temporary work beyond the caller's
destination, separately from retained cache capacity. Zero uses the library's
default budget. It is not a strict process-memory cap: at least one chunk must
fit logically, concurrent reads each have a budget, and single-plane reads
without progress have a documented whole-plane exception.

Use `Image::DecodedChunkBytes` to estimate decoded chunk storage. For a scan
that will not revisit chunks, create a separate zero-byte pool:

```cpp
auto pool = context.NewCachePool(0);
if (!pool) {
    // Report pool.error().
    return;
}
ReadOptions options;
options.control.cache_pool = pool.value();
```

This keeps the scan from evicting an interactive session's working set.

## Progress and cooperative cancellation

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

Spectral and plane-histogram sinks cancel by returning false. Cube histogram
progress has a separate callback and can provide a snapshot of what has been
read so far. Producing that snapshot costs work; request one only when needed.
A small cube histogram can finish without emitting progress callbacks.

## Animation read-ahead

`Image::Prefetch` decodes selected chunks into the selected cache without
producing a pixel destination. `ReadAhead` adds an animation policy:

1. Create it with `ReadAhead::For({{image, options}})` and check the result.
2. After each displayed frame, call `Served` with the frame's start time,
   whether it was late, the shown planes and upcoming frames, nearest first.
3. Keep the controller alive while playing and call `Cancel` when stopping.

The controller prefetches the next distinct run of chunks on its own thread.
The caches must hold at least two runs; creation reports `buffer_too_small`
otherwise. A late frame while a prefetch is running disables further read-ahead.
Cancellation stops scheduling promptly but still respects in-flight storage work.

Read operations can share a handle concurrently. Concurrent reductions share
one worker pool and interleave computation read by read. See the
[read-ahead API](../api/read-ahead.md) for precise controller lifetimes and policy.
