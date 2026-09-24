# A caller's stop is one error, whichever way the caller said it

`ErrorCode::cancelled` answers four different things: a sink that returned false, a read's progress
callback that returned false, a `ReadControl::cancellation_requested` that returned true, and a
`ReadControl::deadline` that passed. We keep them one code. It looks like a conflation worth pulling
apart -- an architecture review proposed reporting a sink's stop as a successful outcome that says
where it stopped, keeping `cancelled` for the rest -- and this records why that was not done.

## Why the stop is not the library's to explain

A sink or a callback returning false tells the library to stop; it does not tell it why, and the why
is the only thing a consumer would want the distinction for. carta-backend's resumable region profile
shows it plainly. Its sink returns false in two cases: at a block boundary once its time slice is
spent, which is a pause it will resume from, and when its own partial callback says the region has
moved, which is an end. Both arrive here as the same false. The flag that tells them apart lives in
the backend because the reason does, and no outcome this library could return would let it go.

Nor would saying where the walk stopped add anything. A reduction hands over its blocks as they are
finished, and a consumer that resumes already counts the channels those blocks covered; that count
is what it resumes from.

So splitting the code changes all four entry points that can be stopped -- `Image::Read`,
`ReduceSpectral`, `ComputeHistogram` and `ComputeCubeHistogram` -- and every consumer's handling of
their results, and removes nothing from either side. The deletion test says the same: a second code
would be one more case to handle, not a complexity that stops existing.

## The deadline is the wrinkle

A deadline that passes is reported as `cancelled` too, and that one is not quite the caller's own
decision. The caller set the deadline but did not see the moment it passed, whereas a callback that
returned false is a line of the caller's own code. carta-backend reads `cancelled` as "the caller
stopped this" and logs nothing for it, which is right for the other three and would hide an overrun
here.

Nothing in carta-backend sets a deadline today, so nothing reaches this. A consumer that starts to
can tell the two apart by reading the clock itself when `cancelled` comes back. If that proves
awkward in practice, a distinct code for an expired deadline is the change to make -- not a split of
the other three.
