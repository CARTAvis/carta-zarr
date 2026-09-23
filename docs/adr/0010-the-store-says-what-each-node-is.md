# The store says what each node is, and its callers decide what to do about it

ADR 0004 put the seam at the transport and gave `Store` everything that turns bytes into meaning:
"node path validation, JSON parsing, consolidated metadata, array metadata, storage layout, and the
caches". The listing was the one place that rule was not followed. `Store::ListNodes` read every
node's document and handed back only the names, and each caller went back to the store to work out
what a node was -- three of them, three ways. The declared size read `node_type` out of the
document. Qualification parsed every node as an array and looked again when the parse failed, to
see whether it had been a group. Flag selection parsed every node, groups included, and ignored the
answer. They agreed only where those three happened to coincide.

`Store::Inventory` replaces it. Every node comes back under one canonical name, sorted, as a
`group`, an `array`, or an `unrecognised` node, and an array comes back with its metadata already
parsed -- a pointer to the one the store holds, so nothing is copied and a group is never parsed as
an array. CONTEXT.md calls it the node inventory.

## Why the parse is carried, not just the kind

Carrying the kind alone would have settled "what is this node" and left "did its metadata parse" for
each caller to find out by asking again. That second question is the one the callers actually
disagree about, so leaving it out would have kept the reconstruction exactly where the divergence
was. With the parse carried, the three states -- a group, an array that parsed, a node that claimed
to be an array and did not -- are in the type, and what is left to each caller is policy.

It costs nothing on the path a consumer takes. Opening a dataset runs discovery, which already
parsed every node, and every other caller walked every array too. The one path that pays is asking a
fresh store about a single image before anything has listed it, which only tests do.

## The policies stay different, on purpose

What a caller does with a node that is not a usable array is its own business:

- **Image discovery diagnoses it**, because it may have been an image, and a store whose sky
  variable disappears without a word reads as a store nobody recognised.
- **Flag selection passes over it**, because a mask candidate that cannot be read is simply not the
  mask. Complaining would make a guess into a finding about the store.
- **The declared size refuses an array it cannot size**, because a total that leaves out an array it
  knows is there is wrong rather than small -- and counts neither a group nor an unrecognised node,
  because neither is an array.

These were recoverable only from three control flows, and nothing said they were meant to differ.
They are written down beside `Inventory` now. Unifying them was considered and rejected: any one
policy applied to all three would be wrong for two.

## One node, one name

`Store` had two rules for what a node is called. `NormalizeNodeName` drops a `.` component, refuses
`..` and absolute paths, and yields the key every read is cached under. `NormalizeMetadataKey`, which
filed the root's consolidated copy, only stripped a leading slash and a trailing `zarr.json`. A block
that spelled a node `./SKY` was filed under `./SKY`, where a read of `SKY` never looked: the node was
read again, cached twice, and listed as `./SKY`.

That became a defect rather than a waste when ADR 0009 made the listing decide what opens. An image
listed as `./SKY` could not be opened as `SKY`. A consolidated key is now stripped of what makes it
name a document, then held to the same rule as every other node name; so are the names a transport
lists.

### Two refusals, and why they are refusals

A consolidated block that lists a key no node path can be, or lists one node under two spellings,
refuses the store. Passing over the offending key looks gentler and is not. For a consolidated store
the block *is* the listing -- nothing enumerates the transport behind it -- so a key dropped here is a
variable missing from the dataset without a word said, and since ADR 0009 it is also a variable that
cannot be opened. A store refused with a reason is the better answer.

ADR 0004 calls consolidated metadata "a copy that saves those reads, not a substitute". That remains
true of the documents. It was never true of the listing, and this is where the difference shows.

### Canonical names are for the store's own references

`FindNode` looks a name up as given and does not respell it, so `OpenImage("./SKY")` answers
`not_found`, while an image declaring `flag: "./MASK_0"` finds its mask. That asymmetry is deliberate.
A flag attribute, a coordinate name, a consolidated key -- those are spellings the store wrote for
itself, and the library has to recognise them. An image id is a name the library handed out, and a
consumer passing it back got it from the listing. Accepting any spelling there would put a
non-canonical id into a descriptor and two cache entries behind one image, to honour a promise nobody
asked for.

## A node whose metadata will not parse is diagnosed

CONTEXT.md already decided this. Under **Diagnostic** it lists "a node skipped because its metadata
would not parse" among the things diagnosed rather than refused, because the rest of the dataset is
still readable. The listing did the opposite: it gave up on the first document that would not parse,
so one stray broken `zarr.json` beside a good image closed the dataset. It was also inconsistent
with its neighbour -- an array whose metadata would not parse as an array was already a diagnostic.

A document that will not parse is an `unrecognised` node now, carrying the parser's error. Only that
failure: a read that failed outright -- a transport error, a node listed and then gone -- says the
hierarchy could not be taken, and still refuses. The line is between "this node could not say what it
is" and "the store could not be read".

Qualification diagnoses it as `unrecognised_node`, not `unreadable_array`. The two say different
things: a node that claimed to be an array and whose metadata will not parse as one, and a node that
did not say what it was, where calling it an array is a guess. A node whose `node_type` names
something Zarr does not define is the second kind as well.

The declared size leaves it out. It might have been an array, so the total may be smaller than it
should be; ADR 0008 already says a declared size is no bound in either direction, and refusing would
leave a dataset that opens with no size at all.

## Consequences

ADR 0009 left "the other readings of a parse failure" undecided, naming `DetermineFlag` and
`TotalArraySizeBytes`. That is settled here. 0009 is left as written, the way 0009 left 0007.

A node whose `node_type` is present and not a string used to throw. Two of the three readers asked
for it as a string with a default, which nlohmann refuses with `type_error.302`, so discovery threw
and the dataset did not open. It is an unrecognised node now, like any other.

carta-backend-2 sees one change: a dataset holding a stray document that will not parse is listed
with its images and a diagnostic, where the file browser used to report an error. Nothing there
matches on a diagnostic code, so `unrecognised_node` is additive.

The store's lock-order DAG gains `inventory → array_metadata`. It stays acyclic, and the comment that
says so in `store.h` changed with it.

## What this does not decide

**The transport still reads `node_type`.** `FilesystemTransport` parses each node's document while
listing, to know whether to descend into it rather than walking millions of chunk files.
`transport.h` says so and why. It is a fourth reader of the field, below the seam, deciding where to
walk and never what anything means; it is not one of the readings this unifies, and should not be
"fixed" into the inventory.

**By-name reads still reach past the inventory.** A declared flag, a beam table, and a coordinate
are read by the name the store wrote, through `ReadArrayMetadata`, which does not consult the
inventory. The snapshot ADR 0009 relies on is enforced for images, which are what a consumer names.
For the store's own references it holds only because the first read of a name is remembered for the
life of the store.

**Diagnostics dropped on `no_match`.** Unchanged from ADR 0009: a store of nothing but unsupported
variables is still reported as unrecognised with no reason attached.
