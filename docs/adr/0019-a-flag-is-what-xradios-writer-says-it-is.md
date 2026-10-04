# A flag is what XRADIO's writer says it is: true where a pixel is bad, linked by data_groups

What a flag looks like on disk, what its values mean and which image it belongs to are XRADIO's to
say, and this library had all three from fixtures written by hand with zarr-python. Each of the
three was wrong against what XRADIO 1.2.3 writes when it converts a CASA image with an internal mask,
and each failed silently, the one way a pixel mask must not: the image opened, and its flagged pixels
counted as valid.

## What XRADIO writes

- **The values are true where the pixel is bad.** XRADIO's CASA reader takes casacore's mask, true
  for a good pixel, and stores its `logical_not`. A flag is a flag, not a mask of good pixels.
  `ApplyPixelMask` dropped the zeros, so it kept exactly the flagged pixels and dropped the rest:
  on a 5 × 4 plane with one masked pixel it read 1 valid pixel instead of 19.
- **The array is an `int8` with `dtype: "bool"` in its attributes.** That is how xarray encodes every
  boolean variable it writes to Zarr, and XRADIO writes through xarray. Only a Zarr `bool` was taken
  as a flag, so such an image had no pixel mask at all and its masked pixel was read as valid: 20
  instead of 19. Only `int8` carries the encoding, because only `int8` is what xarray writes.
- **The link is in the root's `data_groups`.** A group's `flag` is its `sky` image's flag, by XRADIO's
  schema, and the converted image names no flag of its own. ADR 0001 said flags were linked on the
  image's attributes; they are not. With two sky images and two matching flags, the guess from
  dimensions found two candidates for each and masked neither.

## How a flag is found

In this order, the first that answers deciding:

1. The image's own `flag` attribute, which a hand-written store may still use. It is the image's
   own statement and outranks a group's.
2. A `data_groups` entry whose `sky` is this image. Groups that share a sky and name the same flag
   are one declaration; naming two different ones is refused as invalid metadata, not chosen between.
3. A guess from dimensions, as before, leaving out every flag a data group declares for some other
   image. A flag XRADIO writes for a point spread function is not linked by any group, and this is
   what still finds it once the sky's flag is accounted for.

A declared flag, by either route, is binding: one that cannot serve closes the image rather than
opening it unmasked, as ADR 0001's rule for declared flags already said.

## Consequences

`tests/data/images/zarr/xradio/conformance_flagged` is converted by the pinned XRADIO from two masked
CASA images, and the conformance test reads and reduces both images against the pixels their masks
named. The hand-written fixtures were turned over to the same convention, so their good pixels are
the ones they were. The bench generator's flags were turned over too, and its `FORMAT_VERSION` raised
so that a dataset written with the old convention is regenerated rather than read inverted.

Nothing outside the library reads a flag's bytes: a consumer sees a flagged pixel as NaN.
