# Images are identified by variable name and discovered by inspection, not by data group

*Amended by [ADR 0023](0023-only-xradio-image-schema-0-0-2-is-read.md): only schema 0.0.2 is read,
and a dataset naming it must have `data_groups`. Which variables are images is still decided by
inspecting every variable, as below.*

*Amended for 0.3.0 ([#59](https://github.com/CARTAvis/carta-zarr/issues/59)): a variable is an image by
its `type`, which schema 0.0.2 gives every image, rather than by its dimensions alone; and the data
groups' own images are listed first. See Consequences.*

An XRADIO image dataset declares its images twice: implicitly, as data variables carrying the `l`
and `m` coordinates, and explicitly, in the optional root attribute `data_groups`, which maps roles
such as `sky` and `residual` to variable names. We identify an image by its data-variable name
alone, and we discover the set of images by inspecting every variable in the dataset rather than by
reading `data_groups`. (This was written before stores were versioned. Schema 0.0.2 requires
`data_groups` and types every image, and the amendments above say what that changed: an image is
now what its type says, and inspection still finds every one.)

## Considered options

Identifying an image by a `(data group, role)` pair was the alternative, and it is what the XRADIO
schema is built around. We rejected it because CARTA presents images in a flat dropdown, the way it
presents FITS HDUs, and because a variable shared by two data groups would appear twice under that
scheme. Variable names are unique within a dataset, so the name alone is sufficient to address one.

Discovering images from `data_groups` was rejected on evidence rather than on principle: the
XRADIO implementation writes `data_groups` (defaulting to `{"base": {}}`) but populates it sparsely,
commonly with no more than `{"base": {"sky": "SKY"}}`, so a dataset containing a `RESIDUAL` would
report only its `SKY`. Inspection reads what is actually in the store.

## Consequences

An image is a data variable whose `type` is one of schema 0.0.2's images -- `sky`,
`point_spread_function`, `primary_beam` and `mask` on the sky plane; `aperture`, `visibility` and
`uv_sampling` on the aperture plane -- and that carries the whole axis set of that plane (`time`,
`frequency`, `polarization`, and `l` and `m`, or `u` and `v`). Until 0.3.0 the type was not asked:
an image was any real variable over the whole axis set whose `type` was not `"flag"`. That was
chosen before stores were versioned; schema 0.0.2 types every image, and a rule from dimensions alone
offered whatever had the right axes. The whole axis set is still required of an image, and one typed
as an image without it is malformed rather than passed over. Classification is still never driven by
variable names:

- `right_ascension` and `declination`, which XRADIO writes as float64 arrays over `(l, m)` with no
  `type`, are not images, nor are the normalization variables and the beam fit parameters.
- `MASK_0`, despite its name, is not an image. The XRADIO CASA reader names variables after CASA's
  internal masks but tags them `type: "flag"`, so they are pixel masks.
- A variable a data group names in an image role is held to that role's type: a `sky` typed
  `"Intensity"` is the store contradicting itself, and is reported as malformed.

One name is read: an image that is not a sky image, whose flag XRADIO links to it by nothing else,
takes the flag XRADIO names `FLAG_<image>` (ADR 0019). That finds a flag for an image already
classified; it never makes a variable an image or a flag.

The listing puts the data groups' own images first: the base group's, which is the store's main image,
then every other group's by group name (a group's own image is its `sky`, or its `aperture` when it has
no sky). The rest follow, ordered by a known variable-name list so that XRADIO's names come first,
then by name, and they open as before. The first openable image is the default. The name list orders;
it never decides membership.

`data_groups` is carried into `ImageDescriptor` as descriptive metadata. It is not used to locate an
image's beam variable, which XRADIO names on the image variable's own attributes. It is used to locate
an image's flag, which XRADIO's writer names nowhere else; this paragraph used to say otherwise, from
fixtures written by hand rather than by XRADIO. See ADR 0019.
