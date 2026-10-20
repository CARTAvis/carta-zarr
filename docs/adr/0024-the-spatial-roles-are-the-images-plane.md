# The spatial roles are the image's plane, the sky's or the aperture plane's

XRADIO image schema 0.0.2 has two planes. The sky plane carries the direction cosines `l` and `m`.
The aperture plane, its Fourier conjugate, carries `u` and `v` in wavelengths, each with the `crval`
and `cdelt` it was cut from. This library opened only the sky plane, and listed every aperture-plane
variable as not openable.

Some aperture-plane images are ordinary real images. CASA's `ia.fft` writes the amplitude and the
phase of a transform as float32 images on a linear `UU`/`VV` coordinate in `lambda`. CARTA's test
images `UVamp.image` and `UVphase.image` are two, and CARTA opens them through casacore like any
other image. XRADIO 1.2.5 and 1.3.0 convert them to an `APERTURE` typed `"aperture"`:

- the data group is `{"base": {"aperture": "APERTURE"}}`;
- the mask becomes `FLAG_APERTURE`, which no data group names;
- there is no `coordinate_system_info`.

The two versions write the same metadata, compression aside.

## Decision

**A real-valued variable on either plane is an image.** A variable is on the plane whose coordinates it
names among its dimensions: the sky's if it names `l` and `m`, the aperture plane's if it names `u` and
`v` and not the sky's. It is held to the same rules on either:

- five axes;
- extents that agree with its coordinates;
- a declared flag that can mask it.

A complex variable is listed and not openable on either plane, with
`unsupported_data_type`. CARTA opens no complex image, and a visibility grid is one.

**The spatial roles name the image's two plane axes.** `u` plays `spatial_x` and `v` plays `spatial_y`,
as `l` and `m` do. Reads, reductions, histograms and reading ahead reach the axes by role (ADR 0012),
so none of them knows which plane it walks. An image's axes are reported as its plane's spatial pair,
then frequency, polarization and time.

**The descriptor says which plane.** An image of the sky has a `direction`. An image of the aperture
plane has an `aperture` instead: an `ApertureCoordinate` with the unit `u` and `v` share and, for each
of them:

- the reference pixel, 1-based as the direction's;
- the reference value, which is XRADIO's `crval`;
- the increment.

There is never both, and no separate enumerator says which plane: the coordinate present is the answer.
The axis names are the store's (`u`, `v`). `UU` and `VV` are casacore's names, for a consumer
building a casacore coordinate to give.

**Fitted as a direction axis is.** A transform's grid is linear by construction. So the increment is
reported even when the samples come back unevenly spaced, a degenerate axis says so, and nothing is
scaled. The reference pixel is where the samples reach `crval`; `u` runs downward as XRADIO writes it.

**Each plane's coordinates are required of a dataset with an image of it.** The sky plane's coordinates
and `coordinate_system_info` are required only when an image of the sky opens. A dataset of
aperture-plane images alone needs `u`, `v` and the three both planes share. Its having no coordinate
system is not a fault: XRADIO writes none for it.

An aperture-plane variable written without its own `u` or `v` does not open, and is listed as
malformed. The sky images beside it still open. This rule is for hand-written stores, such as the
`minimal` fixture; XRADIO always writes the coordinates.

**The aperture plane has no beams.** XRADIO 1.2.5 and 1.3.0 name a `BEAM_FIT_PARAMS_APERTURE` on the
image and write no such array, so `ReadBeams` answers with none rather than reading the name.

**A flag is found as on the sky.** No data group names `FLAG_APERTURE`, so it is found by ADR 0019's
search by dimensions and shape.

## Consequences

- `ImageDescriptor::aperture` and `ApertureCoordinate` are new public API. This changes the ABI, and
  it is part of the untagged 0.3.0 (ADR 0020).
- `DiagnosticCode::unsupported_coordinate_plane` is no longer emitted. A complex variable on either
  plane is `unsupported_data_type`. The enumerator stays, so the codes after it keep their values.
- A consumer that took `spatial_x` and `spatial_y` to be the sky's must look at `direction` and
  `aperture` first. carta-backend builds a casacore `LinearCoordinate` named `UU`/`VV` in `lambda`
  from `aperture`, where it builds a `DirectionCoordinate` from `direction`.
- The conformance fixture `conformance_aperture` is the amplitude of the cube's transform, masked,
  converted by XRADIO 1.2.5. The test opens it, reads it masked and reduces it by role.
