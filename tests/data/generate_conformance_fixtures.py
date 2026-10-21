#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = [
#   "xradio[zarr,casacore]==1.3.0",
#   "casatools==6.7.6.14",
#   "casaconfig",
#   "astropy==7.1.0",
#   "numpy==2.3.1",
# ]
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate the conformance fixture by converting a FITS image with XRADIO itself.

`generate_zarr_fixtures.py` builds stores with zarr-python directly, which pins our *belief* about
what XRADIO writes. This script pins XRADIO instead: it makes a FITS image with astropy and hands it
to `xradio.image.open_image()` / `write_image()`, so the result is whatever the released XRADIO
actually produces. Regenerate it when the pinned XRADIO version is bumped, and let the conformance
test report which of our assumptions the new version broke.

XRADIO 1.3.0 writes image schema 0.0.2, the only schema carta-zarr reads. It is the release after
1.2.5, and it writes the same fixtures apart from the compressor.

XRADIO's FITS reader needs a CTYPE on every axis and refuses an axis type it does not support. The
header written here also carries DATE-OBS; without it or MJD-OBS, the time is MJD 0.

The flagged and aperture fixtures are converted from CASA images, which are written with casatools
and read by XRADIO through arcae, its casacore extra.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

OUTPUT_DIR = Path(__file__).parent / "images" / "zarr" / "xradio"
FIXTURE_NAME = "conformance"
FLAGGED_FIXTURE_NAME = "conformance_flagged"
# One channel, as a continuum image has, and an image of something other than intensity: XRADIO 1.3.0
# writes the channel's width and the BTYPE, and these are the only two things carta-zarr reads that
# the cube above cannot show.
CONTINUUM_FIXTURE_NAME = "conformance_continuum"
CONTINUUM_BTYPE = "Spectral index"
# The same cube with a pointing center apart from its reference direction, and FITS cards XRADIO
# carries along uninterpreted, one of each kind of value it keeps: XRADIO 1.3.0 writes the first as
# the image's `pointing_center` and the rest, lower-cased, as its `user`. The cube above has neither,
# so its pointing center is its reference direction and its `user` is empty.
KEYWORDS_FIXTURE_NAME = "conformance_keywords"
KEYWORD_CARDS = {
    "OBSRA": 45.5,
    "OBSDEC": 30.25,
    "INSTRUME": "Test instrument",
    "DISTANCE": 101.5,
    "NCOMBINE": 3,
    "PBCOR": True,
}

# The amplitude of the cube's Fourier transform, a real image of the aperture plane as CASA's
# `ia.fft` makes one (UVamp in CARTA's test images is another), with a mask over a few known pixels.
# XRADIO 1.3.0 writes it as APERTURE on (u, v) in wavelengths, its mask as FLAG_APERTURE that no data
# group names, and no coordinate_system_info.
APERTURE_FIXTURE_NAME = "conformance_aperture"
APERTURE_MASKED_PIXELS = [(0, 0, 0, 0), (2, 1, 1, 1)]

# The pixels each CASA image's mask marks bad, as (x, y, stokes, frequency) -- casacore's axis order.
# casacore's mask is true for a good pixel and XRADIO's flag is true for a bad one, so these are the
# pixels the converted flags mark true. Every value of the FITS image is finite and unique, so a
# reader that ignores a flag, inverts it, or applies another image's flag reports a different count
# of valid pixels on a known plane and keeps a value it should have dropped.
MASKED_PIXELS = {
    "sky": [(1, 2, 0, 0), (3, 0, 1, 1)],
    "sky_other": [(4, 3, 2, 0)],
}

AXES = [
    # (CTYPE, CRVAL, CDELT, CUNIT)
    ("RA---SIN", 45.0, -0.001, "deg"),
    ("DEC--SIN", 30.0, 0.001, "deg"),
    ("STOKES", 1.0, 1.0, ""),
    ("FREQ", 1.4e9, 1.0e6, "Hz"),
]


def write_fits(path: Path, channels: int = 2, btype: str = "Intensity", cards: dict | None = None) -> None:
    # FITS axis order is the reverse of the numpy shape: NAXIS1 is the fastest-varying axis.
    values = np.arange(channels * 3 * 4 * 5, dtype=np.float32).reshape(channels, 3, 4, 5)
    hdu = fits.PrimaryHDU(values)
    header = hdu.header
    for index, (ctype, crval, cdelt, cunit) in enumerate(AXES, start=1):
        header[f"CTYPE{index}"] = ctype
        header[f"CRVAL{index}"] = crval
        header[f"CDELT{index}"] = cdelt
        header[f"CRPIX{index}"] = 1.0
        header[f"CUNIT{index}"] = cunit

    header["BUNIT"] = "Jy/beam"
    header["BTYPE"] = btype
    header["OBJECT"] = "Zarr conformance source"
    header["OBSERVER"] = "CARTA"
    header["TELESCOP"] = "Test scope"
    header["RADESYS"] = "FK5"
    header["EQUINOX"] = 2000.0
    header["SPECSYS"] = "LSRK"
    header["RESTFRQ"] = 1.420405751e9
    header["DATE-OBS"] = "2020-05-31T12:00:00.000"
    header["TIMESYS"] = "UTC"
    header["LONPOLE"] = 180.0
    header["LATPOLE"] = 30.0
    header["PC1_1"], header["PC1_2"] = 1.0, 0.0
    header["PC2_1"], header["PC2_2"] = 0.0, 1.0
    header["BMAJ"], header["BMIN"], header["BPA"] = 2.0e-5, 1.0e-5, 0.1
    for keyword, value in (cards or {}).items():
        header[keyword] = value

    hdu.writeto(path, overwrite=True)


def write_masked_casa_image(path: Path, fits_path: Path, masked: list[tuple[int, int, int, int]]) -> None:
    from casatools import image as image_tool

    image = image_tool()
    try:
        image.fromfits(outfile=str(path), infile=str(fits_path), overwrite=True)
        good = np.ones(image.shape(), dtype=bool)
        for pixel in masked:
            good[pixel] = False
        image.putregion(pixelmask=good)
    finally:
        image.done()


def write_flagged(fits_path: Path, scratch: Path) -> None:
    """Two CASA images with internal masks, converted together.

    This is the shape of what XRADIO writes from CASA images: each sky image's flag a Zarr bool, true
    where the pixel is bad, and tied to its image only by the root's `data_groups` -- the image itself
    names no flag. Two of them, of the same shape, so each image's flag can be told only from what its
    data group declares.
    """
    from xradio.image import open_image, write_image

    fixture = OUTPUT_DIR / FLAGGED_FIXTURE_NAME
    shutil.rmtree(fixture, ignore_errors=True)
    stores = {}
    for role, masked in MASKED_PIXELS.items():
        casa_path = scratch / f"{role}.im"
        write_masked_casa_image(casa_path, fits_path, masked)
        stores[role] = str(casa_path)
    xds = open_image(stores)
    write_store(write_image, xds, fixture)


def write_aperture(fits_path: Path, scratch: Path) -> None:
    """The amplitude of the cube's transform over its two direction axes, masked, and converted."""
    from casatools import image as image_tool
    from xradio.image import open_image, write_image

    sky_path = scratch / "aperture_sky.im"
    amplitude_path = scratch / "aperture_amp.im"
    image = image_tool()
    try:
        image.fromfits(outfile=str(sky_path), infile=str(fits_path), overwrite=True)
        image.fft(amp=str(amplitude_path), axes=[0, 1])
    finally:
        image.done()
    amplitude = image_tool()
    try:
        amplitude.open(str(amplitude_path))
        good = np.ones(amplitude.shape(), dtype=bool)
        for pixel in APERTURE_MASKED_PIXELS:
            good[pixel] = False
        amplitude.putregion(pixelmask=good)
    finally:
        amplitude.done()
    write_store(write_image, open_image(str(amplitude_path)), OUTPUT_DIR / APERTURE_FIXTURE_NAME)


def write_store(write_image, xds, fixture: Path) -> None:
    """Writes an image store at `fixture`. XRADIO names an image store `<name>.img.zarr`, which the
    fixture's path is not, so it is written under that name and moved."""
    written = Path(write_image(xds, str(fixture), out_format="zarr", overwrite=True)[0])
    if written != fixture:
        shutil.rmtree(fixture, ignore_errors=True)
        written.rename(fixture)


def write_keywords(scratch: Path) -> None:
    from xradio.image import open_image, write_image

    keywords_path = scratch / f"{KEYWORDS_FIXTURE_NAME}.fits"
    write_fits(keywords_path, cards=KEYWORD_CARDS)
    write_store(write_image, open_image(str(keywords_path)), OUTPUT_DIR / KEYWORDS_FIXTURE_NAME)


def main() -> None:
    from xradio.image import open_image, write_image

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fixture = OUTPUT_DIR / FIXTURE_NAME
    shutil.rmtree(fixture, ignore_errors=True)

    scratch = Path(tempfile.mkdtemp(prefix="carta-zarr-conformance-"))
    fits_path = scratch / f"{FIXTURE_NAME}.fits"
    write_fits(fits_path)
    try:
        xds = open_image(str(fits_path))
        write_store(write_image, xds, fixture)
        write_flagged(fits_path, scratch)

        continuum_path = scratch / f"{CONTINUUM_FIXTURE_NAME}.fits"
        write_fits(continuum_path, channels=1, btype=CONTINUUM_BTYPE)
        write_store(write_image, open_image(str(continuum_path)), OUTPUT_DIR / CONTINUUM_FIXTURE_NAME)

        write_keywords(scratch)
        write_aperture(fits_path, scratch)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()
