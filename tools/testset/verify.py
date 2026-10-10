#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "astropy==7.1.0",
#   "numpy==2.3.1",
#   "zarr==3.2.1",
# ]
#
# [tool.uv]
# python-preference = "only-managed"
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Check that a FITS cube and a Zarr made from it hold the same pixels at the same coordinates.

A comparison between FITS and Zarr is only a comparison of the formats if both hold the same pixels.
What the Zarr means is taken from carta-zarr itself -- `carta-zarr-bench probe --describe`, which
opens it as carta-backend would -- rather than from its metadata read a second way here, so a Zarr
carta-zarr would refuse, or read differently, is refused. That description is held to the FITS
header: the image opens, as SKY, float32, in the FITS unit, in the chunks (and shards) asked for, with
a pixel mask or without one as --flag says; its direction coordinate has the FITS file's projection,
frame, equinox, reference pixel and value, increment, matrix, parameters and native pole; every
channel's frequency, the spectral unit, frame and rest frequency, the Stokes parameters, the time and
the observation date match; and so does the restoring beam of every plane. Then every pixel is
compared, bit for bit, a block of channels at a time, and the pixel mask is checked true exactly where
the FITS cube is NaN, in the same chunks as the pixels.

    verify.py CUBE.fits CUBE.zarr --flag {yes,no} --chunks L,M,F [--shards L,M,F]
              [--bench PATH] [--block-mib 1024]

--bench is the carta-zarr-bench executable, $CARTA_ZARR_BENCH by default. Prints what it compared, and
exits nonzero if anything differs.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import warnings

import numpy as np
import zarr
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

STOKES = {1: "I", 2: "Q", 3: "U", 4: "V"}


def first_difference(a: np.ndarray, b: np.ndarray) -> tuple[int, ...] | None:
    """Where `a` and `b` first differ, without listing everywhere they do."""
    differ = a != b
    if not differ.any():
        return None
    return tuple(int(i) for i in np.unravel_index(int(np.argmax(differ)), differ.shape))


def close(a, b, tolerance: float = 1e-9) -> bool:
    """Equal to a part in 10^9, or within 10^-9 of zero; a NaN or a missing value is never close."""
    try:
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    except (TypeError, ValueError):
        return False
    return a.shape == b.shape and bool(np.all(np.isfinite(a)) and np.allclose(a, b, rtol=tolerance, atol=tolerance))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fits")
    parser.add_argument("zarr")
    parser.add_argument("--flag", choices=("yes", "no"), required=True, help="whether the Zarr must carry a flag")
    parser.add_argument("--chunks", required=True, help="the chunk asked for, as l,m,frequency")
    parser.add_argument("--shards", help="the shard asked for, as l,m,frequency; none means unsharded")
    parser.add_argument("--bench", default=os.environ.get("CARTA_ZARR_BENCH"), help="the carta-zarr-bench executable")
    parser.add_argument("--block-mib", type=int, default=1024, help="how much of each to compare at a time")
    args = parser.parse_args()
    if not args.bench:
        parser.error("--bench or $CARTA_ZARR_BENCH must name carta-zarr-bench")

    failures = []

    def check(what: str, ok: bool) -> None:
        print(f"{what}: {'same' if ok else 'DIFFERENT'}", flush=True)
        if not ok:
            failures.append(what)

    def finish() -> int:
        print("OK" if not failures else f"MISMATCH: {', '.join(failures)}")
        return 1 if failures else 0

    # What carta-zarr sees. A dataset it does not open, or whose default image is not SKY, is refused
    # before anything else is compared.
    probe = subprocess.run([args.bench, "probe", args.zarr, "--describe"], capture_output=True, text=True)
    try:
        report = json.loads(probe.stdout)
    except json.JSONDecodeError:
        report = {"ok": False, "error": (probe.stderr or probe.stdout).strip()}
    image = report.get("image", {})
    check("carta-zarr opens it" + ("" if report.get("ok") else f" ({report.get('error', 'no reason given')})"),
          bool(report.get("ok")) and image.get("id") == "SKY" and report.get("default_image_id") == "SKY")
    if failures:
        return finish()
    described = image["description"]
    axes = {axis["name"]: axis for axis in image["axes"]}

    with fits.open(args.fits, memmap=True) as hdul:
        data = hdul[0].data  # frequency, stokes, y = m, x = l
        header = hdul[0].header.copy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(header)
        wcs.wcs.set()
        n_freq, n_pol, n_m, n_l = data.shape

        check("axes l m frequency polarization time, as the FITS cube's",
              [axis["name"] for axis in image["axes"]] == ["l", "m", "frequency", "polarization", "time"]
              and [axis["length"] for axis in image["axes"]] == [n_l, n_m, n_freq, n_pol, 1])
        if failures:
            return finish()

        # The layout as carta-zarr reads it, in the image's axis order whatever the store's.
        l, m, f = (int(part) for part in args.chunks.split(","))
        chunk = {"l": l, "m": m, "frequency": f, "polarization": 1, "time": 1}
        shard = dict(chunk)
        if args.shards:
            sl, sm, sf = (int(part) for part in args.shards.split(","))
            shard.update(l=sl, m=sm, frequency=sf)
        check(f"chunk {args.chunks}" + (f", shard {args.shards}" if args.shards else ", unsharded"),
              all(axes[name]["chunk"] == chunk[name] and axes[name]["shard"] == shard[name] for name in chunk)
              and image["sharded"] == bool(args.shards))
        check("float32", described["stored_type"] == "float32")
        check("brightness unit", described["unit"] == header.get("BUNIT"))

        mask = described["pixel_mask_id"]
        if args.flag == "yes":
            check("a pixel mask carta-zarr applies", image["has_pixel_mask"] and bool(mask))
        else:
            check("no pixel mask", not image["has_pixel_mask"])

        # The direction coordinate carta-zarr builds, against the FITS file's own.
        direction = described.get("direction") or {}
        pv = {(i, m): value for i, m, value in wcs.wcs.get_pv()}
        check("direction coordinate",
              direction.get("projection") == wcs.wcs.ctype[0][-3:] and wcs.wcs.ctype[1][-3:] == wcs.wcs.ctype[0][-3:]
              and direction.get("reference_frame", "").upper() == wcs.wcs.radesys.upper()
              and (close(direction.get("equinox"), wcs.wcs.equinox)
                   or (direction.get("equinox") is None and math.isnan(wcs.wcs.equinox)))
              and close(direction.get("reference_pixel"), wcs.wcs.crpix[:2])
              and close(direction.get("reference_value"), wcs.wcs.crval[:2])
              and close(direction.get("increment"), wcs.wcs.cdelt[:2])
              and close(direction.get("transformation_matrix"), wcs.wcs.get_pc()[:2, :2])
              and close(direction.get("projection_parameters") or [0.0, 0.0], [pv.get((2, 1), 0.0), pv.get((2, 2), 0.0)])
              and close(direction.get("native_pole_direction"), [wcs.wcs.lonpole, wcs.wcs.latpole]))

        # Every channel's frequency, and what the frequencies mean.
        spectral = described.get("spectral") or {}
        world = wcs.pixel_to_world_values(np.zeros(n_freq), np.zeros(n_freq), np.zeros(n_freq), np.arange(n_freq))[3]
        check("frequency of every channel", close(spectral.get("channel_frequencies"), world))
        check("frequency unit, frame and rest frequency",
              spectral.get("unit") == header.get("CUNIT4") and spectral.get("system", "").upper() == str(header.get("SPECSYS", "")).upper()
              and "RESTFRQ" in header and close(spectral.get("rest_frequency"), header["RESTFRQ"], 1e-12))
        codes = wcs.wcs.crval[2] + wcs.wcs.cdelt[2] * (np.arange(n_pol) + 1 - wcs.wcs.crpix[2])
        check("polarization", described.get("polarization") == [STOKES.get(int(round(code)), "?") for code in codes])

        # The time axis and the observation date, which carta-zarr takes from different places.
        observed = Time(header["DATE-OBS"], scale=str(header.get("TIMESYS", "UTC")).lower()).utc.mjd
        temporal = described.get("temporal") or {}
        check("time", temporal.get("format", "").upper() == "MJD" and temporal.get("scale", "").upper() == "UTC"
              and temporal.get("unit") == "d" and close(temporal.get("values"), [observed], 1e-12))
        observation = described.get("observation") or {}
        check("observation date", observation.get("timesys", "").upper() == "UTC"
              and close(observation.get("mjd_obs"), observed, 1e-12))

        # The restoring beam a Jy/beam image is calibrated by, on every plane; none where FITS has none.
        beams = described.get("beams", [])
        if "BMAJ" in header:
            expected = [math.radians(header["BMAJ"]), math.radians(header["BMIN"]), math.radians(header.get("BPA", 0.0))]
            check("restoring beam of every plane",
                  len(beams) == n_freq * n_pol and all(beam["unit"] == "rad" for beam in beams)
                  and close([[beam["major"], beam["minor"], beam["position_angle"]] for beam in beams], [expected] * len(beams)))
        else:
            check("no restoring beam, as the FITS file has none", not beams)

        # The pixels, and the mask carta-zarr applies, as stored. The mask is read beside the pixels, so
        # its layout is part of what a layout name promises.
        sky = zarr.open_array(f"{args.zarr}/SKY", mode="r")
        flag = zarr.open_array(f"{args.zarr}/{mask}", mode="r") if image["has_pixel_mask"] and mask else None
        if flag is not None:
            check(f"{mask} in the same chunk" + (" and shard" if args.shards else ""),
                  flag.chunks == sky.chunks and flag.shards == sky.shards)
        channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))
        first_pixels = first_flag = None
        for start in range(0, n_freq, channels):
            stop = min(n_freq, start + channels)
            expected_pixels = np.ascontiguousarray(np.asarray(data[start:stop]).transpose(0, 1, 3, 2), dtype=np.float32)
            stored = np.asarray(sky[0, start:stop])  # frequency, polarization, l, m
            if first_pixels is None:
                where = first_difference(expected_pixels.view(np.uint32), stored.view(np.uint32))
                first_pixels = None if where is None else (start + where[0], *where[1:])
            if flag is not None and first_flag is None:
                where = first_difference(np.asarray(flag[0, start:stop], dtype=bool), np.isnan(expected_pixels))
                first_flag = None if where is None else (start + where[0], *where[1:])
        check("every pixel" + (f" (first differs at channel, stokes, x, y {first_pixels})" if first_pixels else ""),
              first_pixels is None)
        if flag is not None:
            check("mask where NaN" + (f" (first differs at {first_flag})" if first_flag else ""), first_flag is None)

    return finish()


if __name__ == "__main__":
    sys.exit(main())
