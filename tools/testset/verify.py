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
channel's frequency (in Hz, whatever unit each is written in), the spectral frame and rest
frequency, the Stokes parameters, the time and the observation date match; and so does the restoring beam of every plane. Then every pixel, as
carta-zarr decodes it (`carta-zarr-bench pixels`), is compared bit for bit, a block of channels at a
time; with the pixel mask applied it must be NaN exactly where the FITS cube is, and the flag itself,
as stored, true exactly there and in the same chunks as the pixels.

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
from astropy import units
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS

STOKES = {1: "I", 2: "Q", 3: "U", 4: "V"}
# The FITS axis types of a direction in each frame carta-zarr reports.
CELESTIAL = {"ICRS": ("RA", "DEC"), "FK5": ("RA", "DEC"), "FK4": ("RA", "DEC"), "GALACTIC": ("GLON", "GLAT")}


def first_difference(a: np.ndarray, b: np.ndarray) -> tuple[int, ...] | None:
    """Where `a` and `b` first differ, without listing everywhere they do."""
    differ = a != b
    if not differ.any():
        return None
    return tuple(int(i) for i in np.unravel_index(int(np.argmax(differ)), differ.shape))


def dimension_names(array: zarr.Array) -> list[str]:
    """The array's dimension names where carta-zarr finds them: the field, or failing that the attribute
    some XRADIO writers keep them in."""
    return list(array.metadata.dimension_names or array.attrs["dimension_names"])


def layout(array: zarr.Array) -> dict[str, list[int]]:
    """Each axis's extent in the outer chunk and in every chunk nested inside it, by name: what decides
    which pixels a chunk holds. A transpose ahead of a sharding codec reorders the shape that codec's
    chunk_shape is given in, so zarr-python's .chunks and .shards can read alike for two different
    footprints; one inside it reorders bytes within a chunk and changes no footprint."""
    names = dimension_names(array)
    metadata = array.metadata.to_dict()
    extents = {name: [extent] for name, extent in zip(names, metadata["chunk_grid"]["configuration"]["chunk_shape"])}

    def descend(chain, order: list[int]) -> None:
        for codec in chain:
            name, configuration = codec.get("name"), codec.get("configuration") or {}
            if name == "transpose":
                order = [order[axis] for axis in configuration["order"]]
            elif name == "sharding_indexed":
                for axis, extent in zip(order, configuration["chunk_shape"]):
                    extents[names[axis]].append(extent)
                descend(configuration.get("codecs", []), order)

    descend(metadata.get("codecs", []), list(range(len(names))))
    return extents


def logical(array: zarr.Array, start: int, stop: int) -> np.ndarray:
    """Channels start to stop of a sky-plane array, as frequency, polarization, l, m whatever order its
    dimensions are stored in -- which carta-zarr reads by name, so a check that indexed by position would
    compare a square cube with l and m swapped against the wrong pixels and pass it."""
    names = dimension_names(array)
    selection = tuple(0 if name == "time" else slice(start, stop) if name == "frequency" else slice(None)
                      for name in names)
    kept = [name for name in names if name != "time"]
    return np.asarray(array[selection]).transpose([kept.index(name) for name in ("frequency", "polarization", "l", "m")])


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
        frame = str(direction.get("reference_frame", "")).upper()
        pv = {(i, m): value for i, m, value in wcs.wcs.get_pv()}
        # The linear part of the FITS WCS whole, increment times matrix, so that CDELT and PC, or CD,
        # alike describe it; carta-zarr's direction is the sky axes' alone, so a FITS file whose sky
        # axes mix with the others -- PC1_4, say -- means something it cannot, and is refused.
        linear = np.diag(wcs.wcs.get_cdelt()) @ wcs.wcs.get_pc()
        # Nor can its Stokes or frequency axis vary along any other: Stokes I in every channel is not
        # Stokes I in the first and V by the thirtieth.
        check("each axis independent of the others but l and m of each other",
              not linear[:2, 2:].any() and not linear[2:, :2].any() and linear[2, 3] == 0 and linear[3, 2] == 0)
        check("direction coordinate",
              direction.get("projection") == wcs.wcs.ctype[0][-3:] and wcs.wcs.ctype[1][-3:] == wcs.wcs.ctype[0][-3:]
              and (frame == "GALACTIC" or frame == wcs.wcs.radesys.upper())
              and [ctype[:4].rstrip("-") for ctype in list(wcs.wcs.ctype)[:2]] == list(CELESTIAL.get(frame, ("?", "?")))
              # Galactic coordinates have no equinox, whatever the dataset carries beside them.
              and (frame == "GALACTIC" or close(direction.get("equinox"), wcs.wcs.equinox)
                   or (direction.get("equinox") is None and math.isnan(wcs.wcs.equinox)))
              and close(direction.get("reference_pixel"), wcs.wcs.crpix[:2])
              and close(direction.get("reference_value"), wcs.wcs.crval[:2])
              and close(np.diag(direction.get("increment") or [np.nan] * 2) @ np.asarray(direction.get("transformation_matrix") or np.nan, dtype=float),
                        linear[:2, :2])
              and close(direction.get("projection_parameters") or [0.0, 0.0], [pv.get((2, 1), 0.0), pv.get((2, 2), 0.0)])
              and close(direction.get("native_pole_direction"), [wcs.wcs.lonpole, wcs.wcs.latpole]))

        # Every channel's frequency, and what the frequencies mean.
        spectral = described.get("spectral") or {}
        world = wcs.pixel_to_world_values(np.zeros(n_freq), np.zeros(n_freq), np.zeros(n_freq), np.arange(n_freq))[3]
        # astropy gives the FITS frequencies in Hz whatever CUNIT4 says, and RESTFRQ is in Hz always;
        # carta-zarr gives both in the dataset's unit.
        try:
            hertz = units.Unit(spectral.get("unit", "")).to(units.Hz)
        except (TypeError, ValueError):
            hertz = math.nan
        check("a FITS frequency axis, in Hz as astropy gives it",
              str(header.get("CTYPE4", "")).startswith("FREQ") and str(wcs.world_axis_units[3]) == "Hz")
        check("frequency of every channel", close(np.asarray(spectral.get("channel_frequencies"), dtype=float) * hertz, world))
        # The linear axis carta-backend builds when carta-zarr gives one, rather than the table. A FITS
        # frequency axis is linear, so carta-zarr must give one -- unless there is one channel, which has
        # no spacing to fit and is read from the table.
        try:
            fitted = spectral["reference_value"] + (np.arange(n_freq) + 1 - spectral["reference_pixel"]) * spectral["increment"]
        except (KeyError, TypeError):
            fitted = None
        if fitted is not None or n_freq > 1:
            check("frequency axis as carta-backend builds it", fitted is not None and close(fitted * hertz, world))
        check("frequency frame and rest frequency",
              math.isfinite(hertz) and spectral.get("system", "").upper() == str(header.get("SPECSYS", "")).upper()
              and "RESTFRQ" in header and spectral.get("rest_frequency") is not None
              and close(spectral["rest_frequency"] * hertz, header["RESTFRQ"], 1e-12))
        codes = wcs.wcs.crval[2] + linear[2, 2] * (np.arange(n_pol) + 1 - wcs.wcs.crpix[2])
        check("polarization", header.get("CTYPE3") == "STOKES" and np.allclose(codes, np.round(codes), rtol=0, atol=1e-9)
              and described.get("polarization") == [STOKES.get(int(round(code)), "?") for code in codes])

        # The time axis and the observation date, which carta-zarr takes from different places.
        observed = Time(header["DATE-OBS"], scale=str(header.get("TIMESYS", "UTC")).lower()).utc.mjd
        temporal = described.get("temporal") or {}
        check("time", temporal.get("format", "").upper() == "MJD" and temporal.get("scale", "").upper() == "UTC"
              and temporal.get("unit") == "d" and close(temporal.get("values"), [observed], 1e-12))
        observation = described.get("observation") or {}
        check("observation date", observation.get("timesys", "").upper() == "UTC"
              and close(observation.get("mjd_obs"), observed, 1e-12))
        # The observatory a spectral frame is placed by: OBSGEO-X, Y and Z, and its name.
        if "OBSGEO-X" in header:
            check("observatory", observation.get("telescope_name") == header.get("TELESCOP")
                  and close(observation.get("observatory_position"), [header["OBSGEO-X"], header["OBSGEO-Y"], header["OBSGEO-Z"]]))
        else:
            check("no observatory, as the FITS file has none", observation.get("observatory_position") is None)

        # The restoring beam a Jy/beam image is calibrated by, on every plane; none where FITS has none.
        beams = described.get("beams", [])
        if "BMAJ" in header:
            expected = [math.radians(header["BMAJ"]), math.radians(header["BMIN"]), math.radians(header.get("BPA", 0.0))]
            planes = sorted((beam["time"], beam["channel"], beam["polarization"]) for beam in beams)
            check("restoring beam of every plane",
                  planes == [(0, c, p) for c in range(n_freq) for p in range(n_pol)]
                  and all(beam["unit"] == "rad" for beam in beams)
                  and close([[beam["major"], beam["minor"], beam["position_angle"]] for beam in beams], [expected] * len(beams)))
        else:
            check("no restoring beam, as the FITS file has none", not beams)

        # The pixels as carta-zarr decodes them, stored and with its pixel mask applied, so a codec it
        # cannot decode, or a mask it applies elsewhere than the NaN, is refused however zarr-python reads
        # them. The flag itself is read as stored too: a flag true exactly where the FITS cube is NaN is
        # what keeps the two formats' statistics alike, and its layout is part of what a layout name
        # promises, since the mask is read beside the pixels.
        def read(start: int, stop: int, masked: bool) -> np.ndarray:
            """Channels start to stop through carta-zarr, as frequency, polarization, l, m."""
            command = [args.bench, "pixels", args.zarr, "--channels", f"{start}:{stop}"] + ([] if masked else ["--unmasked"])
            result = subprocess.run(command, capture_output=True)
            if result.returncode:
                # TensorStore's own source locations and spec follow the reason, at length.
                raise RuntimeError(result.stderr.decode(errors="replace").strip().split(" [source locations=")[0])
            return np.frombuffer(result.stdout, dtype=np.float32).reshape(n_pol, stop - start, n_m, n_l).transpose(1, 0, 3, 2)

        sky = zarr.open_array(f"{args.zarr}/SKY", mode="r")
        flag = zarr.open_array(f"{args.zarr}/{mask}", mode="r") if image["has_pixel_mask"] and mask else None
        if flag is not None:
            check(f"{mask} in the same chunk" + (" and shard" if args.shards else ""), layout(flag) == layout(sky))
        channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))
        first_pixels = first_masked = first_flag = unreadable = None
        for start in range(0, n_freq, channels):
            stop = min(n_freq, start + channels)
            expected_pixels = np.ascontiguousarray(np.asarray(data[start:stop]).transpose(0, 1, 3, 2), dtype=np.float32)
            try:
                stored = np.ascontiguousarray(read(start, stop, masked=False))
                applied = read(start, stop, masked=True) if flag is not None else None
            except RuntimeError as error:
                unreadable = f"channels {start}:{stop}: {error}"
                break
            if first_pixels is None:
                where = first_difference(expected_pixels.view(np.uint32), stored.view(np.uint32))
                first_pixels = None if where is None else (start + where[0], *where[1:])
            if flag is not None and first_masked is None:
                # NaN where the FITS cube is NaN, and the stored pixel, bit for bit, everywhere else.
                nan = np.isnan(applied)
                wrong = (nan != np.isnan(expected_pixels)) | (~nan & (applied.view(np.uint32) != stored.view(np.uint32)))
                where = first_difference(wrong, np.zeros_like(wrong))
                first_masked = None if where is None else (start + where[0], *where[1:])
            if flag is not None and first_flag is None:
                where = first_difference(logical(flag, start, stop).astype(bool), np.isnan(expected_pixels))
                first_flag = None if where is None else (start + where[0], *where[1:])
        check("carta-zarr reads every pixel" + (f" ({unreadable})" if unreadable else ""), unreadable is None)
        if unreadable is None:
            check("every pixel" + (f" (first differs at channel, stokes, x, y {first_pixels})" if first_pixels else ""),
                  first_pixels is None)
            if flag is not None:
                check("mask where NaN" + (f" (first differs at {first_flag})" if first_flag else ""), first_flag is None)
                check("pixels masked where NaN" + (f" (first differs at {first_masked})" if first_masked else ""),
                      first_masked is None)

    return finish()


if __name__ == "__main__":
    sys.exit(main())
