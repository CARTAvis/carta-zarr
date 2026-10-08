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

A comparison between FITS and Zarr is only a comparison of the formats if both hold the same
pixels. This compares every pixel, bit for bit, a block of channels at a time; the frequency of the
first and last channel; and the sky position of the corners and the centre. --flag says whether the
Zarr must carry a flag: one that must is refused without one carta-zarr would apply -- typed "flag",
boolean, shaped as the image, and declared for it or the only such variable -- and its flag must be
true exactly where the FITS cube is NaN; one that must not is refused with one.

    verify.py CUBE.fits CUBE.zarr --flag {yes,no} [--block-mib 1024]

Prints what it compared, and exits nonzero if anything differs.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import zarr
from astropy.io import fits
from astropy.wcs import WCS


def usable_flag(root: Path, image: str) -> tuple[str | None, str]:
    """The flag carta-zarr would apply to `image`, by its rules (src/schema/xradio/flag.h), and why
    there is none when there is none."""
    sky = json.loads((root / image / "zarr.json").read_text())

    def usable(name: str) -> bool:
        try:
            metadata = json.loads((root / name / "zarr.json").read_text())
        except OSError:
            return False
        attributes = metadata.get("attributes", {})
        boolean = metadata.get("data_type") == "bool" or (
            metadata.get("data_type") == "int8" and attributes.get("dtype") == "bool")
        return (attributes.get("type") == "flag" and boolean and metadata.get("node_type") == "array"
                and metadata.get("dimension_names") == sky.get("dimension_names") and metadata.get("shape") == sky.get("shape"))

    groups = json.loads((root / "zarr.json").read_text()).get("attributes", {}).get("data_groups") or {}
    declared = sky.get("attributes", {}).get("flag") or next(
        (group.get("flag") for group in groups.values() if isinstance(group, dict) and group.get("sky") == image), None)
    if declared:
        return (declared, "") if usable(declared) else (None, f"{declared} is declared but is not a usable flag")
    candidates = [entry.name for entry in root.iterdir() if entry.is_dir() and entry.name != image and usable(entry.name)]
    if len(candidates) == 1:
        return candidates[0], ""
    return None, "no flag is declared, and no single variable is typed flag and shaped as the image"


def first_difference(a: np.ndarray, b: np.ndarray) -> tuple[int, ...] | None:
    """Where `a` and `b` first differ, without listing everywhere they do."""
    differ = a != b
    if not differ.any():
        return None
    return tuple(int(i) for i in np.unravel_index(int(np.argmax(differ)), differ.shape))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fits")
    parser.add_argument("zarr")
    parser.add_argument("--flag", choices=("yes", "no"), required=True, help="whether the Zarr must carry a flag")
    parser.add_argument("--block-mib", type=int, default=1024, help="how much of each to compare at a time")
    args = parser.parse_args()

    sky = zarr.open_array(f"{args.zarr}/SKY", mode="r")
    n_time, n_freq, n_pol, n_l, n_m = sky.shape
    failures = []

    def check(what: str, ok: bool) -> None:
        print(f"{what}: {'same' if ok else 'DIFFERENT'}", flush=True)
        if not ok:
            failures.append(what)

    flag_name, why = usable_flag(Path(args.zarr), "SKY")
    flag = zarr.open_array(f"{args.zarr}/{flag_name}", mode="r") if flag_name else None
    if args.flag == "yes":
        check(f"a flag carta-zarr applies{'' if flag_name else f' ({why})'}", flag is not None)
    else:
        check("no flag", flag is None)

    with fits.open(args.fits, memmap=True) as hdul:
        data = hdul[0].data  # frequency, stokes, y = m, x = l
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(hdul[0].header)
        check("shape", n_time == 1 and data.shape == (n_freq, n_pol, n_m, n_l))
        if failures:
            print(f"MISMATCH: {', '.join(failures)}")
            return 1
        channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))
        first_pixels = first_flag = None
        for start in range(0, n_freq, channels):
            stop = min(n_freq, start + channels)
            expected = np.ascontiguousarray(np.asarray(data[start:stop]).transpose(0, 1, 3, 2), dtype=np.float32)
            stored = np.asarray(sky[0, start:stop])  # frequency, polarization, l, m
            if first_pixels is None:
                where = first_difference(expected.view(np.uint32), stored.view(np.uint32))
                first_pixels = None if where is None else (start + where[0], *where[1:])
            if flag is not None and first_flag is None:
                where = first_difference(np.asarray(flag[0, start:stop], dtype=bool), np.isnan(expected))
                first_flag = None if where is None else (start + where[0], *where[1:])
        check("every pixel" + (f" (first differs at channel, stokes, x, y {first_pixels})" if first_pixels else ""),
              first_pixels is None)
        if flag is not None:
            check("flag where NaN" + (f" (first differs at {first_flag})" if first_flag else ""), first_flag is None)

    frequency = zarr.open_array(f"{args.zarr}/frequency", mode="r")[...]
    for f in sorted({0, n_freq - 1}):
        world = wcs.pixel_to_world_values(0, 0, 0, f)[3]
        check(f"frequency of channel {f}", abs(world - frequency[f]) <= 1e-9 * abs(frequency[f]))
    # xradio's converter writes l and m but not the sky position of every pixel, so the position is
    # the SIN projection of l and m from the reference direction, inverted as XRADIO's would be.
    root = json.loads(Path(args.zarr, "zarr.json").read_text())["attributes"]
    ra0, dec0 = root["coordinate_system_info"]["reference_direction"]["data"]
    l_axis = zarr.open_array(f"{args.zarr}/l", mode="r")[...]
    m_axis = zarr.open_array(f"{args.zarr}/m", mode="r")[...]
    for x, y in [(0, 0), (n_l - 1, n_m - 1), (n_l // 2, n_m // 2)]:
        l, m = float(l_axis[x]), float(m_axis[y])
        n = math.sqrt(max(0.0, 1.0 - l * l - m * m))
        dec = math.asin(m * math.cos(dec0) + n * math.sin(dec0))
        ra = ra0 + math.atan2(l, n * math.cos(dec0) - m * math.sin(dec0))
        world_ra, world_dec = wcs.celestial.pixel_to_world_values(x, y)
        close = abs((world_ra - math.degrees(ra) + 180) % 360 - 180) < 1e-8 and abs(world_dec - math.degrees(dec)) < 1e-8
        check(f"sky position of ({x}, {y})", bool(close))

    print("OK" if not failures else f"MISMATCH: {', '.join(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
