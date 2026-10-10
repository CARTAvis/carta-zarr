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
pixels. This checks that the image is float32 over time, frequency, polarization, l and m, in the
chunks (and shards) asked for; compares every pixel, bit for bit, a block of channels at a time; the
frequency of every channel; and the sky position of the corners and the centre, from an unrotated
SIN projection, which is checked to be what the Zarr declares. --flag says whether the Zarr must carry
a flag: one that must is refused without one carta-zarr would apply -- typed "flag", boolean, shaped
as the image, and declared for it or the only such variable not declared for another -- and its flag
must be true exactly where the FITS cube is NaN; one that must not is refused with one, and with a
declaration carta-zarr would refuse to open the image over.

    verify.py CUBE.fits CUBE.zarr --flag {yes,no} --chunks L,M,F [--shards L,M,F] [--block-mib 1024]

Prints what it compared, and exits nonzero if anything differs.
"""

from __future__ import annotations

import argparse
import json
import math
import posixpath
import sys
import warnings
from pathlib import Path

import numpy as np
import zarr
from astropy.io import fits
from astropy.wcs import WCS


AXES = ["time", "frequency", "polarization", "l", "m"]


def canonical(name: str) -> str:
    """A node a document names, as the store files it, the way carta-zarr compares them: "./FLAG" is
    FLAG. A name that leaves the store is kept as written."""
    normal = posixpath.normpath(name.lstrip("/"))
    return name if normal.startswith("..") else normal


def usable_flag(root: Path, image: str) -> tuple[str | None, str, bool]:
    """The flag carta-zarr would apply to `image`, by its rules (src/schema/xradio/flag.h, DeclaredFlag
    and DetermineFlag); why there is none when there is none; and whether that is a declaration
    carta-zarr refuses to open the image over rather than an image it opens unmasked."""
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
    groups = [group for group in groups.values() if isinstance(group, dict) and group.get("flag")]
    # The image's own attribute outranks a group; groups naming different flags for it are refused.
    declared = sky.get("attributes", {}).get("flag")
    if not declared:
        named = {canonical(group["flag"]) for group in groups if canonical(group.get("sky") or "") == image}
        if len(named) > 1:
            return None, f"data groups declare different flags {sorted(named)}", True
        declared = next(iter(named), None)
    if declared:
        declared = canonical(declared)
        if usable(declared):
            return declared, "", False
        return None, f"{declared} is declared but is not a usable flag", True
    # With nothing declared, the one usable variable no group declares for another image.
    others = {canonical(group["flag"]) for group in groups if canonical(group.get("sky") or "") != image}
    candidates = [entry.name for entry in root.iterdir()
                  if entry.is_dir() and entry.name != image and entry.name not in others and usable(entry.name)]
    if len(candidates) == 1:
        return candidates[0], "", False
    return None, "no flag is declared, and no single variable is typed flag and shaped as the image", False


def layout(metadata: dict) -> tuple[list[int], list[int] | None]:
    """An array's inner chunk and, when it is sharded, its shard, as stored."""
    grid = metadata["chunk_grid"]["configuration"]["chunk_shape"]
    for codec in metadata.get("codecs", []):
        if codec.get("name") == "sharding_indexed":
            return codec["configuration"]["chunk_shape"], grid
    return grid, None


def extents(text: str) -> list[int]:
    """L,M,F as the image's own axis order."""
    l, m, f = (int(part) for part in text.split(","))
    return [1, f, 1, l, m]


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
    parser.add_argument("--chunks", required=True, help="the chunk asked for, as l,m,frequency")
    parser.add_argument("--shards", help="the shard asked for, as l,m,frequency; none means unsharded")
    parser.add_argument("--block-mib", type=int, default=1024, help="how much of each to compare at a time")
    args = parser.parse_args()

    sky = zarr.open_array(f"{args.zarr}/SKY", mode="r")
    n_time, n_freq, n_pol, n_l, n_m = sky.shape
    failures = []

    def check(what: str, ok: bool) -> None:
        print(f"{what}: {'same' if ok else 'DIFFERENT'}", flush=True)
        if not ok:
            failures.append(what)

    # The pixels are compared as bits in this axis order, which carta-zarr takes from the names.
    metadata = json.loads(Path(args.zarr, "SKY", "zarr.json").read_text())
    check("float32", metadata.get("data_type") == "float32")
    check(f"axes {' '.join(AXES)}", metadata.get("dimension_names") == AXES)
    # A converter that ignored the layout asked of it would be published under the wrong name.
    inner, shard = layout(metadata)
    check(f"chunk {args.chunks}" + (f", shard {args.shards}" if args.shards else ", unsharded"),
          inner == extents(args.chunks) and shard == (extents(args.shards) if args.shards else None))

    flag_name, why, refused = usable_flag(Path(args.zarr), "SKY")
    flag = zarr.open_array(f"{args.zarr}/{flag_name}", mode="r") if flag_name else None
    if args.flag == "yes":
        check(f"a flag carta-zarr applies{'' if flag_name else f' ({why})'}", flag is not None)
    else:
        check("no flag" + (f" ({why}, so carta-zarr refuses the image)" if refused else ""), flag is None and not refused)

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
    channels = np.arange(n_freq)
    world = wcs.pixel_to_world_values(np.zeros(n_freq), np.zeros(n_freq), np.zeros(n_freq), channels)[3]
    wrong = np.flatnonzero(~(np.isfinite(frequency) & (np.abs(world - frequency) <= 1e-9 * np.abs(frequency))))
    check("frequency of every channel" + (f" (first differs at channel {wrong[0]})" if wrong.size else ""), wrong.size == 0)
    # xradio's converter writes l and m but not the sky position of every pixel, so the position is
    # the SIN projection of l and m from the reference direction, inverted as XRADIO's would be. That
    # holds only for the projection the Zarr declares being an unrotated SIN with no parameters, in the
    # FITS file's frame, so that is checked first.
    root = json.loads(Path(args.zarr, "zarr.json").read_text())["attributes"]
    system = root["coordinate_system_info"]
    frame = system["reference_direction"]["attrs"].get("frame", "")
    wcs.wcs.set()
    pole = [math.degrees(angle) for angle in system.get("native_pole_direction", {}).get("data", [math.nan, math.nan])]
    check("unrotated SIN projection in the FITS frame, about the FITS file's native pole",
          system.get("projection") == "SIN" and system.get("pixel_coordinate_transformation_matrix") == [[1.0, 0.0], [0.0, 1.0]]
          and not any(system.get("projection_parameters", [])) and frame.lower() == wcs.wcs.radesys.lower()
          and all(t.endswith("-SIN") for t in list(wcs.wcs.ctype)[:2])
          and np.allclose(pole, [wcs.wcs.lonpole, wcs.wcs.latpole], rtol=0, atol=1e-9))
    ra0, dec0 = system["reference_direction"]["data"]
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
