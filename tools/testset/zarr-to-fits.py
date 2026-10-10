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
# # As generate.py: never the system interpreter.
# python-preference = "only-managed"
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Write an XRADIO image dataset's sky image as a FITS cube holding the same pixels.

The test set's cubes are synthesized as XRADIO datasets by generate.py, but the Zarr a site reads is
what its FITS-to-Zarr converter made from a FITS file. This writes that FITS file, so that the
test set's Zarr can be made the way a site's is, and so that there is a FITS cube to compare it with.

The axes are RA, Dec, Stokes and frequency (NAXIS1 to 4), the order ASKAPsoft writes, so FITS pixel
(x, y, stokes, channel) is the dataset's (l, m, polarization, frequency). The header carries what
xradio's FITS reader requires of one: CUNIT on every axis, LONPOLE and LATPOLE, PV2_1 and PV2_2, and
DATE-OBS. Channels are written a block at a time, so a cube larger than memory can be written.

    zarr-to-fits.py DATASET OUTPUT.fits [--block-mib 1024]

The file appears at OUTPUT.fits only once it is complete.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import zarr
from astropy.io import fits
from astropy.time import Time

STOKES_CODE = {"I": 1, "Q": 2, "U": 3, "V": 4}
SKY_AXES = ["time", "frequency", "polarization", "l", "m"]


def metadata(path: Path) -> dict[str, Any]:
    return json.loads((path / "zarr.json").read_text())


def nodes_under(root: Path) -> list[tuple[str, dict[str, Any]]]:
    """Every node below the root and its zarr.json, following directory links as a store's reader
    does, each directory once however many ways it is reached."""
    found, seen = [], set()
    for directory, subdirectories, files in os.walk(root, followlinks=True):
        real = os.path.realpath(directory)
        if real in seen:
            subdirectories.clear()
            continue
        seen.add(real)
        if "zarr.json" in files and Path(directory) != root:
            found.append((str(Path(directory).relative_to(root)), json.loads((Path(directory) / "zarr.json").read_text())))
    return found


def header_of(root: Path, image: str, sky: Any) -> fits.Header:
    attributes = metadata(root)["attributes"]
    _, n_freq, n_pol, n_l, n_m = sky.shape
    l = zarr.open_array(str(root / "l"), mode="r")[...]
    m = zarr.open_array(str(root / "m"), mode="r")[...]
    freq = zarr.open_array(str(root / "frequency"), mode="r")[...]
    pol = zarr.open_array(str(root / "polarization"), mode="r")[...]
    system = attributes["coordinate_system_info"]
    ra0, dec0 = system["reference_direction"]["data"]
    projection = system.get("projection", "SIN")
    sky_attributes = sky.attrs.asdict()

    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -32
    header["NAXIS"] = 4
    header["NAXIS1"] = n_l
    header["NAXIS2"] = n_m
    header["NAXIS3"] = n_pol
    header["NAXIS4"] = n_freq
    header["EXTEND"] = True
    header["BUNIT"] = sky_attributes.get("units", "Jy/beam")
    header["OBJECT"] = sky_attributes.get("object_name", "")
    # l is x times the increment from the reference direction, so the reference pixel is where l is 0.
    for axis, values, ctype, crval in ((1, l, "RA---", ra0), (2, m, "DEC--", dec0)):
        increment = float(values[1] - values[0])
        header[f"CTYPE{axis}"] = ctype + projection
        header[f"CRVAL{axis}"] = math.degrees(crval)
        header[f"CDELT{axis}"] = math.degrees(increment)
        header[f"CRPIX{axis}"] = 1.0 - float(values[0]) / increment
        header[f"CUNIT{axis}"] = "deg"
    header["CTYPE3"] = "STOKES"
    header["CRVAL3"] = float(STOKES_CODE[str(pol[0])])
    header["CDELT3"] = 1.0
    header["CRPIX3"] = 1.0
    header["CUNIT3"] = ""
    header["CTYPE4"] = "FREQ"
    header["CRVAL4"] = float(freq[0])
    header["CDELT4"] = float(freq[1] - freq[0]) if n_freq > 1 else 1.0
    header["CRPIX4"] = 1.0
    header["CUNIT4"] = "Hz"
    header["RADESYS"] = "FK5"
    header["EQUINOX"] = 2000.0
    header["LONPOLE"] = 180.0
    header["LATPOLE"] = math.degrees(dec0)
    header["PV2_1"] = 0.0
    header["PV2_2"] = 0.0
    header["SPECSYS"] = "LSRK"
    rest = metadata(root / "frequency").get("attributes", {}).get("rest_frequency", {}).get("data")
    if rest is not None:
        header["RESTFRQ"] = float(rest)
    # ASKAP, as the frequency axis is: xradio wants an observatory and a date to place it. The date is
    # the dataset's own, so that the FITS file and the dataset agree on it as on everything else.
    header["TELESCOP"] = "ASKAP"
    header["DATE-OBS"] = Time(float(zarr.open_array(str(root / "time"), mode="r")[0]), format="mjd", scale="utc").isot
    header["TIMESYS"] = "UTC"
    header["VELREF"] = 257
    header["OBSGEO-X"] = -2.558266717765e06
    header["OBSGEO-Y"] = 5.095672176508e06
    header["OBSGEO-Z"] = -2.849020838078e06
    header["BTYPE"] = "Intensity"
    beam_name = attributes.get("data_groups", {}).get("base", {}).get("beam_fit_params_sky")
    if beam_name and (root / beam_name).exists():
        beam = zarr.open_array(str(root / beam_name), mode="r")[0, 0, 0, :]
        header["BMAJ"] = math.degrees(float(beam[0]))
        header["BMIN"] = math.degrees(float(beam[1]))
        header["BPA"] = math.degrees(float(beam[2]))
    header["ORIGIN"] = "carta-zarr tools/testset/zarr-to-fits.py"
    header["HISTORY"] = f"pixels of {root.name}/{image}"
    return header


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset", help="an XRADIO image dataset")
    parser.add_argument("output", help="the FITS file to write")
    parser.add_argument("--block-mib", type=int, default=1024, help="how much to read and write at a time")
    args = parser.parse_args()

    root = Path(args.dataset)
    image = metadata(root)["attributes"]["data_groups"]["base"]["sky"]
    sky = zarr.open_array(str(root / image), mode="r")
    if list(sky.metadata.dimension_names) != SKY_AXES:
        raise SystemExit(f"{image} is stored as {sky.metadata.dimension_names}; expected {SKY_AXES}")
    if sky.shape[0] != 1:
        raise SystemExit(f"{image} has {sky.shape[0]} times; a FITS cube holds one")
    # A FITS cube marks a flagged pixel by NaN alone, and a flagged dataset may hold finite values under
    # its flag (generate.py --flag does), which would be written as valid pixels. The test set's cubes
    # carry their NaN in the pixels and get a flag only from the converter, so a flag here is refused.
    flags = [node for node, document in nodes_under(root) if document.get("attributes", {}).get("type") == "flag"]
    groups = metadata(root)["attributes"].get("data_groups", {}).values()
    if sky.attrs.get("flag") or any(isinstance(group, dict) and group.get("flag") for group in groups) or flags:
        raise SystemExit(f"{args.dataset} carries a flag ({', '.join(flags) or 'declared'}); FITS would show what it "
                         "hides, so write the cube without one")
    _, n_freq, n_pol, n_l, n_m = sky.shape
    channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))

    out = Path(args.output)
    partial = out.with_name(out.name + ".partial")
    started = last_report = time.monotonic()
    with open(partial, "wb") as f:
        f.write(header_of(root, image, sky).tostring().encode("ascii"))
        written = 0
        for start in range(0, n_freq, channels):
            stop = min(n_freq, start + channels)
            block = np.asarray(sky[0, start:stop])  # frequency, polarization, l, m
            np.ascontiguousarray(block.transpose(0, 1, 3, 2), dtype=">f4").tofile(f)
            written += block.size * 4
            now = time.monotonic()
            if now - last_report >= 5 or stop == n_freq:
                rate = written / (1 << 20) / max(now - started, 1e-9)
                print(f"  {out.name}: {stop}/{n_freq} channels, {written / (1 << 30):.2f} GiB, {rate:.0f} MiB/s",
                      file=sys.stderr, flush=True)
                last_report = now
        f.write(b"\0" * (-written % 2880))
    partial.rename(out)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
