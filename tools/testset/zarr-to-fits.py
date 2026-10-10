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
DATE-OBS. The coordinates are written as carta-zarr reads them -- `carta-zarr-bench probe --describe`,
at --bench or $CARTA_ZARR_BENCH -- and a dataset whose coordinates a FITS header cannot hold is
refused. Channels are written a block at a time, so a cube larger than memory can be written.

    zarr-to-fits.py DATASET OUTPUT.fits [--block-mib 1024] [--bench PATH]

The file appears at OUTPUT.fits only once it is complete.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import zarr
from astropy import units
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


def describe(bench: str, root: Path) -> dict[str, Any]:
    """The default image as carta-zarr reads it: `carta-zarr-bench probe --describe`."""
    probe = subprocess.run([bench, "probe", str(root), "--describe"], capture_output=True, text=True)
    report = json.loads(probe.stdout) if probe.stdout.strip() else {}
    if not report.get("ok"):
        raise SystemExit(f"carta-zarr does not open {root}: {report.get('error') or probe.stderr.strip()}")
    return report["image"]


def header_of(root: Path, image: str, sky: Any, described: dict[str, Any]) -> fits.Header:
    """The FITS header of what carta-zarr reads the image as, so that the FITS file means what the
    dataset does, whatever its coordinates; what a FITS header cannot say -- frequencies that are not
    evenly spaced, Stokes parameters out of sequence, a beam that differs from plane to plane -- is
    refused rather than written as something else."""
    meaning = described["description"]
    lengths = {axis["name"]: axis["length"] for axis in described["axes"]}
    n_l, n_m, n_freq, n_pol = (lengths[name] for name in ("l", "m", "frequency", "polarization"))
    direction, spectral = meaning.get("direction"), meaning.get("spectral")
    if not direction or not spectral:
        raise SystemExit(f"{image} has no direction or no spectral coordinate to write")
    if meaning["stored_type"] != "float32":
        raise SystemExit(f"{image} is {meaning['stored_type']}; the FITS cube is float32")

    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -32
    header["NAXIS"] = 4
    header["NAXIS1"] = n_l
    header["NAXIS2"] = n_m
    header["NAXIS3"] = n_pol
    header["NAXIS4"] = n_freq
    header["EXTEND"] = True
    header["BUNIT"] = meaning["unit"]
    header["OBJECT"] = sky.attrs.get("object_name", "")
    projection = direction["projection"]
    for axis, ctype in ((1, "RA---"), (2, "DEC--")):
        header[f"CTYPE{axis}"] = ctype + projection
        header[f"CRVAL{axis}"] = direction["reference_value"][axis - 1]
        header[f"CDELT{axis}"] = direction["increment"][axis - 1]
        header[f"CRPIX{axis}"] = direction["reference_pixel"][axis - 1]
        header[f"CUNIT{axis}"] = "deg"
    for i in (1, 2):
        for j in (1, 2):
            header[f"PC{i}_{j}"] = direction["transformation_matrix"][i - 1][j - 1]

    labels = meaning.get("polarization") or []
    codes = [STOKES_CODE.get(label) for label in labels]
    if not codes or None in codes or codes != list(range(codes[0], codes[0] + len(codes))):
        raise SystemExit(f"{image}'s Stokes parameters {labels} are not a sequence a FITS axis can hold")
    header["CTYPE3"] = "STOKES"
    header["CRVAL3"] = float(codes[0])
    header["CDELT3"] = 1.0
    header["CRPIX3"] = 1.0
    header["CUNIT3"] = ""

    frequencies = np.asarray(spectral["channel_frequencies"], dtype=float)
    # The increment over the whole axis: one taken from the first two channels is off by their rounding,
    # which 30,000 channels multiply past any tolerance.
    step = float(frequencies[-1] - frequencies[0]) / (n_freq - 1) if n_freq > 1 else 1.0
    if not np.allclose(frequencies, frequencies[0] + step * np.arange(n_freq), rtol=1e-12, atol=0):
        raise SystemExit(f"{image}'s channels are not evenly spaced in frequency, which a FITS axis cannot hold")
    header["CTYPE4"] = "FREQ"
    header["CRVAL4"] = float(frequencies[0])
    header["CDELT4"] = step
    header["CRPIX4"] = 1.0
    header["CUNIT4"] = spectral["unit"]

    header["RADESYS"] = direction["reference_frame"]
    if direction.get("equinox") is not None:
        header["EQUINOX"] = direction["equinox"]
    header["LONPOLE"], header["LATPOLE"] = direction["native_pole_direction"]
    header["PV2_1"], header["PV2_2"] = direction.get("projection_parameters") or [0.0, 0.0]
    header["SPECSYS"] = spectral["system"]
    if spectral.get("rest_frequency") is not None:
        # In the spectral axis's unit as carta-zarr gives it, and in Hz as FITS takes it.
        header["RESTFRQ"] = spectral["rest_frequency"] * units.Unit(spectral["unit"]).to(units.Hz)
    # ASKAP, as the frequency axis is: xradio wants an observatory and a date to place it. The date is
    # the dataset's own, as carta-zarr reads it, and the time axis must agree with it, as it does in
    # a FITS file, which has only the one.
    observation, temporal = meaning.get("observation") or {}, meaning.get("temporal") or {}
    observed = observation.get("mjd_obs")
    scales = {str(observation.get("timesys", "")).upper(), str(temporal.get("scale", "")).upper()}
    if scales != {"UTC"}:
        raise SystemExit(f"{image}'s observation is in {' and '.join(sorted(scales))}; this writes UTC, and a date "
                         "in another scale would be written as the wrong instant")
    if observed is None or temporal.get("values") != [observed]:
        raise SystemExit(f"{image}'s observation date {observed} is not the one time on its time axis "
                         f"{temporal.get('values')}")
    header["TELESCOP"] = "ASKAP"
    header["DATE-OBS"] = Time(observed, format="mjd", scale="utc").isot
    header["TIMESYS"] = "UTC"
    if (velref := {"LSRK": 257, "BARY": 258, "TOPO": 259}.get(spectral["system"].upper())) is not None:
        header["VELREF"] = velref
    header["OBSGEO-X"] = -2.558266717765e06
    header["OBSGEO-Y"] = 5.095672176508e06
    header["OBSGEO-Z"] = -2.849020838078e06
    header["BTYPE"] = "Intensity"
    planes = sorted((beam["time"], beam["channel"], beam["polarization"]) for beam in meaning.get("beams", []))
    if planes and planes != [(0, c, p) for c in range(n_freq) for p in range(n_pol)]:
        raise SystemExit(f"{image}'s restoring beams are not one on each of its planes; a FITS header's beam is "
                         "every plane's")
    beams = {(beam["major"], beam["minor"], beam["position_angle"], beam["unit"]) for beam in meaning.get("beams", [])}
    if len(beams) > 1:
        raise SystemExit(f"{image}'s restoring beam differs from plane to plane; a FITS header holds one")
    if beams:
        major, minor, angle, unit = beams.pop()
        if unit != "rad":
            raise SystemExit(f"{image}'s restoring beam is in {unit}, not rad")
        header["BMAJ"] = math.degrees(major)
        header["BMIN"] = math.degrees(minor)
        header["BPA"] = math.degrees(angle)
    header["ORIGIN"] = "carta-zarr tools/testset/zarr-to-fits.py"
    header["HISTORY"] = f"pixels of {root.name}/{image}"
    return header


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset", help="an XRADIO image dataset")
    parser.add_argument("output", help="the FITS file to write")
    parser.add_argument("--block-mib", type=int, default=1024, help="how much to read and write at a time")
    parser.add_argument("--bench", default=os.environ.get("CARTA_ZARR_BENCH"), help="the carta-zarr-bench executable")
    args = parser.parse_args()
    if not args.bench:
        parser.error("--bench or $CARTA_ZARR_BENCH must name carta-zarr-bench")

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
    described = describe(args.bench, root)
    if described["id"] != image:
        raise SystemExit(f"carta-zarr opens {described['id']} by default, not {image}")
    _, n_freq, n_pol, n_l, n_m = sky.shape
    channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))

    out = Path(args.output)
    partial = out.with_name(out.name + ".partial")
    started = last_report = time.monotonic()
    with open(partial, "wb") as f:
        f.write(header_of(root, image, sky, described).tostring().encode("ascii"))
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
