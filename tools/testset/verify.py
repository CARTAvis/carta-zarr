#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "astropy==7.1.0",
#   "numpy==2.3.1",
# ]
#
# [tool.uv]
# python-preference = "only-managed"
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Check that a Zarr the converter made from a FITS cube is the layout asked for, holding its pixels.

A comparison between FITS and Zarr is only a comparison of the formats if both hold the same pixels,
in the layout the comparison is named after. So this asks carta-zarr, through carta-zarr-bench, as
carta-backend would read it: that the dataset opens, as SKY, with the FITS cube's shape; that it names
XRADIO image schema 0.0.2, the only one carta-zarr reads; that it is in the chunks (and shards) asked
for; that it carries a flag carta-zarr applies, or none, as --flag says; and that every pixel, as
carta-zarr decodes it, is the FITS cube's bit for bit -- with the flag applied, NaN exactly where the
FITS cube is.

What the coordinates, the unit or the beams mean is not compared: the test set exists to measure
reads, and whether carta-zarr reads metadata right is for carta-zarr's own tests to say. The flag is,
because the converter's --compute_mask leaves it undeclared (xradio 1.2.2 to 1.2.5), and a layout
that looks flagged and is not would measure something else.

    verify.py CUBE.fits CUBE.zarr --flag {yes,no} --chunks L,M,F [--shards L,M,F]
              [--bench PATH] [--block-mib 1024]

--bench is the carta-zarr-bench executable, $CARTA_ZARR_BENCH by default. Prints what it compared, and
exits nonzero if anything differs.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
from astropy.io import fits


def first_difference(differ: np.ndarray) -> tuple[int, ...] | None:
    """Where `differ` is first true, without listing everywhere it is."""
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

    probe = subprocess.run([args.bench, "probe", args.zarr], capture_output=True, text=True)
    try:
        report = json.loads(probe.stdout)
    except json.JSONDecodeError:
        report = {"ok": False, "error": (probe.stderr or probe.stdout).strip()}
    image = report.get("image", {})
    check("carta-zarr opens it" + ("" if report.get("ok") else f" ({report.get('error', 'no reason given')})"),
          bool(report.get("ok")) and image.get("id") == "SKY" and report.get("default_image_id") == "SKY")
    if failures:
        return finish()
    check("XRADIO image schema 0.0.2", report.get("schema", {}).get("version") == "0.0.2")
    axes = {axis["name"]: axis for axis in image["axes"]}

    with fits.open(args.fits, memmap=True) as hdul:
        data = hdul[0].data  # frequency, stokes, m, l
        n_freq, n_pol, n_m, n_l = data.shape
        check("shape l, m, frequency, polarization, time, as the FITS cube's",
              [axis["name"] for axis in image["axes"]] == ["l", "m", "frequency", "polarization", "time"]
              and [axis["length"] for axis in image["axes"]] == [n_l, n_m, n_freq, n_pol, 1])
        if failures:
            return finish()

        l, m, f = (int(part) for part in args.chunks.split(","))
        chunk = {"l": l, "m": m, "frequency": f, "polarization": 1, "time": 1}
        shard = dict(chunk)
        if args.shards:
            sl, sm, sf = (int(part) for part in args.shards.split(","))
            shard.update(l=sl, m=sm, frequency=sf)
        check(f"chunk {args.chunks}" + (f", shard {args.shards}" if args.shards else ", unsharded"),
              all(axes[name]["chunk"] == chunk[name] and axes[name]["shard"] == shard[name] for name in chunk)
              and image["sharded"] == bool(args.shards))
        flagged = args.flag == "yes"
        check("a flag carta-zarr applies" if flagged else "no flag carta-zarr applies", image["has_pixel_mask"] == flagged)

        def read(start: int, stop: int, masked: bool) -> np.ndarray:
            """Channels start to stop through carta-zarr, as frequency, polarization, m, l."""
            command = [args.bench, "pixels", args.zarr, "--channels", f"{start}:{stop}"] + ([] if masked else ["--unmasked"])
            result = subprocess.run(command, capture_output=True)
            if result.returncode:
                # TensorStore's own source locations and spec follow the reason, at length.
                raise RuntimeError(result.stderr.decode(errors="replace").strip().split(" [source locations=")[0])
            return np.frombuffer(result.stdout, dtype=np.float32).reshape(n_pol, stop - start, n_m, n_l).transpose(1, 0, 2, 3)

        channels = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))
        first_pixel = first_masked = unreadable = None
        for start in range(0, n_freq, channels):
            stop = min(n_freq, start + channels)
            expected = np.asarray(data[start:stop], dtype=np.float32)
            try:
                stored = read(start, stop, masked=False)
                applied = read(start, stop, masked=True) if flagged else None
            except RuntimeError as error:
                unreadable = f"channels {start}:{stop}: {error}"
                break
            if first_pixel is None and (where := first_difference(expected.view(np.uint32) != stored.view(np.uint32))):
                first_pixel = (start + where[0], *where[1:])
            if applied is not None and first_masked is None:
                # NaN where the FITS cube is NaN, and its pixel, bit for bit, everywhere else.
                nan = np.isnan(applied)
                wrong = (nan != np.isnan(expected)) | (~nan & (applied.view(np.uint32) != expected.view(np.uint32)))
                if where := first_difference(wrong):
                    first_masked = (start + where[0], *where[1:])
        check("carta-zarr reads every pixel" + (f" ({unreadable})" if unreadable else ""), unreadable is None)
        if unreadable is None:
            check("every pixel" + (f" (first differs at channel, stokes, y, x {first_pixel})" if first_pixel else ""),
                  first_pixel is None)
            if flagged:
                check("flagged where NaN" + (f" (first differs at {first_masked})" if first_masked else ""),
                      first_masked is None)

    return finish()


if __name__ == "__main__":
    sys.exit(main())
