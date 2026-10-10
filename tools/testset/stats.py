#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
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

"""The statistics a synthetic cube is calibrated against a real one with, as one JSON line.

  compression    decoded bytes over bytes on disk, of the chunks on disk, each counted whole
  nan            the median, over channels, of a plane's NaN fraction
  rms            the median, over channels, of a plane's robust rms (1.4826 x MAD)
  above, below   the median fraction of a plane's finite pixels beyond +5 and -5 robust rms
  nan_channels   how many channels are wholly NaN
  nan_chunks     how many chunks are wholly NaN (or not written at all)

    stats.py DATASET [--every N] [--block-mib 1024] [--bench PATH]

The pixels, the shape and the chunk are taken as carta-zarr reads them -- `carta-zarr-bench`, at
--bench or $CARTA_ZARR_BENCH -- so they are by axis name, of the chunk a read decodes, whatever order
the dimensions are stored in and whatever transpose precedes the shards; the pixel mask is applied.
Every Nth channel is sampled for the per-plane medians (8 by default); every channel for the NaN
counts. A statistic with nothing to take it of -- no chunk on disk, or no finite pixel in any plane
sampled -- is null.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


def compression(array: Path) -> float | None:
    """Decoded bytes over bytes on disk, over the chunks that are on disk. A writer leaves out a chunk
    that is all fill, and counting one as decoded bytes would credit the codec with what the writer
    saved. A chunk is counted whole, the padding of one at the image's edge included, because that is
    what the codec compressed and what a read decodes; in a shard, only the chunks its index lists."""
    metadata = json.loads((array / "zarr.json").read_text())
    itemsize = {"float32": 4, "float64": 8}[metadata["data_type"]]
    chunk = metadata["chunk_grid"]["configuration"]["chunk_shape"]
    sharding = next((codec["configuration"] for codec in metadata.get("codecs", [])
                     if codec.get("name") == "sharding_indexed"), None)
    decoded = stored = 0
    # Every file of the array but its metadata is a chunk or a shard, however chunk_key_encoding names
    # them: c/0/0, c.0.0, 0.0 or 0/0.
    for entry in array.rglob("*"):
        if not entry.is_file() or entry == array / "zarr.json" or entry.name.startswith("."):
            continue
        stored += entry.stat().st_size
        if sharding is None:
            decoded += itemsize * math.prod(chunk)
            continue
        inner = sharding["chunk_shape"]
        # A transpose ahead of the shards gives the inner chunk in its own axis order, so the count is
        # taken of the volumes, which no order changes.
        count = math.prod(chunk) // math.prod(inner)
        # The index is an (offset, length) pair of little-endian uint64 a chunk, then a crc32c when
        # its codecs say so; an absent chunk is all ones.
        checksum = any(codec.get("name") == "crc32c" for codec in sharding.get("index_codecs", []))
        size = 16 * count + (4 if checksum else 0)
        with entry.open("rb") as shard:
            if sharding.get("index_location", "end") == "end":
                shard.seek(-size, 2)
            index = np.frombuffer(shard.read(16 * count), dtype="<u8").reshape(count, 2)
        present = int((index[:, 0] != np.iinfo(np.uint64).max).sum())
        decoded += present * itemsize * math.prod(inner)
    # A dataset whose every chunk is fill has nothing on disk to compress.
    return decoded / stored if stored else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset")
    parser.add_argument("--every", type=int, default=8)
    parser.add_argument("--block-mib", type=int, default=1024, help="how much to read at a time")
    parser.add_argument("--bench", default=os.environ.get("CARTA_ZARR_BENCH"), help="the carta-zarr-bench executable")
    args = parser.parse_args()
    if not args.bench:
        parser.error("--bench or $CARTA_ZARR_BENCH must name carta-zarr-bench")

    root = Path(args.dataset)
    probe = subprocess.run([args.bench, "probe", str(root)], capture_output=True, text=True)
    report = json.loads(probe.stdout) if probe.stdout.strip() else {}
    if not report.get("ok"):
        raise SystemExit(f"carta-zarr does not open {root}: {report.get('error') or probe.stderr.strip()}")
    axes = {axis["name"]: axis for axis in report["image"]["axes"]}
    n_l, n_m, n_freq, n_pol = (axes[name]["length"] for name in ("l", "m", "frequency", "polarization"))
    chunk_l, chunk_m, depth = (axes[name]["chunk"] for name in ("l", "m", "frequency"))

    def first_polarization(start: int, stop: int) -> np.ndarray:
        """Channels start to stop of the first polarization through carta-zarr, as frequency, l, m."""
        result = subprocess.run([args.bench, "pixels", str(root), "--channels", f"{start}:{stop}"], capture_output=True)
        if result.returncode:
            raise SystemExit(result.stderr.decode(errors="replace").strip())
        return np.frombuffer(result.stdout, dtype=np.float32).reshape(n_pol, stop - start, n_m, n_l)[0].transpose(0, 2, 1)

    # Read a bounded number of planes at a time, whatever the chunk depth, and keep for each chunk
    # only whether any of its pixels so far is finite.
    planes = max(1, (args.block_mib << 20) // (n_pol * n_l * n_m * 4))
    grid = (math.ceil(n_l / chunk_l), math.ceil(n_m / chunk_m))
    padded = (grid[0] * chunk_l, grid[1] * chunk_m)
    chunk_finite = np.zeros((math.ceil(n_freq / depth), *grid), dtype=bool)
    nan_fraction, rms, above, below = [], [], [], []
    channel_nan = np.zeros(n_freq, dtype=bool)
    for start in range(0, n_freq, planes):
        block = first_polarization(start, min(n_freq, start + planes))
        for offset, plane in enumerate(block):
            f = start + offset
            finite = np.isfinite(plane)
            channel_nan[f] = not finite.any()
            whole = np.zeros(padded, dtype=bool)
            whole[:n_l, :n_m] = finite
            chunk_finite[f // depth] |= whole.reshape(grid[0], chunk_l, grid[1], chunk_m).any(axis=(1, 3))
            if f % args.every:
                continue
            values = plane[finite]
            nan_fraction.append(1 - values.size / plane.size)
            if values.size == 0:
                continue
            median = np.median(values)
            sigma = 1.4826 * np.median(np.abs(values - median))
            rms.append(sigma)
            above.append(np.mean(values > median + 5 * sigma))
            below.append(np.mean(values < median - 5 * sigma))
    nan_chunks = int((~chunk_finite).sum())

    def median(values: list[float], digits: int) -> float | None:
        """The median, or None when no plane sampled had a finite pixel to take one of."""
        return round(float(np.median(values)), digits) if values else None

    ratio = compression(root / "SKY")
    print(json.dumps({
        "dataset": root.name,
        "shape": [n_l, n_m, n_freq],
        "chunk": [chunk_l, chunk_m, depth],
        "compression": None if ratio is None else round(ratio, 3),
        "nan": median(nan_fraction, 4),
        "rms": float(f"{np.median(rms):.4g}") if rms else None,
        "above": median(above, 5),
        "below": median(below, 5),
        "nan_channels": int(channel_nan.sum()),
        "nan_chunks": nan_chunks,
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
