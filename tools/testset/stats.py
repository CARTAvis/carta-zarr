#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
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

"""The statistics a synthetic cube is calibrated against a real one with, as one JSON line.

  compression    decoded bytes over bytes on disk, of the chunks on disk
  nan            the median, over channels, of a plane's NaN fraction
  rms            the median, over channels, of a plane's robust rms (1.4826 x MAD)
  above, below   the median fraction of a plane's finite pixels beyond +5 and -5 robust rms
  nan_channels   how many channels are wholly NaN
  nan_chunks     how many chunks are wholly NaN (or not written at all)

    stats.py DATASET [--every N] [--block-mib 1024]

Every Nth channel is sampled for the per-plane medians (8 by default); every channel for the NaN
counts.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import zarr


def compression(array: Path) -> float:
    """Decoded bytes over bytes on disk, over the chunks (or shards) that are on disk. A writer leaves
    out a chunk that is all fill, and counting one as decoded bytes would credit the codec with what
    the writer saved."""
    metadata = json.loads((array / "zarr.json").read_text())
    shape = metadata["shape"]
    outer = metadata["chunk_grid"]["configuration"]["chunk_shape"]
    itemsize = {"float32": 4, "float64": 8}[metadata["data_type"]]
    separator = metadata.get("chunk_key_encoding", {}).get("configuration", {}).get("separator", "/")
    decoded = stored = 0
    for entry in (array / "c").rglob("*"):
        if not entry.is_file():
            continue
        relative = entry.relative_to(array / "c")
        index = [int(part) for part in (relative.parts if separator == "/" else relative.name.split(separator))]
        decoded += itemsize * math.prod(min(extent, length - i * extent) for i, extent, length in zip(index, outer, shape))
        stored += entry.stat().st_size
    return decoded / stored


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset")
    parser.add_argument("--every", type=int, default=8)
    parser.add_argument("--block-mib", type=int, default=1024, help="how much to read at a time")
    args = parser.parse_args()

    root = Path(args.dataset)
    sky = zarr.open_array(str(root / "SKY"), mode="r")
    _, n_freq, _, n_l, n_m = sky.shape
    chunk = sky.chunks  # the inner chunk, sharded or not
    depth = chunk[1]

    # Read a bounded number of planes at a time, whatever the chunk depth, and keep for each chunk
    # only whether any of its pixels so far is finite.
    planes = max(1, (args.block_mib << 20) // (n_l * n_m * 4))
    grid = (math.ceil(n_l / chunk[3]), math.ceil(n_m / chunk[4]))
    padded = (grid[0] * chunk[3], grid[1] * chunk[4])
    chunk_finite = np.zeros((math.ceil(n_freq / depth), *grid), dtype=bool)
    nan_fraction, rms, above, below = [], [], [], []
    channel_nan = np.zeros(n_freq, dtype=bool)
    for start in range(0, n_freq, planes):
        block = np.asarray(sky[0, start : min(n_freq, start + planes), 0])  # frequency, l, m
        for offset, plane in enumerate(block):
            f = start + offset
            finite = np.isfinite(plane)
            channel_nan[f] = not finite.any()
            whole = np.zeros(padded, dtype=bool)
            whole[:n_l, :n_m] = finite
            chunk_finite[f // depth] |= whole.reshape(grid[0], chunk[3], grid[1], chunk[4]).any(axis=(1, 3))
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

    print(json.dumps({
        "dataset": root.name,
        "shape": [n_l, n_m, n_freq],
        "chunk": [chunk[3], chunk[4], depth],
        "compression": round(compression(root / "SKY"), 3),
        "nan": round(float(np.median(nan_fraction)), 4),
        "rms": float(f"{np.median(rms):.4g}"),
        "above": round(float(np.median(above)), 5),
        "below": round(float(np.median(below)), 5),
        "nan_channels": int(channel_nan.sum()),
        "nan_chunks": nan_chunks,
    }))
    return 0


if __name__ == "__main__":
    sys.exit(main())
