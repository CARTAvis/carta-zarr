#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "xradio[zarr]==1.3.0",
# ]
#
# [tool.uv]
# python-preference = "only-managed"
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convert a FITS cube to an XRADIO image Zarr with XRADIO's own open_image and write_image.

write_image writes XRADIO image schema 0.0.2 in the chunks, shards and compressor asked for, and
--compute_mask gives it a flag, true where a pixel is NaN, which the image declares. The FITS file is
read in whole rows of l, and along every other axis one output chunk deep, or one shard deep when
sharded: a FITS file is row-major, so narrower blocks read it in short strided pieces, and blocks of
whole chunks are written without being rechunked.

    convert.py CUBE.fits --chunks T,F,P,L,M [--shards T,F,P,L,M] [--compute_mask]
               [--compressor blosc-zstd|zstd] --output CUBE.zarr

--chunks and --shards are lengths along time, frequency, polarization, l and m, 0 for the whole axis.
--compressor is blosc with zstd at level 5 and byte shuffle, XRADIO's default, or zstd at level 5.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from xradio.image import open_image, write_image
from zarr.codecs import BloscCodec, ZstdCodec

DIMS = ("time", "frequency", "polarization", "l", "m")


def lengths(text: str) -> dict[str, int]:
    values = [int(value) for value in text.split(",")]
    if len(values) != len(DIMS) or any(value < 0 for value in values):
        raise SystemExit(f"expected {len(DIMS)} lengths of at least 0 along {', '.join(DIMS)}, got {text!r}")
    # write_image takes -1 for the whole axis.
    return {dim: value or -1 for dim, value in zip(DIMS, values)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a FITS cube to an XRADIO image Zarr.")
    parser.add_argument("fits")
    parser.add_argument("--chunks", required=True)
    parser.add_argument("--shards")
    parser.add_argument("--compute_mask", action="store_true")
    parser.add_argument("--compressor", choices=("blosc-zstd", "zstd"), default="blosc-zstd")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    chunks = lengths(args.chunks)
    shards = lengths(args.shards) if args.shards else None
    written = shards or chunks
    # open_image reads a dimension it is not given whole, which is what l needs.
    read = {dim: length for dim, length in written.items() if dim != "l" and length > 0}
    xds = open_image(args.fits, chunks=read, compute_mask=args.compute_mask)

    if args.compressor == "blosc-zstd":
        compressor = BloscCodec(cname="zstd", clevel=5, shuffle="shuffle")
    else:
        compressor = ZstdCodec(level=5)
    output = Path(args.output)
    # write_image names a store *.img.zarr, so it is written under that name and moved to the one
    # asked for.
    store = output.with_name(output.name + ".img.zarr")
    shutil.rmtree(store, ignore_errors=True)
    [written_store] = write_image(
        xds, str(store), out_format="zarr", chunks=chunks, shards=shards, compressor=compressor
    )
    shutil.rmtree(output, ignore_errors=True)
    Path(written_store).rename(output)


if __name__ == "__main__":
    main()
