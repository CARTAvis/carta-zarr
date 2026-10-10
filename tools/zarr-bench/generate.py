#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "numpy==2.3.1",
#   "zarr==3.2.1",
# ]
#
# [tool.uv]
# # Never the system interpreter. Ubuntu 24.04's Python 3.12.3 segfaults the moment zarr starts its
# # event-loop thread, where uv's own 3.12 build does not -- and a measurement tool should not
# # depend on which Python a server happens to ship.
# python-preference = "only-managed"
# ///

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Write one XRADIO image dataset in a chosen storage layout, for measuring that layout.

A storage system answers "which chunk shape, which shard shape, which compressor" differently from
the disk these were first measured on, and the only way to know is to read the same pixels laid out
each way and time it. This writes those copies. carta-zarr-bench reads them; sweep.py drives both.

Two sources of pixels:

  --source PATH    Rewrite an existing XRADIO image dataset. The image (and its flag, when it has
                   one) gets the new layout; every other array is copied as it was. --crop takes a
                   window of it, which is how a cube too large to rewrite many times becomes one
                   that is not. This is the source to prefer: a layout's cost depends on how well
                   the pixels compress, and real pixels are the only ones that compress like real
                   pixels.

  --synthetic      A continuum-subtracted HI cube from a mosaic, like ASKAP's, by default: noise
                   that varies by channel and rises towards the footprint's edge, faint line
                   sources, NaN outside the footprint and runs of flagged channels; continuum
                   sources on request. The coordinates are built from a template XRADIO wrote (the
                   conformance fixture by default), stretched to --shape, on ASKAP's frequency axis
                   unless --frequency-start and --channel-width say otherwise.
                   Every pixel is a function of its position and --seed alone, so two layouts of
                   the same synthetic cube hold the same pixels, and only the layout differs
                   between them.

The layout is named per axis, in XRADIO's axis names, and an axis left out is 1:

  --chunk l=512,m=512,frequency=16 --shard frequency=256 --codec blosc:zstd:5:shuffle

Each dataset carries bench-manifest.json beside its zarr.json: what it was made from, the layout,
the striping it got, and what it came to on disk. The manifest is written last, so a dataset with
one is complete. With --output-root the directory is named after a hash of everything that decides
its bytes, and an existing complete one is reused rather than rewritten, unless --force asks for it
to be written again.

The last line on stdout is the dataset's path. Everything else goes to stderr.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import dataclasses
import datetime
import functools
import hashlib
import json
import math
import multiprocessing
import os
import platform
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import zarr

from fingerprint import reached, source_content

MANIFEST_NAME = "bench-manifest.json"
# Bumped whenever the same arguments would produce different bytes, so that a dataset written by an
# older generator is not mistaken for one this one would write.
FORMAT_VERSION = 4

AXES = ("time", "frequency", "polarization", "l", "m")
STOKES = ("I", "Q", "U", "V")
SPEED_OF_LIGHT = 299_792_458.0

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TEMPLATE = REPO_ROOT / "tests" / "data" / "images" / "zarr" / "xradio" / "conformance"

# The square of noise one seed draws at a time. Fixed, and independent of any layout, because it is
# what makes a synthetic pixel depend on where it is rather than on which block wrote it.
NOISE_TILE = 128

# A synthetic cube's defaults are an ASKAP HI cube's: its first channel, its channel width, HI's rest
# frequency, and the noise, in Jy/beam, that with the rise towards the footprint's edge comes to
# ASKAP Hydra's median plane rms of 2.1 mJy/beam.
# Where ASKAP is, as the site's converter writes it from a FITS file's OBSGEO-X, Y and Z: a spectral
# frame is placed by it.
ASKAP_TELESCOPE = {
    "name": "ASKAP",
    "direction": {"attrs": {"coordinate_system": "geocentric", "frame": "ITRF", "origin_object_name": "earth",
                            "type": "location", "units": "rad"},
                  "data": [2.0360801614255952, -0.46338342174681973], "dims": ["ellipsoid_dir_label"],
                  "coords": {"ellipsoid_dir_label": {"dims": ["ellipsoid_dir_label"], "data": ["lon", "lat"]}}},
    "distance": {"attrs": {"coordinate_system": "geocentric", "frame": "ITRF", "origin_object_name": "earth",
                           "type": "location", "units": "m"},
                 "data": [6373972.330145822], "dims": ["ellipsoid_dis_label"],
                 "coords": {"ellipsoid_dis_label": {"dims": ["ellipsoid_dis_label"], "data": ["dist"]}}},
}
ASKAP_FREQUENCY_START = 1.2955e9
ASKAP_CHANNEL_WIDTH = 18518.518518518518
HI_REST_FREQUENCY = 1.420405751786e9
ASKAP_NOISE = 1.7e-3
# One line source for every this many pixels, when --line-sources does not say.
PIXELS_PER_LINE_SOURCE = 10**8

# The shape of a synthetic cube, chosen so that its statistics match ASKAP Hydra's (see
# tools/testset/README.md): the footprint's corners (a superellipse of this power), how much the noise
# rises at its edge, how many periods of ripple the channel rms has across the band, and the fraction
# and width of the noise drawn from a wider Gaussian.
FOOTPRINT_POWER = 6
EDGE_GAIN = 1.0
RIPPLE_PERIODS = 3
TAIL_FRACTION = 0.044
TAIL_SCALE = 6.0

# Seeds of the tables a synthetic cube consults, kept apart from the noise's, which is seeded by
# position.
STREAM_CHANNELS = 0xC4A1
STREAM_FOOTPRINT = 0xF007
STREAM_FLAGGED = 0xF1A6
STREAM_LINES = 0x11E5
STREAM_CONTINUUM = 0x5EED


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


# -- Arguments ------------------------------------------------------------------------------------


def parse_axis_map(text: str | None, what: str) -> dict[str, int]:
    if not text:
        return {}
    result: dict[str, int] = {}
    for item in text.split(","):
        name, _, value = item.partition("=")
        name = name.strip()
        if name not in AXES:
            raise SystemExit(f"{what}: unknown axis {name!r}; the axes are {', '.join(AXES)}")
        try:
            result[name] = int(value)
        except ValueError:
            raise SystemExit(f"{what}: {item!r} is not axis=integer") from None
        if result[name] <= 0:
            raise SystemExit(f"{what}: {name} must be positive")
    return result


def parse_crop(text: str | None) -> dict[str, tuple[int, int]]:
    if not text:
        return {}
    result: dict[str, tuple[int, int]] = {}
    for item in text.split(","):
        name, _, window = item.partition("=")
        name = name.strip()
        if name not in AXES:
            raise SystemExit(f"--crop: unknown axis {name!r}")
        start, sep, stop = window.partition(":")
        if not sep:
            raise SystemExit(f"--crop: {item!r} is not axis=start:stop")
        result[name] = (int(start), int(stop))
    return result


def parse_size(text: str) -> int:
    """A byte count with an optional binary suffix: 4M, 1MiB, 512k, 1048576."""
    value = text.strip().lower().removesuffix("ib").removesuffix("b")
    scale = {"k": 1 << 10, "m": 1 << 20, "g": 1 << 30}.get(value[-1:], 1)
    digits = value[:-1] if scale > 1 else value
    try:
        return int(digits) * scale
    except ValueError:
        raise SystemExit(f"{text!r} is not a size") from None


@dataclasses.dataclass(frozen=True)
class Codec:
    """A compressor, as a canonical name and the codec JSON zarr v3 stores for it."""

    name: str
    level: int | None = None
    blosc_cname: str | None = None
    blosc_shuffle: str | None = None

    @staticmethod
    def parse(text: str) -> Codec:
        parts = text.split(":")
        kind = parts[0]
        if kind == "none" and len(parts) == 1:
            return Codec("none")
        if kind in ("zstd", "gzip") and len(parts) <= 2:
            default = 3 if kind == "zstd" else 6
            return Codec(kind, int(parts[1]) if len(parts) == 2 else default)
        if kind == "blosc" and len(parts) <= 4:
            cname = parts[1] if len(parts) > 1 else "zstd"
            level = int(parts[2]) if len(parts) > 2 else 5
            shuffle = parts[3] if len(parts) > 3 else "shuffle"
            if shuffle not in ("noshuffle", "shuffle", "bitshuffle"):
                raise SystemExit(f"--codec: blosc shuffle must be noshuffle, shuffle or bitshuffle, not {shuffle!r}")
            return Codec("blosc", level, cname, shuffle)
        raise SystemExit(
            f"--codec: {text!r} is not none, zstd[:level], gzip[:level] or blosc[:cname[:level[:shuffle]]]"
        )

    def spelling(self) -> str:
        if self.name == "none":
            return "none"
        if self.name == "blosc":
            return f"blosc:{self.blosc_cname}:{self.level}:{self.blosc_shuffle}"
        return f"{self.name}:{self.level}"

    def json(self, itemsize: int) -> list[dict[str, Any]]:
        if self.name == "none":
            return []
        if self.name == "zstd":
            return [{"name": "zstd", "configuration": {"level": self.level, "checksum": False}}]
        if self.name == "gzip":
            return [{"name": "gzip", "configuration": {"level": self.level}}]
        return [
            {
                "name": "blosc",
                "configuration": {
                    "cname": self.blosc_cname,
                    "clevel": self.level,
                    "shuffle": self.blosc_shuffle,
                    "typesize": itemsize,
                    "blocksize": 0,
                },
            }
        ]


@dataclasses.dataclass(frozen=True)
class Stripe:
    """Striping to set on the output directory before anything is written into it."""

    filesystem: str  # lustre or beegfs
    count: int
    size: int  # bytes

    @staticmethod
    def parse(text: str) -> Stripe:
        filesystem, _, rest = text.partition(":")
        if filesystem not in ("lustre", "beegfs"):
            raise SystemExit(f"--stripe: {filesystem!r} is not lustre or beegfs")
        settings = dict(item.split("=", 1) for item in rest.split(",") if item)
        unknown = set(settings) - {"count", "size"}
        if unknown or len(settings) != 2:
            raise SystemExit(f"--stripe: {filesystem} takes count= and size= and nothing else")
        return Stripe(filesystem, int(settings["count"]), parse_size(settings["size"]))

    def spelling(self) -> str:
        return f"{self.filesystem}:count={self.count},size={self.size}"

    def size_with_unit(self, kibi: str, mebi: str) -> str:
        """The size in the largest binary unit it is a whole number of, spelt as the tool wants."""
        if self.size % (1 << 20) == 0:
            return f"{self.size >> 20}{mebi}"
        if self.size % (1 << 10) == 0:
            return f"{self.size >> 10}{kibi}"
        return str(self.size)


# -- Layout ---------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Layout:
    chunk: tuple[int, ...]  # one per stored axis, in the image's dimension order
    shard: tuple[int, ...] | None
    codec: Codec
    keep_bits: int | None

    def outer(self) -> tuple[int, ...]:
        """The unit one file holds: the shard when there is one, the chunk otherwise."""
        return self.shard or self.chunk


def resolve_layout(args: argparse.Namespace, dims: list[str], shape: list[int]) -> Layout:
    chunk_request = parse_axis_map(args.chunk, "--chunk")
    shard_request = parse_axis_map(args.shard, "--shard")
    for name in list(chunk_request) + list(shard_request):
        if name not in dims:
            raise SystemExit(f"the image has no {name} axis; its axes are {', '.join(dims)}")

    # A chunk longer than its axis buys nothing but padding, so it is clipped to the axis.
    chunk = tuple(min(chunk_request.get(name, 1), length) for name, length in zip(dims, shape))
    shard = None
    if args.shard:
        outer = []
        for name, length, inner in zip(dims, shape, chunk):
            requested = shard_request.get(name, inner)
            if requested % inner:
                raise SystemExit(f"--shard: {name}={requested} is not a multiple of its chunk, {inner}")
            # Clipped to the whole chunks that cover the axis, which keeps it a multiple of one.
            outer.append(min(requested, math.ceil(length / inner) * inner))
        shard = tuple(outer)
    keep_bits = args.keep_bits
    if keep_bits is not None and not 1 <= keep_bits <= 23:
        raise SystemExit("--keep-bits: a float32 mantissa keeps between 1 and 23 bits")
    return Layout(chunk, shard, Codec.parse(args.codec), keep_bits)


def image_codecs(layout: Layout, itemsize: int) -> list[dict[str, Any]]:
    serializer = {"name": "bytes", "configuration": {"endian": "little"}} if itemsize > 1 else {"name": "bytes"}
    inner = [serializer, *layout.codec.json(itemsize)]
    if layout.shard is None:
        return inner
    return [
        {
            "name": "sharding_indexed",
            "configuration": {
                "chunk_shape": list(layout.chunk),
                "codecs": inner,
                "index_codecs": [{"name": "bytes", "configuration": {"endian": "little"}}, {"name": "crc32c"}],
                "index_location": "end",
            },
        }
    ]


def with_layout(metadata: dict[str, Any], shape: list[int], layout: Layout, itemsize: int) -> dict[str, Any]:
    result = json.loads(json.dumps(metadata))
    result["shape"] = list(shape)
    result["chunk_grid"] = {"name": "regular", "configuration": {"chunk_shape": list(layout.outer())}}
    result["codecs"] = image_codecs(layout, itemsize)
    result.pop("storage_transformers", None)
    return result


def write_array_metadata(path: Path, metadata: dict[str, Any]) -> None:
    path.mkdir(parents=True, exist_ok=False)
    (path / "zarr.json").write_text(json.dumps(metadata, indent=2) + "\n")


def read_metadata(path: Path) -> dict[str, Any]:
    return json.loads((path / "zarr.json").read_text())


# -- Pixels ---------------------------------------------------------------------------------------


def round_mantissa(values: np.ndarray, keep_bits: int) -> np.ndarray:
    """Round float32 values to `keep_bits` mantissa bits, ties to even, leaving NaN and inf alone.

    The trailing bits become zero, which is what lets a compressor find something to compress in
    data that is otherwise noise. It is what numcodecs' BitRound does, and what a pipeline that
    quantises its products would hand a reader.
    """
    if keep_bits >= 23:
        return values
    bits = values.view(np.uint32)
    drop = 23 - keep_bits
    half = np.uint32((1 << (drop - 1)) - 1)
    rounded = (bits + half + ((bits >> np.uint32(drop)) & np.uint32(1))) & np.uint32(~((1 << drop) - 1) & 0xFFFFFFFF)
    return np.where(np.isfinite(values), rounded.view(np.float32), values)


@dataclasses.dataclass(frozen=True)
class Synthetic:
    """A cube whose every pixel is a function of its position and the seed.

    By default it looks like a continuum-subtracted HI cube from a mosaic, such as ASKAP's: noise
    whose rms varies from channel to channel and rises towards the edge of the footprint, with heavy
    tails; faint line sources, each present in a few hundred km/s of the band; NaN outside an
    irregular footprint; and runs of channels flagged in part or in whole. Continuum sources, point
    and extended, can be added.

    The noise is drawn a NOISE_TILE square at a time from a generator seeded by that tile's position,
    and everything else is a fixed table each block consults, so the value of a pixel does not depend
    on the block that wrote it. That is what makes two layouts of one synthetic cube comparable.
    """

    shape: tuple[int, int, int, int, int]  # time, frequency, polarization, l, m
    frequencies: tuple[float, ...]
    seed: int
    noise: float
    line_sources: int
    point_sources: int
    extended_sources: int
    footprint_fill: float | None  # the fraction of a plane inside the footprint; None for no footprint
    flagged_channels: float  # the fraction of channels in flagged runs
    flag: bool
    flagged_range: tuple[int, int] | None = None  # channels flagged whole besides the runs, start and stop

    @functools.cached_property
    def channel_rms(self) -> np.ndarray:
        """Each channel's noise rms: a slow ripple across the band and a jitter from channel to channel,
        together about 15 % either way."""
        rng = np.random.default_rng([self.seed, STREAM_CHANNELS])
        count = self.shape[1]
        position = np.arange(count) / max(count, 1)
        ripple = 0.08 * np.sin(2 * math.pi * (RIPPLE_PERIODS * position + rng.uniform()))
        jitter = rng.normal(0.0, 0.08, count)
        return (self.noise * np.exp(ripple + jitter)).astype(np.float32)

    @functools.cached_property
    def footprint(self) -> dict[str, Any]:
        """A rounded rectangle whose radius wanders with angle -- a mosaic of beams, roughly -- scaled
        so that `footprint_fill` of the plane is inside it."""
        rng = np.random.default_rng([self.seed, STREAM_FOOTPRINT])
        orders = np.arange(3, 9)
        amplitudes = rng.uniform(0.0, 0.02, orders.size)
        phases = rng.uniform(0.0, 2 * math.pi, orders.size)
        # Inside |u|^p + |v|^p <= R^p is 4 R^2 G(1 + 1/p)^2 / G(1 + 2/p) of the plane's 4.
        p = FOOTPRINT_POWER
        radius = math.sqrt((self.footprint_fill or 1.0) * math.gamma(1 + 2 / p) / math.gamma(1 + 1 / p) ** 2)
        return {"orders": orders, "amplitudes": amplitudes, "phases": phases, "radius": radius}

    def within(self, l_range: range, m_range: range) -> np.ndarray | None:
        """How far into the footprint each pixel is: 0 at its centre, 1 at its edge, over 1 outside.
        None when there is no footprint."""
        if self.footprint_fill is None:
            return None
        _, _, _, length_l, length_m = self.shape
        u = (np.asarray(l_range, dtype=np.float64)[:, None] - (length_l - 1) / 2) / (length_l / 2)
        v = (np.asarray(m_range, dtype=np.float64)[None, :] - (length_m - 1) / 2) / (length_m / 2)
        shape = self.footprint
        angle = np.arctan2(v, u)
        wander = sum(a * np.cos(k * angle + phase)
                     for k, a, phase in zip(shape["orders"], shape["amplitudes"], shape["phases"]))
        p = FOOTPRINT_POWER
        return (np.abs(u) ** p + np.abs(v) ** p) ** (1 / p) / (shape["radius"] * (1 + wander))

    @functools.cached_property
    def flagged(self) -> dict[str, np.ndarray]:
        """Runs of flagged channels, of lengths spread evenly in log from 1 to 128 (or to the number of
        channels to flag, when that is fewer): half flagged whole, half on one side of a line across
        the plane. `flagged_range`, when there is one, is flagged whole on top of them: a run placed
        to cover whole chunks, which random runs of a few per cent of a short cube never do."""
        rng = np.random.default_rng([self.seed, STREAM_FLAGGED])
        count = self.shape[1]
        whole = np.zeros(count, dtype=bool)
        angle = np.full(count, np.nan)
        offset = np.zeros(count)
        target = round(self.flagged_channels * count)
        longest = max(1, min(128, target))
        while (whole | np.isfinite(angle)).sum() < target:
            length = min(int(math.exp(rng.uniform(0.0, math.log(longest + 1)))), longest)
            start = int(rng.integers(0, count - length + 1))
            if rng.uniform() < 0.5:
                whole[start : start + length] = True
            else:
                angle[start : start + length] = rng.uniform(0.0, 2 * math.pi)
                offset[start : start + length] = rng.uniform(-0.3, 0.3)
        if self.flagged_range is not None:
            whole[self.flagged_range[0] : self.flagged_range[1]] = True
        return {"whole": whole, "angle": angle, "offset": offset}

    def dropped(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        """Which pixels of a block are not data: outside the footprint, or in a flagged channel."""
        _, f0, _, l0, m0 = start
        _, f1, _, l1, m1 = stop
        dropped = np.zeros([b - a for a, b in zip(start, stop)], dtype=bool)
        within = self.within(range(l0, l1), range(m0, m1))
        if within is not None:
            dropped[..., within > 1] = True
        flagged = self.flagged
        _, _, _, length_l, length_m = self.shape
        u = (np.arange(l0, l1, dtype=np.float64)[:, None] - (length_l - 1) / 2) / (length_l / 2)
        v = (np.arange(m0, m1, dtype=np.float64)[None, :] - (length_m - 1) / 2) / (length_m / 2)
        for f in range(f0, f1):
            if flagged["whole"][f]:
                dropped[:, f - f0] = True
            elif math.isfinite(flagged["angle"][f]):
                side = u * math.cos(flagged["angle"][f]) + v * math.sin(flagged["angle"][f]) > flagged["offset"][f]
                dropped[:, f - f0][..., side] = True
        return dropped

    @functools.cached_property
    def lines(self) -> dict[str, np.ndarray]:
        """The line sources: an elliptical Gaussian on the sky with a Gaussian or double-horned profile
        50 to 500 km/s wide, its peak 1 to 20 times the noise."""
        rng = np.random.default_rng([self.seed, STREAM_LINES])
        n = self.line_sources
        _, count, _, length_l, length_m = self.shape
        centre = rng.uniform(-0.05 * count, 1.05 * count, n)
        width_kms = rng.uniform(50.0, 500.0, n)
        frequencies = np.asarray(self.frequencies)
        step = abs(frequencies[1] - frequencies[0]) if count > 1 else 1.0
        at = frequencies[np.clip(centre.astype(int), 0, count - 1)]
        width = width_kms * 1e3 / SPEED_OF_LIGHT * at / step  # in channels
        return {
            "l": rng.uniform(0, length_l, n),
            "m": rng.uniform(0, length_m, n),
            "major": np.exp(rng.uniform(math.log(1.5), math.log(12.0), n)),
            "ratio": rng.uniform(0.3, 1.0, n),
            "angle": rng.uniform(0.0, math.pi, n),
            "peak": self.noise * 10 ** rng.uniform(0.0, 1.3, n),
            "centre": centre,
            "width": width,
            "horn": np.where(rng.uniform(size=n) < 0.5, 0.0, rng.uniform(0.2, 0.5, n)),
        }

    @functools.cached_property
    def continuum(self) -> dict[str, np.ndarray]:
        """Point and extended continuum sources, each with a spectral index: elliptical Gaussians, the
        point ones a beam or two across, the extended ones tens of pixels."""
        rng = np.random.default_rng([self.seed, STREAM_CONTINUUM])
        _, _, _, length_l, length_m = self.shape
        points, extended = self.point_sources, self.extended_sources
        n = points + extended
        major = np.concatenate([rng.uniform(1.5, 4.0, points), np.exp(rng.uniform(math.log(20), math.log(100), extended))])
        return {
            "l": rng.uniform(0, length_l, n),
            "m": rng.uniform(0, length_m, n),
            "major": major,
            "ratio": np.concatenate([np.ones(points), rng.uniform(0.3, 1.0, extended)]),
            "angle": rng.uniform(0.0, math.pi, n),
            "peak": np.concatenate([10 ** rng.uniform(-3, 0, points), 10 ** rng.uniform(-3, -1.5, extended)]),
            "index": rng.normal(-0.7, 0.3, n),
        }

    def noise_block(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        block = np.empty([b - a for a, b in zip(start, stop)], dtype=np.float32)
        t0, f0, p0, l0, m0 = start
        t1, f1, p1, l1, m1 = stop
        tile = NOISE_TILE
        for t in range(t0, t1):
            for f in range(f0, f1):
                for p in range(p0, p1):
                    plane = block[t - t0, f - f0, p - p0]
                    for tl in range(l0 // tile, (l1 - 1) // tile + 1):
                        for tm in range(m0 // tile, (m1 - 1) // tile + 1):
                            rng = np.random.default_rng([self.seed, t, f, p, tl, tm])
                            noise = rng.standard_normal((tile, tile), dtype=np.float32)
                            # A few pixels from a wider distribution: the symmetric excess beyond
                            # five sigma that sidelobes and calibration errors leave in real cubes.
                            wide = rng.random((tile, tile), dtype=np.float32) < TAIL_FRACTION
                            noise[wide] *= np.float32(TAIL_SCALE)
                            la, lb = max(l0, tl * tile), min(l1, (tl + 1) * tile)
                            ma, mb = max(m0, tm * tile), min(m1, (tm + 1) * tile)
                            plane[la - l0 : lb - l0, ma - m0 : mb - m0] = noise[
                                la - tl * tile : lb - tl * tile, ma - tm * tile : mb - tm * tile
                            ]
        block *= self.channel_rms[f0:f1][None, :, None, None, None]
        # Primary-beam correction raises the noise towards the edge of the footprint.
        within = self.within(range(l0, l1), range(m0, m1))
        if within is not None:
            block *= (1 + EDGE_GAIN * np.minimum(within, 1.0) ** 8).astype(np.float32)
        return block

    def add_gaussians(self, block: np.ndarray, start: tuple[int, ...], stop: tuple[int, ...],
                      table: dict[str, np.ndarray], spectra) -> None:
        """Add elliptical Gaussians, each clipped to five of its sigma before it is clipped to the
        block, so that a source reaches the same pixels whatever the block. `spectra(index, f0, f1)`
        is the source's value in channels f0 to f1, or None where it has none."""
        _, f0, _, l0, m0 = start
        _, f1, _, l1, m1 = stop
        # Each source's own pixels, in whole pixels, decided before the block is looked at: a source
        # is in a block exactly when these meet it, so no block can leave out a pixel another adds.
        reach = 5 * table["major"]
        first_l, last_l = np.floor(table["l"] - reach), np.ceil(table["l"] + reach)
        first_m, last_m = np.floor(table["m"] - reach), np.ceil(table["m"] + reach)
        near = np.flatnonzero((last_l >= l0) & (first_l < l1) & (last_m >= m0) & (first_m < m1))
        for index in near:
            spectrum = spectra(index, f0, f1)
            if spectrum is None:
                continue
            centre_l, centre_m = table["l"][index], table["m"][index]
            la, lb = max(l0, int(first_l[index])), min(l1, int(last_l[index]) + 1)
            ma, mb = max(m0, int(first_m[index])), min(m1, int(last_m[index]) + 1)
            dl = np.arange(la, lb)[:, None] - centre_l
            dm = np.arange(ma, mb)[None, :] - centre_m
            cos, sin = math.cos(table["angle"][index]), math.sin(table["angle"][index])
            major = table["major"][index]
            minor = major * table["ratio"][index]
            x, y = dl * cos + dm * sin, -dl * sin + dm * cos
            sky = np.exp(-0.5 * ((x / major) ** 2 + (y / minor) ** 2))
            # Only the channels the source has: adding its zeros elsewhere would turn a -0 of the
            # noise into +0 in a block that holds some of its band and not in one that holds none.
            band = np.flatnonzero(spectrum)
            if band.size == 0:
                continue
            fa, fb = int(band[0]), int(band[-1]) + 1
            blob = spectrum[fa:fb, None, None] * sky[None, :, :]
            block[:, fa:fb, 0, la - l0 : lb - l0, ma - m0 : mb - m0] += blob.astype(np.float32)[None]

    def line_spectrum(self, index: int, f0: int, f1: int) -> np.ndarray | None:
        lines = self.lines
        centre, width, horn = lines["centre"][index], lines["width"][index], lines["horn"][index]
        if horn == 0:
            sigma = width / 2.355
            half = 5 * sigma
        else:
            edge = max(width / 10, 0.5)
            half = width / 2 + 5 * edge
        a, b = max(f0, math.floor(centre - half)), min(f1, math.ceil(centre + half) + 1)
        if a >= b:
            return None
        spectrum = np.zeros(f1 - f0)
        x = np.arange(a, b) - centre
        if horn == 0:
            profile = np.exp(-0.5 * (x / sigma) ** 2)
        else:
            # A boxcar with soft edges, dipped in the middle: the two horns of a rotating disc.
            erf = np.vectorize(math.erf)
            box = 0.5 * (erf((x + width / 2) / (math.sqrt(2) * edge)) - erf((x - width / 2) / (math.sqrt(2) * edge)))
            profile = box * (1 - horn + horn * np.minimum(1.0, (2 * x / width) ** 2))
        spectrum[a - f0 : b - f0] = lines["peak"][index] * profile
        return spectrum

    def continuum_spectrum(self, index: int, f0: int, f1: int) -> np.ndarray:
        base = np.asarray(self.frequencies[f0:f1], dtype=np.float64) / self.frequencies[0]
        continuum = self.continuum
        return continuum["peak"][index] * base ** continuum["index"][index]

    def pixels(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        block = self.noise_block(start, stop)
        # Sources are in Stokes I alone.
        if start[2] == 0:
            if self.line_sources:
                self.add_gaussians(block, start, stop, self.lines, self.line_spectrum)
            if self.point_sources or self.extended_sources:
                self.add_gaussians(block, start, stop, self.continuum, self.continuum_spectrum)
        if not self.flag:
            block[self.dropped(start, stop)] = np.nan
        return block

    def flags(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        # True is a flagged pixel: XRADIO's flag says which pixels to drop. See src/pixel_mask.h.
        return self.dropped(start, stop)


# -- Writing --------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class SkyDirection:
    """A synthetic cube's right ascension and declination: the SIN projection inverted, as XRADIO's
    are, over the cube's l and m. Each is as large as a plane -- 8 GiB apiece at 32768 square -- so
    it is computed a block at a time, as the pixels are, rather than whole."""

    reference: tuple[float, float]  # right ascension and declination of the projection centre
    lengths: tuple[int, int]  # of l and m
    increments: tuple[float, float]  # of l and m

    def axis(self, which: int, start: int, stop: int) -> np.ndarray:
        """l or m from `start` to `stop`, exactly as the cube's l and m arrays hold them."""
        return (np.arange(start, stop) - self.lengths[which] // 2) * self.increments[which]

    def block(self, kind: str, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        ra0, dec0 = self.reference
        l, m = self.axis(0, start[0], stop[0])[:, None], self.axis(1, start[1], stop[1])[None, :]
        n = np.sqrt(np.maximum(0.0, 1.0 - l * l - m * m))
        if kind == "declination":
            return np.arcsin(m * math.cos(dec0) + n * math.sin(dec0))
        return ra0 + np.arctan2(l, n * math.cos(dec0) - m * math.sin(dec0))


# The coordinates as large as a plane, written block by block through fill().
SKY_DIRECTION = ("right_ascension", "declination")


@dataclasses.dataclass(frozen=True)
class Job:
    """One array to fill, block by block: everything a worker process needs to do its share."""

    target: str  # path of the array being written
    kind: str  # "copy", "pixels", "flags", or one of SKY_DIRECTION
    source: str | None = None  # array read from, for a copy
    offset: tuple[int, ...] = ()  # where the target's origin sits in the source, for a crop
    synthetic: Synthetic | None = None
    keep_bits: int | None = None
    direction: SkyDirection | None = None


_OPEN: dict[str, Any] = {}


def _open(path: str, mode: str) -> Any:
    key = f"{mode}:{path}"
    if key not in _OPEN:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _OPEN[key] = zarr.open_array(path, mode=mode)
    return _OPEN[key]


_JOB: Job | None = None


def _take_job(job: Job) -> None:
    # Handed to each worker once rather than with every block: a synthetic job carries the whole
    # frequency axis, and a large cube is tens of thousands of blocks.
    global _JOB
    _JOB = job


def _fill_taken(start: tuple[int, ...], stop: tuple[int, ...]) -> int:
    return fill_block(_JOB, start, stop)


def fill_block(job: Job, start: tuple[int, ...], stop: tuple[int, ...]) -> int:
    window = tuple(slice(a, b) for a, b in zip(start, stop))
    if job.kind == "pixels":
        values = job.synthetic.pixels(start, stop)
    elif job.kind == "flags":
        values = job.synthetic.flags(start, stop)
    elif job.kind in SKY_DIRECTION:
        values = job.direction.block(job.kind, start, stop)
    else:
        shifted = tuple(slice(a + o, b + o) for a, b, o in zip(start, stop, job.offset))
        values = np.asarray(_open(job.source, "r")[shifted])
    if job.keep_bits is not None and values.dtype == np.float32:
        values = round_mantissa(values, job.keep_bits)
    _open(job.target, "r+")[window] = values
    return values.nbytes


def plan_blocks(shape: list[int], unit: tuple[int, ...], itemsize: int, target_bytes: int) -> list[tuple]:
    """Blocks of whole units -- whole shards, or whole chunks when there are none -- so that no two
    workers ever write the same file, grown last axis first until each is about `target_bytes`."""
    multiple = [1] * len(shape)
    units_per_axis = [math.ceil(length / extent) for length, extent in zip(shape, unit)]
    unit_bytes = math.prod(unit) * itemsize
    for axis in reversed(range(len(shape))):
        while unit_bytes * math.prod(multiple) < target_bytes and multiple[axis] < units_per_axis[axis]:
            multiple[axis] = min(units_per_axis[axis], multiple[axis] * 2)
    block = [extent * factor for extent, factor in zip(unit, multiple)]
    starts = [range(0, length, step) for length, step in zip(shape, block)]
    blocks = []
    for origin in np.ndindex(*[len(axis) for axis in starts]):
        start = tuple(axis[i] for axis, i in zip(starts, origin))
        stop = tuple(min(a + step, length) for a, step, length in zip(start, block, shape))
        blocks.append((start, stop))
    return blocks


def fill(job: Job, shape: list[int], unit: tuple[int, ...], itemsize: int, args: argparse.Namespace) -> None:
    blocks = plan_blocks(shape, unit, itemsize, args.block_mib << 20)
    total = math.prod(shape) * itemsize
    name = Path(job.target).name
    started = last_report = time.monotonic()
    written = 0
    if args.workers <= 1 or len(blocks) == 1:
        results = (fill_block(job, start, stop) for start, stop in blocks)
        pool = None
    else:
        # spawn rather than fork: zarr runs an event loop on a thread of its own, and a forked child
        # inherits the loop's state without the thread.
        pool = concurrent.futures.ProcessPoolExecutor(
            args.workers, mp_context=multiprocessing.get_context("spawn"), initializer=_take_job, initargs=(job,)
        )
        futures = [pool.submit(_fill_taken, start, stop) for start, stop in blocks]
        results = (future.result() for future in concurrent.futures.as_completed(futures))
    try:
        for count, nbytes in enumerate(results, start=1):
            written += nbytes
            now = time.monotonic()
            if now - last_report >= 5 or count == len(blocks):
                rate = written / max(now - started, 1e-9) / (1 << 20)
                log(f"  {name}: {count}/{len(blocks)} blocks, {written / (1 << 30):.2f} of "
                    f"{total / (1 << 30):.2f} GiB, {rate:.0f} MiB/s")
                last_report = now
    finally:
        if pool is not None:
            pool.shutdown(cancel_futures=True)


# -- Dataset assembly -----------------------------------------------------------------------------


def holds_booleans(metadata: dict[str, Any]) -> bool:
    """A Zarr bool, or xarray's encoding of one -- an int8 marked dtype "bool" -- which is what XRADIO
    writes for every flag. HoldsBooleans in src/schema/xradio/flag.h."""
    if metadata.get("data_type") == "bool":
        return True
    return metadata.get("data_type") == "int8" and metadata.get("attributes", {}).get("dtype") == "bool"


def flags_of_data_groups(root_metadata: dict[str, Any], image: str) -> tuple[set[str], set[str]]:
    """The flags the root's data groups give this image as its sky, and those they give any other."""
    own, others = set(), set()
    groups = root_metadata.get("attributes", {}).get("data_groups")
    for group in (groups.values() if isinstance(groups, dict) else ()):
        flag = group.get("flag") if isinstance(group, dict) else None
        if isinstance(flag, str) and flag:
            (own if group.get("sky") == image else others).add(flag)
    return own, others


def find_flag(root: Path, image: str, image_metadata: dict[str, Any], names: list[str],
              root_metadata: dict[str, Any]) -> str | None:
    """The flag carta-zarr would mask this image with, by the rule of src/schema/xradio/flag.cc: the
    one its own `flag` attribute names, else the one data_groups give it as its sky, else the only
    flag-typed boolean array with its dimensions and shape that no data group gives another image.
    Two groups giving it different flags is refused, as carta-zarr refuses to open it."""
    declared = image_metadata.get("attributes", {}).get("flag")
    if isinstance(declared, str) and declared:
        return declared
    own, others = flags_of_data_groups(root_metadata, image)
    if len(own) > 1:
        raise SystemExit(f"data_groups name more than one flag for {image}: {', '.join(sorted(own))}")
    if own:
        return next(iter(own))
    matches = []
    for name in names:
        if name == image or name in others:
            continue
        metadata = read_metadata(root / name)
        if (
            metadata.get("node_type") == "array"
            and metadata.get("attributes", {}).get("type") == "flag"
            and holds_booleans(metadata)
            and metadata.get("dimension_names") == image_metadata.get("dimension_names")
            and metadata.get("shape") == image_metadata.get("shape")
        ):
            matches.append(name)
    return matches[0] if len(matches) == 1 else None


def default_image(root_metadata: dict[str, Any]) -> str:
    groups = root_metadata.get("attributes", {}).get("data_groups", {})
    base = groups.get("base", {}) if isinstance(groups, dict) else {}
    return base.get("sky", "SKY") if isinstance(base, dict) else "SKY"


def member_names(root: Path) -> list[str]:
    return sorted(entry.name for entry in root.iterdir() if (entry / "zarr.json").is_file())


def rewrite_source(args: argparse.Namespace, out: Path) -> dict[str, Any]:
    source = Path(args.source).resolve()
    root_metadata = read_metadata(source)
    names = member_names(source)
    image = args.image or default_image(root_metadata)
    if image not in names:
        raise SystemExit(f"{source} has no array {image!r}")
    image_metadata = read_metadata(source / image)
    dims = image_metadata.get("dimension_names") or []
    crop = parse_crop(args.crop)
    for name, (start, stop) in crop.items():
        if name not in dims:
            raise SystemExit(f"--crop: the image has no {name} axis")
        length = image_metadata["shape"][dims.index(name)]
        if not 0 <= start < stop <= length:
            raise SystemExit(f"--crop: {name}={start}:{stop} is not inside 0:{length}")

    def cropped(metadata: dict[str, Any]) -> tuple[list[int], tuple[int, ...]]:
        array_dims = metadata.get("dimension_names") or [None] * len(metadata["shape"])
        shape, offset = [], []
        for name, length in zip(array_dims, metadata["shape"]):
            start, stop = crop.get(name, (0, length))
            shape.append(stop - start)
            offset.append(start)
        return shape, tuple(offset)

    shape, _ = cropped(image_metadata)
    layout = resolve_layout(args, dims, shape)
    flag = find_flag(source, image, image_metadata, names, root_metadata)
    rewritten = {image} | ({flag} if flag else set())

    zarr.create_group(store=str(out), zarr_format=3, attributes=root_metadata.get("attributes", {}))
    for name in names:
        metadata = read_metadata(source / name)
        if metadata.get("node_type") != "array":
            shutil.copytree(source / name, out / name)
            continue
        new_shape, offset = cropped(metadata)
        if name in rewritten:
            itemsize = np.dtype(metadata["data_type"]).itemsize
            write_array_metadata(out / name, with_layout(metadata, new_shape, layout, itemsize))
            log(f"rewriting {name} {new_shape} as {layout.outer()}")
            job = Job(str(out / name), "copy", str(source / name), offset, keep_bits=layout.keep_bits)
            fill(job, new_shape, layout.outer(), itemsize, args)
        elif new_shape == metadata["shape"]:
            shutil.copytree(source / name, out / name)
        else:
            # Cropped but not rewritten: the same encoding, over the window. Chunks longer than the
            # new shape are left as they were; zarr pads the edge, as it does for any edge chunk.
            changed = json.loads(json.dumps(metadata))
            changed["shape"] = new_shape
            write_array_metadata(out / name, changed)
            unit = tuple(changed["chunk_grid"]["configuration"]["chunk_shape"])
            log(f"cropping {name} to {new_shape}")
            job = Job(str(out / name), "copy", str(source / name), offset)
            itemsize = _open(str(source / name), "r").dtype.itemsize
            fill(job, new_shape, unit, itemsize, args)

    return {
        "source": {"kind": "rewrite", "path": str(source), "image": image, "crop": {k: list(v) for k, v in crop.items()}},
        "image": image,
        "flag": flag,
        "dims": dims,
        "shape": shape,
        "itemsize": np.dtype(image_metadata["data_type"]).itemsize,
        "layout": layout,
        "consolidate": not args.no_consolidate,
    }


def parse_range(text: str | None, length: int) -> tuple[int, int] | None:
    if not text:
        return None
    start, sep, stop = text.partition(":")
    try:
        bounds = (int(start), int(stop))
    except ValueError:
        raise SystemExit(f"--flagged-range: {text!r} is not START:STOP") from None
    if not sep or not 0 <= bounds[0] < bounds[1] <= length:
        raise SystemExit(f"--flagged-range: {text!r} is not a range of the {length} channels")
    return bounds


def line_sources(args: argparse.Namespace, shape: list[int]) -> int:
    if args.line_sources is not None:
        return args.line_sources
    return round(math.prod(shape) / PIXELS_PER_LINE_SOURCE)


def synthesize(args: argparse.Namespace, out: Path) -> dict[str, Any]:
    template = Path(args.template).resolve()
    root_metadata = read_metadata(template)
    names = member_names(template)
    image = default_image(root_metadata)
    image_metadata = read_metadata(template / image)
    dims = image_metadata["dimension_names"]
    if tuple(dims) != AXES:
        raise SystemExit(f"the template's {image} is stored as {dims}; synthesis expects {list(AXES)}")

    requested = parse_axis_map(args.shape, "--shape")
    if "time" in requested and requested["time"] != 1:
        raise SystemExit("--shape: a synthetic cube has one time; CARTA displays no other kind")
    shape = [1, requested.get("frequency", 1), requested.get("polarization", 1), requested.get("l", 1), requested.get("m", 1)]
    if shape[2] > len(STOKES):
        raise SystemExit(f"--shape: at most {len(STOKES)} polarizations, {', '.join(STOKES)}")
    layout = resolve_layout(args, dims, shape)

    def values(name: str) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.asarray(zarr.open_array(str(template / name), mode="r")[...])

    _, length_f, length_p, length_l, length_m = shape
    frequencies = args.frequency_start + args.channel_width * np.arange(length_f)
    l_values, m_values = values("l"), values("m")
    l_axis = (np.arange(length_l) - length_l // 2) * (l_values[1] - l_values[0])
    m_axis = (np.arange(length_m) - length_m // 2) * (m_values[1] - m_values[0])

    rest = args.rest_frequency
    reference = root_metadata["attributes"]["coordinate_system_info"]["reference_direction"]["data"]

    direction = SkyDirection((float(reference[0]), float(reference[1])), (length_l, length_m),
                             (float(l_values[1] - l_values[0]), float(m_values[1] - m_values[0])))

    beam = values("BEAM_FIT_PARAMS_SKY") if "BEAM_FIT_PARAMS_SKY" in names else None
    coordinates: dict[str, np.ndarray] = {
        "time": values("time"),
        "frequency": frequencies,
        "velocity": SPEED_OF_LIGHT * (1.0 - frequencies / rest),
        "polarization": np.asarray(STOKES[:length_p], dtype="<U1"),
        "l": l_axis,
        "m": m_axis,
        "beam_params_label": values("beam_params_label") if "beam_params_label" in names else None,
    }
    if beam is not None:
        coordinates["BEAM_FIT_PARAMS_SKY"] = np.broadcast_to(beam[:, :1, :1, :], (1, length_f, length_p, beam.shape[3])).copy()
    unknown = [name for name in names if name != image and name not in coordinates and name not in SKY_DIRECTION]
    if unknown:
        raise SystemExit(f"the template has arrays this generator does not know how to stretch: {', '.join(unknown)}")

    synthetic = Synthetic(
        shape=tuple(shape),
        frequencies=tuple(float(value) for value in frequencies),
        seed=args.seed,
        noise=args.noise,
        line_sources=line_sources(args, shape),
        point_sources=args.point_sources,
        extended_sources=args.extended_sources,
        footprint_fill=None if args.footprint_fill >= 1 else args.footprint_fill,
        flagged_channels=args.flagged_channels,
        flag=args.flag,
        flagged_range=parse_range(args.flagged_range, shape[1]),
    )

    attributes = json.loads(json.dumps(root_metadata.get("attributes", {})))
    flag_name = "FLAG_SKY" if args.flag else None
    if flag_name:
        attributes.setdefault("data_groups", {}).setdefault("base", {})["flag"] = flag_name
    zarr.create_group(store=str(out), zarr_format=3, attributes=attributes)

    for name in names:
        if name == image:
            continue
        metadata = read_metadata(template / name)
        array_dims = metadata.get("dimension_names") or []
        if name in SKY_DIRECTION:
            if array_dims != ["l", "m"]:
                raise SystemExit(f"the template's {name} is stored as {array_dims}; synthesis expects ['l', 'm']")
            data_shape = [length_l, length_m]
        else:
            data_shape = list(coordinates[name].shape)
        metadata["shape"] = data_shape
        if name == "frequency":
            attributes = metadata.setdefault("attributes", {})
            if "rest_frequency" in attributes:
                attributes["rest_frequency"]["data"] = rest
            if "reference_frequency" in attributes:
                attributes["reference_frequency"]["data"] = float(frequencies[0])
        # Coordinates over l and m are as large as a plane, so they are chunked like the image's
        # planes; the rest are small enough to be one chunk, which is what XRADIO writes.
        chunk = [layout.chunk[AXES.index(d)] if d in ("l", "m") else n for d, n in zip(array_dims, data_shape)]
        metadata["chunk_grid"] = {"name": "regular", "configuration": {"chunk_shape": chunk or data_shape}}
        write_array_metadata(out / name, metadata)
        if name in SKY_DIRECTION:
            log(f"synthesizing {name}")
            fill(Job(str(out / name), name, direction=direction), data_shape, tuple(chunk),
                 np.dtype(metadata["data_type"]).itemsize, args)
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            zarr.open_array(str(out / name), mode="r+")[...] = coordinates[name]

    sky = with_layout(image_metadata, shape, layout, 4)
    sky["attributes"]["object_name"] = "carta-zarr-bench synthetic"
    sky["attributes"]["telescope"] = ASKAP_TELESCOPE
    if flag_name:
        sky["attributes"]["flag"] = flag_name
    write_array_metadata(out / image, sky)
    log(f"synthesizing {image} {shape} as {layout.outer()}")
    fill(Job(str(out / image), "pixels", synthetic=synthetic, keep_bits=layout.keep_bits), shape, layout.outer(), 4, args)

    if flag_name:
        flag_metadata = with_layout(
            {
                "zarr_format": 3,
                "node_type": "array",
                "data_type": "bool",
                "fill_value": False,
                "chunk_key_encoding": {"name": "default", "configuration": {"separator": "/"}},
                "dimension_names": list(AXES),
                "attributes": {"type": "flag"},
            },
            shape,
            layout,
            1,
        )
        write_array_metadata(out / flag_name, flag_metadata)
        log(f"synthesizing {flag_name}")
        fill(Job(str(out / flag_name), "flags", synthetic=synthetic), shape, layout.outer(), 1, args)

    return {
        "source": {
            "kind": "synthetic",
            "template": str(template),
            "seed": args.seed,
            "noise": args.noise,
            "frequency_start": args.frequency_start,
            "channel_width": args.channel_width,
            "rest_frequency": args.rest_frequency,
            "line_sources": synthetic.line_sources,
            "point_sources": args.point_sources,
            "extended_sources": args.extended_sources,
            "footprint_fill": synthetic.footprint_fill,
            "flagged_channels": args.flagged_channels,
            "flagged_range": args.flagged_range,
        },
        "image": image,
        "flag": flag_name,
        "dims": list(dims),
        "shape": shape,
        "itemsize": 4,
        "layout": layout,
        "consolidate": not args.no_consolidate,
    }


# -- Storage --------------------------------------------------------------------------------------


def filesystem_of(path: Path) -> str:
    """The type of the mount `path` is on, from /proc/mounts; "unknown" where there is none."""
    mounts = Path("/proc/mounts")
    if not mounts.is_file():
        return "unknown"
    resolved = str(path.resolve())
    best, kind = "", "unknown"
    for line in mounts.read_text().splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        mountpoint = fields[1].replace("\\040", " ")
        inside = resolved == mountpoint or resolved.startswith(mountpoint.rstrip("/") + "/")
        if inside and len(mountpoint) > len(best):
            best, kind = mountpoint, fields[2]
    return kind


def run(command: list[str]) -> str:
    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError:
        raise SystemExit(f"{command[0]} is not installed here, so {' '.join(command)} cannot run") from None
    if completed.returncode != 0:
        raise SystemExit(f"{' '.join(command)} failed ({completed.returncode}): {completed.stderr.strip()}")
    return completed.stdout


def stripe_commands(stripe: Stripe, directory: Path) -> tuple[list[str], list[str]]:
    """How to set this striping on `directory`, and how to read back what it got."""
    if stripe.filesystem == "lustre":
        size = stripe.size_with_unit("K", "M")
        return (["lfs", "setstripe", "-c", str(stripe.count), "-S", size, str(directory)],
                ["lfs", "getstripe", "-d", str(directory)])
    if shutil.which("beegfs-ctl"):  # BeeGFS 7
        size = stripe.size_with_unit("k", "m")
        return (["beegfs-ctl", "--setpattern", f"--numtargets={stripe.count}", f"--chunksize={size}", str(directory)],
                ["beegfs-ctl", "--getentryinfo", str(directory)])
    # BeeGFS 8 replaced beegfs-ctl with the beegfs tool.
    size = stripe.size_with_unit("KiB", "MiB")
    return (["beegfs", "entry", "set", f"--num-targets={stripe.count}", f"--chunk-size={size}", str(directory)],
            ["beegfs", "entry", "info", str(directory)])


def current_striping(filesystem: str, directory: Path) -> str | None:
    """What the directory's striping is, said by the filesystem's own tool, when there is one."""
    if filesystem == "lustre" and shutil.which("lfs"):
        command = ["lfs", "getstripe", "-d", str(directory)]
    elif filesystem == "beegfs" and shutil.which("beegfs-ctl"):
        command = ["beegfs-ctl", "--getentryinfo", str(directory)]
    elif filesystem == "beegfs" and shutil.which("beegfs"):
        command = ["beegfs", "entry", "info", str(directory)]
    else:
        return None
    completed = subprocess.run(command, capture_output=True, text=True)
    return completed.stdout.strip() if completed.returncode == 0 else None


def tree_size(path: Path) -> tuple[int, int]:
    total = files = 0
    for directory, _, entries in os.walk(path):
        for entry in entries:
            total += os.lstat(os.path.join(directory, entry)).st_size
            files += 1
    return total, files


# -- The source's own layout ----------------------------------------------------------------------


def codec_spelling(codecs: list[dict[str, Any]], where: str) -> str:
    """The --codec spelling of a chain of zarr v3 codecs this generator could have written: bytes,
    then at most one compressor. Anything else is refused rather than approximated, since a layout
    labelled the source's own that is not would make every comparison with it wrong."""
    names = [codec.get("name") for codec in codecs]
    if not names or names[0] != "bytes" or len(names) > 2:
        raise SystemExit(f"--layout-from-source: {where} has codecs {names}, which this generator cannot write")
    if len(names) == 1:
        return "none"
    compressor = codecs[1]
    config = compressor.get("configuration", {})
    if compressor["name"] == "zstd" and not config.get("checksum", False):
        return f"zstd:{config.get('level', 3)}"
    if compressor["name"] == "gzip":
        return f"gzip:{config.get('level', 6)}"
    if compressor["name"] == "blosc" and config.get("shuffle") in ("noshuffle", "shuffle", "bitshuffle"):
        return f"blosc:{config['cname']}:{config['clevel']}:{config['shuffle']}"
    raise SystemExit(f"--layout-from-source: {where} is compressed as {compressor}, which this generator cannot write")


def take_layout_from_source(args: argparse.Namespace) -> None:
    """Fill in --chunk, --shard, --codec and --no-consolidate as the source has them, so that the
    rewrite is the source's own layout over the crop -- what a data producer already writes, to
    compare every other layout with -- and has the identity it would have had spelt out by hand."""
    source = Path(args.source).resolve()
    root = read_metadata(source)
    image = args.image or default_image(root)
    metadata = read_metadata(source / image)
    dims = metadata.get("dimension_names") or []
    outer = metadata["chunk_grid"]["configuration"]["chunk_shape"]
    codecs = metadata["codecs"]
    shard = None
    chunk = outer
    if codecs and codecs[0].get("name") == "sharding_indexed":
        config = codecs[0]["configuration"]
        shard, chunk, codecs = outer, config["chunk_shape"], config["codecs"]
    args.chunk = ",".join(f"{name}={length}" for name, length in zip(dims, chunk))
    args.shard = ",".join(f"{name}={length}" for name, length in zip(dims, shard)) if shard else None
    args.codec = codec_spelling(codecs, f"{source.name}/{image}")
    args.no_consolidate = root.get("consolidated_metadata") is None
    log(f"the source's own layout: --chunk {args.chunk}"
        + (f" --shard {args.shard}" if args.shard else "")
        + f" --codec {args.codec}" + (" --no-consolidate" if args.no_consolidate else ""))


# -- Identity and reuse ---------------------------------------------------------------------------


def identity(args: argparse.Namespace) -> dict[str, Any]:
    """Everything that decides the bytes of the dataset, and nothing that does not (workers, block
    size, where it is written)."""
    if args.source:
        source = Path(args.source).resolve()
        image_metadata = read_metadata(source / (args.image or default_image(read_metadata(source))))
        origin = {
            "kind": "rewrite",
            "path": str(source),
            "image": args.image,
            "image_metadata": hashlib.sha256(json.dumps(image_metadata, sort_keys=True).encode()).hexdigest(),
            "content": source_content(source),
            "crop": args.crop,
        }
    else:
        template = Path(args.template).resolve()
        origin = {
            "kind": "synthetic",
            "template": hashlib.sha256((template / "zarr.json").read_bytes()).hexdigest(),
            # The coordinates and their metadata come from the template's arrays, not only its root
            # document: rechunked along frequency, it writes another dataset under the same root.
            "template_content": source_content(template),
            "shape": args.shape,
            "seed": args.seed,
            "noise": args.noise,
            "frequency_start": args.frequency_start,
            "channel_width": args.channel_width,
            "rest_frequency": args.rest_frequency,
            "line_sources": args.line_sources,
            "point_sources": args.point_sources,
            "extended_sources": args.extended_sources,
            "footprint_fill": args.footprint_fill,
            "flagged_channels": args.flagged_channels,
            "flagged_range": args.flagged_range,
            "flag": args.flag,
        }
    return {
        "format": FORMAT_VERSION,
        "origin": origin,
        "chunk": args.chunk,
        "shard": args.shard,
        "codec": Codec.parse(args.codec).spelling(),
        "keep_bits": args.keep_bits,
        "consolidated": not args.no_consolidate,
        "stripe": Stripe.parse(args.stripe).spelling() if args.stripe else None,
    }


def read_manifest(directory: Path) -> dict[str, Any] | None:
    try:
        return json.loads((directory / MANIFEST_NAME).read_text())
    except (OSError, ValueError):
        return None


def git_commit() -> str | None:
    try:
        completed = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], capture_output=True, text=True)
    except OSError:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


# -- Main -----------------------------------------------------------------------------------------


def parse_arguments(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    origin = parser.add_mutually_exclusive_group(required=True)
    origin.add_argument("--source", help="an XRADIO image dataset to rewrite")
    origin.add_argument("--synthetic", action="store_true", help="synthesize the pixels instead")

    where = parser.add_mutually_exclusive_group(required=True)
    where.add_argument("--output", help="the dataset to write")
    where.add_argument("--output-root", help="a directory to write it in, named after what decides its bytes")

    layout = parser.add_argument_group("layout")
    layout.add_argument("--chunk", help="chunk shape per axis, e.g. l=512,m=512,frequency=16; others are 1")
    layout.add_argument("--shard", help="shard shape per axis, a multiple of the chunk; others are the chunk")
    layout.add_argument("--codec", help="none, zstd[:level], gzip[:level], blosc[:cname[:level[:shuffle]]] (zstd:3)")
    layout.add_argument("--keep-bits", type=int, help="round pixels to this many float32 mantissa bits")
    layout.add_argument("--no-consolidate", action="store_true", help="leave out the root's consolidated metadata")
    layout.add_argument("--stripe", help="lustre:count=N,size=S or beegfs:count=N,size=S, set before writing")
    layout.add_argument("--layout-from-source", action="store_true",
                        help="the source's own chunk, shard, codec and consolidation, in place of those four")

    rewrite = parser.add_argument_group("rewriting a source")
    rewrite.add_argument("--image", help="the image to rewrite; the dataset's base sky image by default")
    rewrite.add_argument("--crop", help="a window per axis, e.g. frequency=0:1000,l=1024:3072")

    synthetic = parser.add_argument_group("synthesizing")
    synthetic.add_argument("--template", default=str(DEFAULT_TEMPLATE), help="an XRADIO dataset to take metadata from")
    synthetic.add_argument("--shape", help="frequency=F,polarization=P,l=L,m=M")
    synthetic.add_argument("--seed", type=int, default=1)
    synthetic.add_argument("--noise", type=float, default=ASKAP_NOISE,
                           help="noise rms at the footprint's centre, in the image's units")
    synthetic.add_argument("--frequency-start", type=float, default=ASKAP_FREQUENCY_START,
                           help="the first channel's frequency, in Hz (ASKAP's 1295.5 MHz)")
    synthetic.add_argument("--channel-width", type=float, default=ASKAP_CHANNEL_WIDTH,
                           help="in Hz (ASKAP's 18.5 kHz)")
    synthetic.add_argument("--rest-frequency", type=float, default=HI_REST_FREQUENCY, help="in Hz (HI's)")
    synthetic.add_argument("--line-sources", type=int,
                           help=f"HI-like line sources; one per {PIXELS_PER_LINE_SOURCE:.0e} pixels by default")
    synthetic.add_argument("--point-sources", type=int, default=0, help="continuum point sources")
    synthetic.add_argument("--extended-sources", type=int, default=0, help="extended continuum sources")
    synthetic.add_argument("--footprint-fill", type=float, default=0.76,
                           help="the fraction of a plane inside the footprint, NaN outside it; 1 for none")
    synthetic.add_argument("--flagged-channels", type=float, default=0.02,
                           help="the fraction of channels in flagged runs, flagged whole or in part")
    synthetic.add_argument("--flagged-range", metavar="START:STOP",
                           help="channels to flag whole besides those runs, e.g. 128:192 to cover whole chunks")
    synthetic.add_argument("--flag", action="store_true",
                           help="flag what is not data in a flag variable instead of writing NaN there")

    execution = parser.add_argument_group("execution")
    execution.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    execution.add_argument("--block-mib", type=int, default=64, help="how much each worker writes at a time")
    execution.add_argument("--force", action="store_true",
                           help="replace whatever dataset is at the output, this one included, rather than reuse it")

    args = parser.parse_args(argv)
    if args.synthetic and not args.shape:
        parser.error("--synthetic needs --shape")
    if not 0 < args.footprint_fill:
        parser.error("--footprint-fill must be positive")
    if not 0 <= args.flagged_channels <= 0.5:
        parser.error("--flagged-channels must be between 0 and 0.5")
    if args.source and args.shape:
        parser.error("--shape is for --synthetic; a rewrite takes its shape from the source, and --crop")
    if args.layout_from_source:
        if not args.source:
            parser.error("--layout-from-source needs --source: a synthetic cube has no layout of its own")
        given = [flag for flag, value in (("--chunk", args.chunk), ("--shard", args.shard), ("--codec", args.codec),
                                          ("--no-consolidate", args.no_consolidate), ("--keep-bits", args.keep_bits))
                 if value]
        if given:
            parser.error(f"--layout-from-source takes the layout from the source; leave out {', '.join(given)}")
        take_layout_from_source(args)
    elif not args.chunk:
        parser.error("--chunk is required, unless --layout-from-source takes it from the source")
    if args.codec is None:
        args.codec = "zstd:3"
    return args


def check_stripe(argv: list[str]) -> int:
    """Whether this striping can be set under the output root, tried on an empty directory that is
    removed again: so that sweep.py --dry-run finds a BeeGFS that reserves striping for root before
    the sweep has spent hours on the layouts that do not need it."""
    parser = argparse.ArgumentParser(prog="generate.py --check-stripe")
    parser.add_argument("--check-stripe", required=True, metavar="STRIPE")
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args(argv)
    stripe = Stripe.parse(args.check_stripe)
    root = Path(args.output_root).resolve()
    filesystem = filesystem_of(root if root.exists() else root.parent)
    if filesystem not in ("unknown", stripe.filesystem):
        raise SystemExit(f"--stripe is for {stripe.filesystem}, but {root} is on {filesystem}")
    root.mkdir(parents=True, exist_ok=True)
    probe = root / f".stripe-check-{os.getpid()}"
    probe.mkdir()
    try:
        set_command, _ = stripe_commands(stripe, probe)
        run(set_command)
    finally:
        shutil.rmtree(probe, ignore_errors=True)
    log(f"{stripe.spelling()} can be set under {root}")
    return 0


def overlaps(a: Path, b: Path) -> bool:
    """Whether either resolved path is the other or lies beneath it.

    Compared by file identity as well as by spelling, so that a hard-linked or case-folded alias of
    a directory is the directory, as resolve() alone does not settle on a case-insensitive disk.
    """

    def within(inner: Path, outer: Path) -> bool:
        for candidate in (inner, *inner.parents):
            if candidate == outer:
                return True
            if candidate.exists() and outer.exists() and os.path.samefile(candidate, outer):
                return True
        return False

    return within(a, b) or within(b, a)


def inputs_read(args: argparse.Namespace) -> list[tuple[str, Path]]:
    """The datasets a run reads, each with the flag that named it: the source it rewrites, or the
    template a synthetic cube takes its metadata from."""
    if args.source:
        return [("--source", Path(args.source).resolve())]
    return [("--template", Path(args.template).resolve())]


@dataclasses.dataclass(frozen=True)
class Output:
    """The dataset a run writes, claimed before anything at its path is deleted or made.

    Claiming is the only way this script deletes or creates an output, and it refuses an output over
    any dataset the run reads before it does either: --force deletes what is at the output, so an
    output that is the source, an alias of it or a directory around it would take the source with
    it, and one inside it would be written into the dataset being read. That used to be two
    statements of main() that had to stay in that order.
    """

    path: Path
    # Already this dataset, finished: nothing is to be written.
    reused: bool

    @classmethod
    def claim(cls, path: Path, reads: list[tuple[str, Path]], digest: str, force: bool) -> Output:
        """Refuse an output over `reads`; with `force`, replace whatever is there; without it, reuse a
        finished dataset whose manifest says `digest` and refuse any other; and leave an empty
        directory to write in.

        `force` replaces even this dataset: it is how to write one again when what the identity can
        see of its origin is not all there is to it."""
        for flag, read in reads:
            # What a link in the dataset points at is the dataset too, however far from its own
            # directory: an output that holds it would delete those bytes or write among them.
            for held in (read, *(entry.resolve() for entry in reached(read) if entry.is_symlink())):
                if overlaps(path, held):
                    where = "" if held == read else f" through a link to {held}"
                    raise SystemExit(f"{path} overlaps {flag} {read}{where}; write the dataset somewhere else")
        if path.exists():
            manifest = read_manifest(path)
            if not force and manifest and manifest.get("identity_hash") == digest and manifest.get("complete"):
                return cls(path, reused=True)
            if not force:
                raise SystemExit(f"{path} exists and is not this dataset; --force replaces it")
            shutil.rmtree(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.mkdir()
        return cls(path, reused=False)

    @contextlib.contextmanager
    def writing(self) -> Iterator[Path]:
        """The output to write in. Anything short of finishing removes it again, an interrupt or a
        refusal part-way included: an unfinished dataset is of no use to a measurement, and left
        behind it would make the same command refuse to run a second time."""
        try:
            yield self.path
        except BaseException:
            shutil.rmtree(self.path, ignore_errors=True)
            raise


def main(argv: list[str] | None = None) -> int:
    arguments = sys.argv[1:] if argv is None else argv
    if any(argument.split("=", 1)[0] == "--check-stripe" for argument in arguments):
        return check_stripe(arguments)
    args = parse_arguments(argv)
    stripe = Stripe.parse(args.stripe) if args.stripe else None
    key = identity(args)
    digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]
    if args.output:
        out = Path(args.output).resolve()
    else:
        stem = "synthetic" if args.synthetic else Path(args.source).resolve().name.removesuffix(".zarr")
        out = Path(args.output_root).resolve() / f"{stem}-{digest}"
    output = Output.claim(out, inputs_read(args), digest, args.force)
    if output.reused:
        log(f"{out} is already this dataset; reusing it")
        print(out)
        return 0

    with output.writing():
        filesystem = filesystem_of(out)
        if stripe:
            if filesystem not in ("unknown", stripe.filesystem):
                raise SystemExit(f"--stripe is for {stripe.filesystem}, but {out} is on {filesystem}")
            set_command, _ = stripe_commands(stripe, out)
            log(f"striping: {' '.join(set_command)}")
            run(set_command)

        started = time.monotonic()
        made = rewrite_source(args, out) if args.source else synthesize(args, out)
        if made["consolidate"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                zarr.consolidate_metadata(str(out))
        seconds = time.monotonic() - started

    layout: Layout = made["layout"]
    image_bytes, image_files = tree_size(out / made["image"])
    total_bytes, total_files = tree_size(out)
    itemsize = made["itemsize"]
    uncompressed = math.prod(made["shape"]) * itemsize
    manifest = {
        "manifest_version": 1,
        "complete": True,
        "identity_hash": digest,
        "identity": key,
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "generator": {
            "format": FORMAT_VERSION,
            "carta_zarr_commit": git_commit(),
            "zarr_python": zarr.__version__,
            "numpy": np.__version__,
            "python": platform.python_version(),
            "host": platform.node(),
        },
        "source": made["source"],
        "image": {
            "variable": made["image"],
            "flag": made["flag"],
            "dimension_names": made["dims"],
            "shape": made["shape"],
            "itemsize": itemsize,
        },
        "layout": {
            "chunk_shape": list(layout.chunk),
            "shard_shape": list(layout.shard) if layout.shard else None,
            "codec": layout.codec.spelling(),
            "keep_bits": layout.keep_bits,
            "consolidated": made["consolidate"],
            "chunk_bytes": math.prod(layout.chunk) * itemsize,
            "file_bytes": math.prod(layout.outer()) * itemsize,
        },
        "storage": {
            "filesystem": filesystem,
            "stripe_requested": stripe.spelling() if stripe else None,
            "stripe_effective": current_striping(filesystem, out),
        },
        "result": {
            "uncompressed_bytes": uncompressed,
            "image_stored_bytes": image_bytes,
            "image_files": image_files,
            "compression_ratio": uncompressed / image_bytes if image_bytes else None,
            "stored_bytes": total_bytes,
            "files": total_files,
            "write_seconds": round(seconds, 3),
        },
    }
    # Written last and renamed into place, so that a dataset with a manifest is a finished one.
    partial = out / (MANIFEST_NAME + ".partial")
    partial.write_text(json.dumps(manifest, indent=2) + "\n")
    partial.replace(out / MANIFEST_NAME)
    ratio = manifest["result"]["compression_ratio"]
    log(f"wrote {out} in {seconds:.1f} s: {image_files} files, compression ratio "
        f"{ratio:.2f}" if ratio else f"wrote {out} in {seconds:.1f} s")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
