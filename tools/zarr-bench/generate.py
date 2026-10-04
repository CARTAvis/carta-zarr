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

  --synthetic      Noise, point sources with a spectral index, and NaN outside a circle -- a
                   primary-beam-corrected cube, roughly. The coordinates are built from a template
                   XRADIO wrote (the conformance fixture by default) and stretched to --shape.
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

from fingerprint import source_content

MANIFEST_NAME = "bench-manifest.json"
# Bumped whenever the same arguments would produce different bytes, so that a dataset written by an
# older generator is not mistaken for one this one would write.
FORMAT_VERSION = 2

AXES = ("time", "frequency", "polarization", "l", "m")
STOKES = ("I", "Q", "U", "V")
SPEED_OF_LIGHT = 299_792_458.0

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TEMPLATE = REPO_ROOT / "tests" / "data" / "images" / "zarr" / "xradio" / "conformance"

# The square of noise one seed draws at a time. Fixed, and independent of any layout, because it is
# what makes a synthetic pixel depend on where it is rather than on which block wrote it.
NOISE_TILE = 128


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

    The noise is drawn a NOISE_TILE square at a time from a generator seeded by that tile's position,
    and the sources are a fixed table each block consults, so the value of a pixel does not depend on
    the block that wrote it. That is what makes two layouts of one synthetic cube comparable.
    """

    shape: tuple[int, int, int, int, int]  # time, frequency, polarization, l, m
    frequencies: tuple[float, ...]
    seed: int
    noise: float
    sources: int
    nan_radius: float | None
    flag: bool

    def table(self) -> dict[str, np.ndarray]:
        rng = np.random.default_rng([self.seed, 0x5EED])
        _, _, _, length_l, length_m = self.shape
        return {
            "l": rng.uniform(0, length_l, self.sources),
            "m": rng.uniform(0, length_m, self.sources),
            "flux": 10 ** rng.uniform(-3, 0, self.sources),
            "sigma": rng.uniform(1.5, 4.0, self.sources),
            "index": rng.normal(-0.7, 0.3, self.sources),
        }

    def outside(self, l_range: range, m_range: range) -> np.ndarray | None:
        if self.nan_radius is None:
            return None
        _, _, _, length_l, length_m = self.shape
        radius = self.nan_radius * min(length_l, length_m) / 2
        l = np.asarray(l_range, dtype=np.float64)[:, None] - (length_l - 1) / 2
        m = np.asarray(m_range, dtype=np.float64)[None, :] - (length_m - 1) / 2
        return l * l + m * m > radius * radius

    def pixels(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
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
                            noise *= np.float32(self.noise)
                            la, lb = max(l0, tl * tile), min(l1, (tl + 1) * tile)
                            ma, mb = max(m0, tm * tile), min(m1, (tm + 1) * tile)
                            plane[la - l0 : lb - l0, ma - m0 : mb - m0] = noise[
                                la - tl * tile : lb - tl * tile, ma - tm * tile : mb - tm * tile
                            ]

        # Sources are in Stokes I alone, and every source reaches the same pixels whatever the block,
        # because its extent is clipped to five sigma before it is clipped to the block.
        if p0 == 0 and self.sources:
            table = self.table()
            spectrum_base = np.asarray(self.frequencies[f0:f1], dtype=np.float64) / self.frequencies[0]
            for index in range(self.sources):
                sigma = table["sigma"][index]
                centre_l, centre_m = table["l"][index], table["m"][index]
                la = max(l0, math.floor(centre_l - 5 * sigma))
                lb = min(l1, math.ceil(centre_l + 5 * sigma) + 1)
                ma = max(m0, math.floor(centre_m - 5 * sigma))
                mb = min(m1, math.ceil(centre_m + 5 * sigma) + 1)
                if la >= lb or ma >= mb:
                    continue
                profile_l = np.exp(-0.5 * ((np.arange(la, lb) - centre_l) / sigma) ** 2)
                profile_m = np.exp(-0.5 * ((np.arange(ma, mb) - centre_m) / sigma) ** 2)
                spectrum = table["flux"][index] * spectrum_base ** table["index"][index]
                blob = spectrum[:, None, None] * np.outer(profile_l, profile_m)[None, :, :]
                block[:, :, 0, la - l0 : lb - l0, ma - m0 : mb - m0] += blob.astype(np.float32)[None]

        if not self.flag:
            outside = self.outside(range(l0, l1), range(m0, m1))
            if outside is not None:
                block[..., outside] = np.nan
        return block

    def flags(self, start: tuple[int, ...], stop: tuple[int, ...]) -> np.ndarray:
        # True is a flagged pixel: XRADIO's flag says which pixels to drop. See src/pixel_mask.h.
        block = np.zeros([b - a for a, b in zip(start, stop)], dtype=bool)
        outside = self.outside(range(start[3], stop[3]), range(start[4], stop[4]))
        if outside is not None:
            block[..., outside] = True
        return block


# -- Writing --------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Job:
    """One array to fill, block by block: everything a worker process needs to do its share."""

    target: str  # path of the array being written
    kind: str  # "copy", "pixels" or "flags"
    source: str | None = None  # array read from, for a copy
    offset: tuple[int, ...] = ()  # where the target's origin sits in the source, for a crop
    synthetic: Synthetic | None = None
    keep_bits: int | None = None


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


def find_flag(root: Path, image: str, image_metadata: dict[str, Any], names: list[str]) -> str | None:
    """The flag carta-zarr would mask this image with: the one it declares, or else the only
    flag-typed boolean array with its dimensions and shape. The same rule as src/schema/xradio/flag.cc."""
    declared = image_metadata.get("attributes", {}).get("flag")
    if isinstance(declared, str) and declared:
        return declared
    matches = []
    for name in names:
        if name == image:
            continue
        metadata = read_metadata(root / name)
        if (
            metadata.get("node_type") == "array"
            and metadata.get("attributes", {}).get("type") == "flag"
            and metadata.get("data_type") == "bool"
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
    flag = find_flag(source, image, image_metadata, names)
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
    frequency = values("frequency")
    increment = frequency[1] - frequency[0] if frequency.size > 1 else 1.0e6
    frequencies = frequency[0] + increment * np.arange(length_f)
    l_values, m_values = values("l"), values("m")
    l_axis = (np.arange(length_l) - length_l // 2) * (l_values[1] - l_values[0])
    m_axis = (np.arange(length_m) - length_m // 2) * (m_values[1] - m_values[0])

    frequency_attributes = read_metadata(template / "frequency").get("attributes", {})
    rest = frequency_attributes.get("rest_frequency", {}).get("data", frequencies[0])
    reference = root_metadata["attributes"]["coordinate_system_info"]["reference_direction"]["data"]

    def sky_direction() -> tuple[np.ndarray, np.ndarray]:
        # The SIN projection inverted, as XRADIO's right_ascension and declination are.
        ra0, dec0 = reference
        l, m = l_axis[:, None], m_axis[None, :]
        n = np.sqrt(np.maximum(0.0, 1.0 - l * l - m * m))
        dec = np.arcsin(m * math.cos(dec0) + n * math.sin(dec0))
        ra = ra0 + np.arctan2(l, n * math.cos(dec0) - m * math.sin(dec0))
        return ra, dec

    beam = values("BEAM_FIT_PARAMS_SKY") if "BEAM_FIT_PARAMS_SKY" in names else None
    ra, dec = sky_direction()
    coordinates: dict[str, np.ndarray] = {
        "time": values("time"),
        "frequency": frequencies,
        "velocity": SPEED_OF_LIGHT * (1.0 - frequencies / rest),
        "polarization": np.asarray(STOKES[:length_p], dtype="<U1"),
        "l": l_axis,
        "m": m_axis,
        "right_ascension": ra,
        "declination": dec,
        "beam_params_label": values("beam_params_label") if "beam_params_label" in names else None,
    }
    if beam is not None:
        coordinates["BEAM_FIT_PARAMS_SKY"] = np.broadcast_to(beam[:, :1, :1, :], (1, length_f, length_p, beam.shape[3])).copy()
    unknown = [name for name in names if name != image and name not in coordinates]
    if unknown:
        raise SystemExit(f"the template has arrays this generator does not know how to stretch: {', '.join(unknown)}")

    synthetic = Synthetic(
        shape=tuple(shape),
        frequencies=tuple(float(value) for value in frequencies),
        seed=args.seed,
        noise=args.noise,
        sources=args.sources,
        nan_radius=None if args.nan_radius <= 0 else args.nan_radius,
        flag=args.flag,
    )

    attributes = json.loads(json.dumps(root_metadata.get("attributes", {})))
    flag_name = "FLAG_SKY" if args.flag else None
    if flag_name:
        attributes.setdefault("data_groups", {}).setdefault("base", {})["flag"] = flag_name
    zarr.create_group(store=str(out), zarr_format=3, attributes=attributes)

    for name in names:
        if name == image:
            continue
        data = coordinates[name]
        metadata = read_metadata(template / name)
        metadata["shape"] = list(data.shape)
        array_dims = metadata.get("dimension_names") or []
        # Coordinates over l and m are as large as a plane, so they are chunked like the image's
        # planes; the rest are small enough to be one chunk, which is what XRADIO writes.
        chunk = [layout.chunk[AXES.index(d)] if d in ("l", "m") else n for d, n in zip(array_dims, data.shape)]
        metadata["chunk_grid"] = {"name": "regular", "configuration": {"chunk_shape": chunk or list(data.shape)}}
        write_array_metadata(out / name, metadata)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            zarr.open_array(str(out / name), mode="r+")[...] = data

    sky = with_layout(image_metadata, shape, layout, 4)
    sky["attributes"]["object_name"] = "carta-zarr-bench synthetic"
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
            "sources": args.sources,
            "nan_radius": synthetic.nan_radius,
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
            "shape": args.shape,
            "seed": args.seed,
            "noise": args.noise,
            "sources": args.sources,
            "nan_radius": args.nan_radius,
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
    synthetic.add_argument("--noise", type=float, default=1.0e-3, help="noise sigma, in the image's units")
    synthetic.add_argument("--sources", type=int, default=200)
    synthetic.add_argument("--nan-radius", type=float, default=1.0,
                           help="NaN outside this fraction of the half-width; 0 for none")
    synthetic.add_argument("--flag", action="store_true",
                           help="flag the outside of that circle in a flag variable instead of writing NaN there")

    execution = parser.add_argument_group("execution")
    execution.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    execution.add_argument("--block-mib", type=int, default=64, help="how much each worker writes at a time")
    execution.add_argument("--force", action="store_true",
                           help="replace whatever dataset is at the output, this one included, rather than reuse it")

    args = parser.parse_args(argv)
    if args.synthetic and not args.shape:
        parser.error("--synthetic needs --shape")
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
            if overlaps(path, read):
                raise SystemExit(f"{path} overlaps {flag} {read}; write the dataset somewhere else")
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
