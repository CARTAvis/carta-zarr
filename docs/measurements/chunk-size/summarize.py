#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas>=2.2"]
# ///
"""Reduces the bench results to medians.csv: one row per machine, cube, layout, user count, cache and mode.

    summarize.py MACHINE=RESULTS_DIR ... > medians.csv

Each RESULTS_DIR holds a CUBE.csv per cube, which carta-zarr-bench wrote with --csv and labelled
--label "LAYOUT|users=N|cache=SIZE", with "|warm" after the planes read without emptying any cache.
plane, spectrum and region take the median of their first touches -- operations that read no chunk
another operation of the trial read before or alongside them -- as sweep.py ranks them, and leave it
empty when fewer than three operations, one a trial, were first touches; so does plane_warm, the
Mac's check that purge emptied its caches. An animation reports the share of its frames that were
late, and the medians of its first, typical and slowest frames.
"""

import sys
from pathlib import Path

import pandas as pd

MIN_FIRST_TOUCHES = 3


def load(machine: str, results: Path) -> pd.DataFrame:
    frames = []
    for path in sorted(results.glob("*.csv")):
        frame = pd.read_csv(path, low_memory=False)
        frame["machine"], frame["cube"] = machine, path.stem
        frames.append(frame)
    return pd.concat(frames)


def axis(shape: pd.Series, name: str) -> pd.Series:
    return shape.fillna("").str.extract(rf"{name}=(\d+)")[0]


def summarize(rows: pd.DataFrame) -> pd.DataFrame:
    label = rows.label.str.split("|", expand=True)
    rows = rows.assign(
        users=label[1].str.removeprefix("users=").astype(int),
        cache=label[2].str.removeprefix("cache="),
        mode=rows["mode"].where(label[3].fillna("") != "warm", "plane_warm"),
        chunk=axis(rows.chunk_shape, "l") + "x" + axis(rows.chunk_shape, "m") + "x"
        + axis(rows.chunk_shape, "frequency"),
        shard=(axis(rows.shard_shape, "l") + "x" + axis(rows.shard_shape, "m") + "x"
               + axis(rows.shard_shape, "frequency")).fillna(""),
        chunk_mib=rows.chunk_decoded_bytes / 2**20,
    )
    keys = ["machine", "cube", "chunk", "shard", "chunk_mib", "users", "cache", "mode"]
    out = []
    for key, group in rows.groupby(keys, sort=False):
        row = dict(zip(keys, key))
        ok = group[group.status == "ok"]
        row["ops"], row["failed"] = len(group), len(group) - len(ok)
        if row["mode"] == "animation":
            row["late_pct"] = 100 * ok.late_frames.sum() / max(ok.frames_played.sum(), 1)
            row["frame_first_s"] = ok.frame_first_s.median()
            row["frame_median_s"] = ok.frame_median_s.median()
            row["frame_max_s"] = ok.frame_max_s.median()
            row["prefetches"] = ok.prefetches.median()
        else:
            touched = ok[~ok.shares_chunks.astype(bool)]
            row["first_touches"] = len(touched)
            if len(touched) >= MIN_FIRST_TOUCHES:
                row["median_s"], row["p90_s"] = touched.seconds.median(), touched.seconds.quantile(0.9)
        out.append(row)
    columns = keys + ["ops", "failed", "first_touches", "median_s", "p90_s", "late_pct", "frame_first_s",
                      "frame_median_s", "frame_max_s", "prefetches"]
    return pd.DataFrame(out).reindex(columns=columns)


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    frames = []
    for argument in sys.argv[1:]:
        machine, _, results = argument.partition("=")
        frames.append(load(machine, Path(results)))
    medians = summarize(pd.concat(frames))
    medians["first_touches"] = medians.first_touches.astype("Int64")
    medians.to_csv(sys.stdout, index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
