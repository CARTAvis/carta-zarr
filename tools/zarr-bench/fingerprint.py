"""What a source dataset holds, as generate.py and sweep.py both have to say it.

One function in a module of its own because the two must agree: the generator decides by it whether
a rewrite it already made can be reused, and the sweep whether runs it already measured were of the
same data. Imported from the scripts' own directory, which Python puts first on the path of a script
it runs, and which needs no package.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterator
from pathlib import Path


def reached(root: Path) -> Iterator[Path]:
    """Every file and directory under `root`, spelt as reached through it, with directory links
    followed: a dataset's arrays are as often as not links to where their bytes are, and what is
    behind them is the dataset as much as what is not. Each directory is entered once, so a link back
    up the tree is listed and not followed round again."""
    entered = {os.path.realpath(root)}
    for directory, directories, files in os.walk(root, followlinks=True):
        here = Path(directory)
        kept = []
        for name in sorted(directories):
            yield here / name
            real = os.path.realpath(here / name)
            if real not in entered:
                entered.add(real)
                kept.append(name)
        directories[:] = kept
        for name in sorted(files):
            yield here / name


def source_content(source: Path) -> str:
    """What every file of a source dataset is, as far as the filesystem says without reading it: its
    path, size and modification time. A rewrite reads pixels, flags and coordinates as well as
    metadata, and a source edited in place keeps its path and, as often as not, its metadata, so
    without this the rewrite of what it used to hold is reused.

    Stated rather than read because a source can be terabytes: walking it costs a stat per file, and
    hashing it would cost reading all of it. A file copied without its times looks changed, which
    costs a rewrite and never a wrong measurement."""
    digest = hashlib.sha256()
    for path in sorted(entry for entry in reached(source) if entry.is_file()):
        stat = path.stat()
        digest.update(f"{path.relative_to(source).as_posix()}\0{stat.st_size}\0{stat.st_mtime_ns}\n".encode())
    return digest.hexdigest()
