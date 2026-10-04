"""What a source dataset holds, as generate.py and sweep.py both have to say it.

One function in a module of its own because the two must agree: the generator decides by it whether
a rewrite it already made can be reused, and the sweep whether runs it already measured were of the
same data. Imported from the scripts' own directory, which Python puts first on the path of a script
it runs, and which needs no package.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


def source_content(source: Path) -> str:
    """What every file of a source dataset is, as far as the filesystem says without reading it: its
    path, size and modification time. A rewrite reads pixels, flags and coordinates as well as
    metadata, and a source edited in place keeps its path and, as often as not, its metadata, so
    without this the rewrite of what it used to hold is reused.

    Stated rather than read because a source can be terabytes: walking it costs a stat per file, and
    hashing it would cost reading all of it. A file copied without its times looks changed, which
    costs a rewrite and never a wrong measurement."""
    digest = hashlib.sha256()
    for path in sorted(entry for entry in source.rglob("*") if entry.is_file()):
        stat = path.stat()
        digest.update(f"{path.relative_to(source).as_posix()}\0{stat.st_size}\0{stat.st_mtime_ns}\n".encode())
    return digest.hexdigest()
