#!/usr/bin/env python3

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Empty the caches between this host and a Lustre file system's disks, for every file under DIR,
without root.

  drop-lustre-cache.py DIR

For carta-zarr-bench --drop-cache-cmd, or sweep.py's measure.drop_cache_cmd with DIR its paths.work,
where a sweep keeps one dataset at a time. Two caches stand between a read and the disks, and each
needs its own call:

  - this client's page cache, emptied by posix_fadvise(DONTNEED) on each file, which is what
    carta-zarr-bench --cold fadvise does;
  - each OSS's read cache, which no client call reaches but `lfs ladvise -a dontneed`, a request to
    the servers that Lustre 2.9 and later accept from the file's owner.

The second matters: on a Lustre 2.15 system, reading a whole cube took 1.2 s warm, 2.0 s with the
client's cache emptied, and 5.5 s with both -- the OSS cache alone made it look nearly three times
faster than a first read. Other users' reads of the same servers are untouched, and so are the
servers' disks' own caches.

Written for the system python3 rather than as a uv script like the tools beside it: it runs before
every trial, where uv's start-up would be paid each time, and uses nothing but the standard library
of Python 3.6 on.
"""

import os
import shutil
import subprocess
import sys

# Files per lfs ladvise, to stay well inside the argument-length limit.
BATCH = 500


def main(argv):
    if len(argv) != 2 or argv[1] in ("-h", "--help"):
        print(__doc__.strip(), file=sys.stderr)
        return 2
    root = argv[1]
    if not os.path.isdir(root):
        print("drop-lustre-cache: {} is not a directory".format(root), file=sys.stderr)
        return 1
    lfs = shutil.which("lfs")
    if lfs is None:
        print("drop-lustre-cache: lfs is not on the path; is this a Lustre client?", file=sys.stderr)
        return 1

    paths = []
    for directory, _, names in os.walk(root):
        for name in names:
            path = os.path.join(directory, name)
            try:
                fd = os.open(path, os.O_RDONLY)
            except OSError:
                continue
            try:
                # Dirty pages are not dropped, and a dataset just written may still have some.
                os.fdatasync(fd)
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(fd)
            paths.append(path)

    for start in range(0, len(paths), BATCH):
        completed = subprocess.run([lfs, "ladvise", "-a", "dontneed"] + paths[start:start + BATCH])
        if completed.returncode != 0:
            print("drop-lustre-cache: lfs ladvise failed ({})".format(completed.returncode), file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
