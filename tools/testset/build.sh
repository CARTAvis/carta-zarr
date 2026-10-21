#!/usr/bin/env bash
# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# Build the standard read-path test set in OUTPUT: each cube as a FITS file, written by generate.py,
# and each of its Zarr layouts made from that FITS by XRADIO's own open_image and write_image, through
# convert.py. XRADIO_CONVERTER names another converter that takes the same arguments.
#
#   CARTA_ZARR_BENCH=/path/to/carta-zarr-bench tools/testset/build.sh OUTPUT
#
# What exists is skipped, so a run that stopped part-way is finished by running it again; a FITS
# file or a Zarr appears under its name only once it is complete, and a Zarr only once verify.py has
# found it to hold the FITS file's pixels. See README.md for what the cubes are and why.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
generate="$here/../zarr-bench/generate.py"
output=${1:?usage: CARTA_ZARR_BENCH=... $0 OUTPUT}
converter=${XRADIO_CONVERTER:-$here/convert.py}
[ -f "$converter" ] || { echo "no converter at $converter" >&2; exit 1; }
# verify.py reads each Zarr through carta-zarr itself, as carta-backend would.
bench=${CARTA_ZARR_BENCH:?set CARTA_ZARR_BENCH to carta-zarr-bench, from a build with -DCARTA_ZARR_BUILD_BENCH=ON}
[ -x "$bench" ] || { echo "no carta-zarr-bench at $bench" >&2; exit 1; }
export PATH="$HOME/.local/bin:$PATH"
mkdir -p "$output"

# name | shape (l, m, frequency) | generator options beyond the defaults
cubes=(
    "pancake|7763,4742,256|--flagged-range 128:192"
    "cigar|512,512,30000|"
)

# cube | chunk (l, m, frequency) | shard (l, m, frequency), or - | flag: yes or no |
# compressor: blosc-zstd or zstd
layouts=(
    "pancake|256,256,16|-|no|blosc-zstd"
    "pancake|512,512,4|-|no|blosc-zstd"
    "pancake|256,256,16|-|yes|blosc-zstd"
    "cigar|128,128,64|-|no|blosc-zstd"
    "cigar|64,64,256|-|no|blosc-zstd"
)

say() { echo "== $(date '+%F %T') $*" >&2; }

for entry in "${cubes[@]}"; do
    IFS='|' read -r cube shape options <<< "$entry"
    IFS=',' read -r l m f <<< "$shape"
    fits="$output/$cube.fits"
    [ -e "$fits" ] && continue
    say "$cube: synthesizing $l x $m x $f as $fits"
    # shellcheck disable=SC2086  # the options are words
    uv run --quiet "$generate" --synthetic --shape "frequency=$f,polarization=1,l=$l,m=$m" $options \
        --fits "$fits" > /dev/null
done

for entry in "${layouts[@]}"; do
    IFS='|' read -r cube chunk shard flag compressor <<< "$entry"
    IFS=',' read -r cl cm cf <<< "$chunk"
    name="${cube}_c${cl}x${cm}x${cf}"
    arguments=(--chunks "0,$cf,0,$cl,$cm" --compressor "$compressor")
    expected=(--chunks "$chunk")
    if [ "$shard" != - ]; then
        IFS=',' read -r sl sm sf <<< "$shard"
        if (( sl % cl || sm % cm || sf % cf )); then
            echo "$name: shard $shard is not a whole number of chunks $chunk" >&2
            exit 1
        fi
        name+="_s${sl}x${sm}x${sf}"
        arguments+=(--shards "0,$sf,0,$sl,$sm")
        expected+=(--shards "$shard")
    fi
    if [ "$flag" = yes ]; then
        name+="_flag"
        arguments+=(--compute_mask)
    fi
    [ "$compressor" = zstd ] && name+="_zstd"
    zarr="$output/$name.zarr"
    [ -e "$zarr" ] && continue
    partial="$output/$name.partial.zarr"
    rm -rf "$partial"
    say "$name: converting $cube.fits"
    uv run --quiet --python-preference only-managed "$converter" "$output/$cube.fits" "${arguments[@]}" \
        --output "$partial" 2> >(grep -v -i warning >&2)
    say "$name: verifying"
    uv run --quiet "$here/verify.py" "$output/$cube.fits" "$partial" --flag "$flag" "${expected[@]}" --bench "$bench" \
        > "$output/$name.verify.txt" || {
        cat "$output/$name.verify.txt" >&2
        exit 1
    }
    mv "$partial" "$zarr"
done

say "done"
