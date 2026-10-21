# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# generate.py writes a synthetic cube as FITS, and verify.py checks a Zarr made from it. The test set's
# Zarr comes from convert.py, which needs XRADIO, but generate.py writes the same cube as a dataset
# too, so verifying the two against each other checks both: the order of the pixels, and what
# verify.py holds a layout to. A Zarr in other chunks than asked for, with a pixel changed, or without
# the flag it was meant to carry -- or with one carta-zarr would not apply, because the image does not
# declare it, or one over a pixel that is data -- must then fail verification.

if(NOT DEFINED UV OR NOT DEFINED SOURCE_DIR OR NOT DEFINED OUTPUT_DIR OR NOT DEFINED BENCH)
    message(FATAL_ERROR "UV, SOURCE_DIR, OUTPUT_DIR and BENCH must be set")
endif()
# verify.py reads every Zarr through carta-zarr-bench.
set(ENV{CARTA_ZARR_BENCH} "${BENCH}")

file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

function(run what)
    execute_process(COMMAND "${UV}" run --quiet --script ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(result)
        message(FATAL_ERROR "${what} failed: ${result}\n${out}\n${err}")
    endif()
endfunction()

# Refused, and for the reason given: a check whose name starts so, reported DIFFERENT. So that a
# fixture that fails some other way does not pass unnoticed.
function(refused_for what reason)
    execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/verify.py" ${ARGN}
        RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
    string(REGEX MATCH "(^|\n)${reason}[^\n]*: DIFFERENT" found "${out}")
    if(NOT result OR NOT found)
        message(FATAL_ERROR "verify.py did not refuse ${what} for ${reason}: ${result}\n${out}\n${err}")
    endif()
endfunction()

# A Zarr changed by `code`, Python run with `root` set to it, and consolidated again as a converter
# would, so that what is changed is what carta-zarr reads.
function(change dataset code)
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                python -c "import json, numpy as np, zarr\nroot = '${dataset}'\n${code}\nzarr.consolidate_metadata(root)"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "changing ${dataset} failed: ${result}")
    endif()
endfunction()

set(generate "${SOURCE_DIR}/tools/zarr-bench/generate.py")
set(verify "${SOURCE_DIR}/tools/testset/verify.py")
set(fits "${OUTPUT_DIR}/cube.fits")
set(synthetic --synthetic --shape frequency=40,polarization=1,l=90,m=70 --flagged-channels 0.1 --workers 1)
set(cube "${OUTPUT_DIR}/cube.zarr")
run("generate.py --fits" "${generate}" ${synthetic} --fits "${fits}")
run("generate.py" "${generate}" ${synthetic} --chunk l=32,m=32,frequency=8 --output "${cube}")
run("verify.py" "${verify}" "${fits}" "${cube}" --flag no --chunks 32,32,8)

refused_for("a Zarr in chunks it was not asked for" "chunk" "${fits}" "${cube}" --flag no --chunks 64,32,8)
refused_for("a Zarr unsharded where a shard was asked for" "chunk" "${fits}" "${cube}" --flag no --chunks 32,32,8
            --shards 64,64,8)
refused_for("a Zarr without the flag it was meant to carry" "a flag carta-zarr applies" "${fits}" "${cube}" --flag yes
            --chunks 32,32,8)
set(sharded "${OUTPUT_DIR}/sharded.zarr")
run("generate.py sharded" "${generate}" ${synthetic} --chunk l=32,m=32,frequency=8 --shard l=64,m=64,frequency=8
    --output "${sharded}")
run("verify.py sharded" "${verify}" "${fits}" "${sharded}" --flag no --chunks 32,32,8 --shards 64,64,8)

# One pixel changed, away from any plane or spectrum one might sample.
file(COPY "${cube}/" DESTINATION "${OUTPUT_DIR}/changed.zarr")
change("${OUTPUT_DIR}/changed.zarr" "zarr.open_array(root + '/SKY', mode='r+')[0, 7, 0, 61, 13] = 1.0")
refused_for("a Zarr whose pixels differ from the FITS cube" "every pixel" "${fits}" "${OUTPUT_DIR}/changed.zarr" --flag no
            --chunks 32,32,8)

# Copies of the cube with a flag beside it, true where the pixel is NaN, as --compute_mask computes
# one. Typed and declared it is the flag carta-zarr applies, and verifies; untyped and undeclared,
# carta-zarr ignores it, and so must verify.py.
function(with_flag name typed declared)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${cube}/" DESTINATION "${copy}")
    change("${copy}" "sky = zarr.open_array(root + '/SKY', mode='r')
flag = zarr.create_array(root + '/FLAG_SKY', shape=sky.shape, chunks=sky.chunks, dtype=bool,
                         dimension_names=list(sky.metadata.dimension_names), fill_value=False)
flag[...] = np.isnan(sky[...])
if ${typed}:
    flag.attrs['type'] = 'flag'
if ${declared}:
    meta = json.load(open(root + '/zarr.json'))
    meta['attributes']['data_groups']['base']['flag'] = 'FLAG_SKY'
    json.dump(meta, open(root + '/zarr.json', 'w'))")
endfunction()
with_flag(flagged True True)
run("verify.py with a flag" "${verify}" "${fits}" "${OUTPUT_DIR}/flagged.zarr" --flag yes --chunks 32,32,8)
refused_for("a Zarr with a flag it was not meant to carry" "no flag carta-zarr applies" "${fits}"
            "${OUTPUT_DIR}/flagged.zarr" --flag no --chunks 32,32,8)
with_flag(undeclared False False)
refused_for("a flag carta-zarr does not apply" "a flag carta-zarr applies" "${fits}" "${OUTPUT_DIR}/undeclared.zarr"
            --flag yes --chunks 32,32,8)
run("verify.py beside a flag carta-zarr ignores" "${verify}" "${fits}" "${OUTPUT_DIR}/undeclared.zarr" --flag no
    --chunks 32,32,8)
# A flag over a pixel that is data, which the pixels as stored do not show.
file(COPY "${OUTPUT_DIR}/flagged.zarr/" DESTINATION "${OUTPUT_DIR}/overflagged.zarr")
change("${OUTPUT_DIR}/overflagged.zarr" "sky = zarr.open_array(root + '/SKY', mode='r')
assert np.isfinite(sky[0, 7, 0, 45, 35])
zarr.open_array(root + '/FLAG_SKY', mode='r+')[0, 7, 0, 45, 35] = True")
refused_for("a flag over a pixel that is data" "flagged where NaN" "${fits}" "${OUTPUT_DIR}/overflagged.zarr" --flag yes
            --chunks 32,32,8)

# stats.py counts a chunk whole, the padding at the image's edge included, since that is what the codec
# compressed; and of a shard only the chunks its index lists. Uncompressed, a chunk is then as large on
# disk as it decodes to: shape 5 in chunks of 4 is two chunks of 16 bytes, 1.0, where counting only
# the pixels inside the image gave 0.625.
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
            python -c "import importlib.util, pathlib, numpy as np, zarr
spec = importlib.util.spec_from_file_location('stats', '${SOURCE_DIR}/tools/testset/stats.py')
stats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stats)
root = pathlib.Path('${OUTPUT_DIR}/stats')
plain = zarr.create_array(str(root / 'plain'), shape=(5,), chunks=(4,), dtype='float32', compressors=None, fill_value=np.nan)
plain[:] = 1.0
assert stats.compression(root / 'plain') == 1.0, stats.compression(root / 'plain')
sharded = zarr.create_array(str(root / 'sharded'), shape=(20,), chunks=(4,), shards=(8,), dtype='float32', compressors=None,
                            fill_value=np.nan)
sharded[:5] = 1.0
# One shard on disk holding two of its chunks, 32 bytes, behind an index of 2 x 16 bytes and a checksum.
assert stats.compression(root / 'sharded') == 32 / 68, stats.compression(root / 'sharded')
sparse = zarr.create_array(str(root / 'sparse'), shape=(20,), chunks=(4,), shards=(8,), dtype='float32', compressors=None,
                           fill_value=np.nan)
sparse[:4] = 1.0
# The same shard with its second chunk left out: 16 bytes, not the shard's 32.
assert stats.compression(root / 'sparse') == 16 / 52, stats.compression(root / 'sparse')
# Nothing on disk at all, as a cube flagged whole leaves: no ratio to give, rather than a division by zero.
zarr.create_array(str(root / 'empty'), shape=(20,), chunks=(4,), dtype='float32', compressors=None, fill_value=np.nan)
assert stats.compression(root / 'empty') is None
# The same chunks named by the other encodings, c.0 and 0, are found as well.
for name, encoding in (('dotted', {'name': 'default', 'separator': '.'}), ('v2', {'name': 'v2', 'separator': '.'})):
    other = zarr.create_array(str(root / name), shape=(5,), chunks=(4,), dtype='float32', compressors=None, fill_value=np.nan,
                              chunk_key_encoding=encoding)
    other[:] = 1.0
    assert stats.compression(root / name) == 1.0, (name, stats.compression(root / name))
# And nothing that is not a chunk, beside the chunks or among a shard's.
for name in ('plain', 'sharded', 'v2'):
    (root / name / 'README.txt').write_text('x')
    (root / name / 'zarr.json~').write_text('x')
(root / 'sharded' / 'c' / 'notes').write_text('x')
assert stats.compression(root / 'plain') == 1.0
assert stats.compression(root / 'sharded') == 32 / 68
assert stats.compression(root / 'v2') == 1.0"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(result)
    message(FATAL_ERROR "stats.py compression failed: ${result}\n${out}\n${err}")
endif()
