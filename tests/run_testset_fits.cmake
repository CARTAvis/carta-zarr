# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# tools/testset/zarr-to-fits.py writes a synthetic cube as FITS, and verify.py checks a FITS cube against
# a Zarr. The test set's Zarr comes from a converter this repository does not have, but generate.py's
# own dataset holds the same pixels and coordinates, so writing it as FITS and verifying the two
# against each other checks both scripts: the axes, the order of the pixels, and the world coordinates.
# A Zarr whose pixels differ in any one place must then fail verification, as must one without the
# flag it was meant to carry, or with one carta-zarr would not apply; and one whose axes, chunks, data
# type, frequencies or projection are not what the comparison takes them to be.

if(NOT DEFINED UV OR NOT DEFINED SOURCE_DIR OR NOT DEFINED OUTPUT_DIR)
    message(FATAL_ERROR "UV, SOURCE_DIR and OUTPUT_DIR must be set")
endif()

file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

function(run what)
    execute_process(COMMAND "${UV}" run --quiet --script ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(result)
        message(FATAL_ERROR "${what} failed: ${result}\n${out}\n${err}")
    endif()
endfunction()

set(cube "${OUTPUT_DIR}/cube.zarr")
run("generate.py" "${SOURCE_DIR}/tools/zarr-bench/generate.py" --synthetic --shape frequency=40,polarization=1,l=90,m=70
    --chunk l=32,m=32,frequency=8 --flagged-channels 0.1 --workers 1 --output "${cube}")
run("zarr-to-fits.py" "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${cube}" "${OUTPUT_DIR}/cube.fits" --block-mib 1)
run("verify.py" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag no --chunks 32,32,8)

function(refused what)
    execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/verify.py" ${ARGN}
        RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
    if(NOT result)
        message(FATAL_ERROR "verify.py passed ${what}")
    endif()
endfunction()

refused("a Zarr without the flag it was meant to carry" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag yes --chunks 32,32,8)
refused("a Zarr in chunks it was not asked for" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag no --chunks 64,32,8)
refused("a Zarr unsharded where a shard was asked for" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag no --chunks 32,32,8
        --shards 64,64,8)

# A copy of the cube changed by `code`, Python run with `root` set to the copy.
function(changed name code)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${cube}/" DESTINATION "${copy}")
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                python -c "import json, numpy as np, zarr\nroot = '${copy}'\n${code}"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "changing ${name} failed: ${result}")
    endif()
endfunction()

# carta-zarr takes the axes by name, so l and m swapped in the names is a transposed image.
changed(swapped "meta = json.load(open(root + '/SKY/zarr.json'))
meta['dimension_names'] = ['time', 'frequency', 'polarization', 'm', 'l']
json.dump(meta, open(root + '/SKY/zarr.json', 'w'))")
refused("a Zarr whose l and m are named the other way round" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/swapped.zarr"
        --flag no --chunks 32,32,8)
changed(middle "f = zarr.open_array(root + '/frequency', mode='r+')
f[20] = f[20] + 1e8")
refused("a Zarr with one channel's frequency wrong" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/middle.zarr"
        --flag no --chunks 32,32,8)
changed(tan "meta = json.load(open(root + '/zarr.json'))
meta['attributes']['coordinate_system_info']['projection'] = 'TAN'
json.dump(meta, open(root + '/zarr.json', 'w'))")
refused("a Zarr declaring a projection the sky check does not invert" "${OUTPUT_DIR}/cube.fits"
        "${OUTPUT_DIR}/tan.zarr" --flag no --chunks 32,32,8)
changed(pole "meta = json.load(open(root + '/zarr.json'))
meta['attributes']['coordinate_system_info']['native_pole_direction']['data'][0] = 0.0
json.dump(meta, open(root + '/zarr.json', 'w'))")
refused("a Zarr about another native pole" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/pole.zarr" --flag no --chunks 32,32,8)
changed(rotated "meta = json.load(open(root + '/zarr.json'))
meta['attributes']['coordinate_system_info']['pixel_coordinate_transformation_matrix'] = [[0.0, -1.0], [1.0, 0.0]]
json.dump(meta, open(root + '/zarr.json', 'w'))")
changed(equinox "meta = json.load(open(root + '/zarr.json'))
meta['attributes']['coordinate_system_info']['reference_direction']['attrs']['equinox'] = 'j1950.0'
json.dump(meta, open(root + '/zarr.json', 'w'))")
refused("a Zarr in another equinox" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/equinox.zarr" --flag no --chunks 32,32,8)
# The same equinox in every form carta-zarr reads.
foreach(form 2000.0 "'2000.0'" "'J2000'" "'b2000'" "'J2e3'" "'+2000'")
    string(MAKE_C_IDENTIFIER "equinox_${form}" name)
    changed(${name} "meta = json.load(open(root + '/zarr.json'))
meta['attributes']['coordinate_system_info']['reference_direction']['attrs']['equinox'] = ${form}
json.dump(meta, open(root + '/zarr.json', 'w'))")
    run("verify.py with the equinox written ${form}" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
        "${OUTPUT_DIR}/${name}.zarr" --flag no --chunks 32,32,8)
endforeach()
refused("a Zarr rotated on the sky" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/rotated.zarr" --flag no --chunks 32,32,8)
changed(nonfinite "f = zarr.open_array(root + '/frequency', mode='r+')
f[0] = np.nan")
refused("a Zarr with a channel's frequency not a number" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/nonfinite.zarr"
        --flag no --chunks 32,32,8)
changed(integer "meta = json.load(open(root + '/SKY/zarr.json'))
meta['data_type'] = 'uint32'
meta['fill_value'] = 0
json.dump(meta, open(root + '/SKY/zarr.json', 'w'))")
refused("a Zarr of integers with the pixels' bits" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/integer.zarr" --flag no
        --chunks 32,32,8)

# Copies of the cube with a flag beside it, true where the pixel is NaN, as the converter's
# --compute_mask writes one: declared and typed it is the flag carta-zarr applies, and verifies; with
# no type, as xradio 1.2.2's FITS reader leaves it, carta-zarr ignores it, and verify.py must too;
# untyped and declared, carta-zarr refuses to open the image, and verify.py refuses it even unflagged.
function(with_flag name typed declared)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${cube}/" DESTINATION "${copy}")
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                python -c "
import json, numpy as np, zarr
root = '${copy}'
sky = zarr.open_array(root + '/SKY', mode='r')
flag = zarr.create_array(root + '/FLAG_SKY', shape=sky.shape, chunks=sky.chunks, dtype=bool,
                         dimension_names=list(sky.metadata.dimension_names), fill_value=False)
flag[...] = np.isnan(sky[...])
if ${typed}:
    flag.attrs['type'] = 'flag'
if ${declared}:
    meta = json.load(open(root + '/zarr.json'))
    meta['attributes']['data_groups']['base']['flag'] = 'FLAG_SKY'
    json.dump(meta, open(root + '/zarr.json', 'w'))
"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "adding a flag to ${name} failed: ${result}")
    endif()
endfunction()

with_flag(flagged True True)
run("verify.py with a flag" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flagged.zarr" --flag yes
    --chunks 32,32,8)
refused("a Zarr whose flag it was not meant to carry" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flagged.zarr" --flag no
        --chunks 32,32,8)
with_flag(untyped False False)
refused("a flag carta-zarr would not apply" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/untyped.zarr" --flag yes --chunks 32,32,8)
run("verify.py beside a flag carta-zarr ignores" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/untyped.zarr" --flag no --chunks 32,32,8)
with_flag(declared_untyped False True)
refused("a declared flag carta-zarr refuses the image over" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/declared_untyped.zarr"
        --flag no --chunks 32,32,8)

# Declarations as carta-zarr reads them: names compared as the store files them, groups naming two
# flags refused, and a flag another image's group declares never taken for this one's.
function(declared name code)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${OUTPUT_DIR}/flagged.zarr/" DESTINATION "${copy}")
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                python -c "import json, shutil\nroot = '${copy}'\nmeta = json.load(open(root + '/zarr.json'))\ngroups = meta['attributes']['data_groups']\n${code}\njson.dump(meta, open(root + '/zarr.json', 'w'))"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "declaring ${name} failed: ${result}")
    endif()
endfunction()
declared(aliased "groups['other'] = {'sky': 'SKY', 'flag': './FLAG_SKY'}")
run("verify.py with one flag spelt two ways" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/aliased.zarr" --flag yes --chunks 32,32,8)
declared(conflicting "shutil.copytree(root + '/FLAG_SKY', root + '/OTHER_FLAG')
groups['other'] = {'sky': 'SKY', 'flag': 'OTHER_FLAG'}")
refused("data groups naming two usable flags for the image" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/conflicting.zarr"
        --flag yes --chunks 32,32,8)
declared(absolute "groups['base']['flag'] = '/FLAG_SKY'")
refused("a flag declared by an absolute path" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/absolute.zarr" --flag yes
        --chunks 32,32,8)
declared(climbing "shutil.copytree(root + '/FLAG_SKY', root + '/unused')
groups['base']['flag'] = 'unused/../FLAG_SKY'")
refused("a flag declared through '..'" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/climbing.zarr" --flag yes --chunks 32,32,8)
declared(owned "del groups['base']['flag']
groups['other'] = {'sky': 'OTHER', 'flag': './FLAG_SKY'}")
refused("a flag another image's group declares, taken for this one" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/owned.zarr"
        --flag yes --chunks 32,32,8)

# The same cube with one pixel changed, away from any plane or spectrum one might sample.
set(changed "${OUTPUT_DIR}/changed.zarr")
file(COPY "${cube}/" DESTINATION "${changed}")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
            python -c "import zarr; a = zarr.open_array('${changed}/SKY', mode='r+'); a[0, 7, 0, 61, 13] = 1.0"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "changing a pixel failed: ${result}")
endif()
refused("a Zarr whose pixels differ from the FITS cube" "${OUTPUT_DIR}/cube.fits" "${changed}" --flag no --chunks 32,32,8)

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
assert stats.compression(root / 'sparse') == 16 / 52, stats.compression(root / 'sparse')"
    RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(result)
    message(FATAL_ERROR "stats.py compression failed: ${result}\n${out}\n${err}")
endif()
