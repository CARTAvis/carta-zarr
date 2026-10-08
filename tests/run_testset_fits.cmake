# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# tools/testset/zarr-to-fits.py writes a synthetic cube as FITS, and verify.py checks a FITS cube against
# a Zarr. The test set's Zarr comes from a converter this repository does not have, but generate.py's
# own dataset holds the same pixels and coordinates, so writing it as FITS and verifying the two
# against each other checks both scripts: the axes, the order of the pixels, and the world coordinates.
# A Zarr whose pixels differ in any one place must then fail verification, as must one without the
# flag it was meant to carry, or with one carta-zarr would not apply.

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
run("verify.py" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag no)

function(refused what)
    execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/verify.py" ${ARGN}
        RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
    if(NOT result)
        message(FATAL_ERROR "verify.py passed ${what}")
    endif()
endfunction()

refused("a Zarr without the flag it was meant to carry" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag yes)

# Copies of the cube with a flag beside it, true where the pixel is NaN, as the converter's
# --compute_mask writes one: declared and typed it is the flag carta-zarr applies, and verifies; with
# no type, as xradio 1.2.2's FITS reader leaves it, carta-zarr ignores it, and verify.py must too.
function(with_flag name typed)
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
    meta = json.load(open(root + '/zarr.json'))
    meta['attributes']['data_groups']['base']['flag'] = 'FLAG_SKY'
    json.dump(meta, open(root + '/zarr.json', 'w'))
"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "adding a flag to ${name} failed: ${result}")
    endif()
endfunction()

with_flag(flagged True)
run("verify.py with a flag" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flagged.zarr" --flag yes)
refused("a Zarr whose flag it was not meant to carry" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flagged.zarr" --flag no)
with_flag(untyped False)
refused("a flag carta-zarr would not apply" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/untyped.zarr" --flag yes)

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
refused("a Zarr whose pixels differ from the FITS cube" "${OUTPUT_DIR}/cube.fits" "${changed}" --flag no)
