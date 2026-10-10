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

# A copy of the cube changed by `code`, Python run with `root` set to the copy, and consolidated again
# as a converter would, so that what is changed is what carta-zarr reads.
function(changed name code)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${cube}/" DESTINATION "${copy}")
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                --with astropy==7.1.0 python -c "import json, numpy as np, zarr\nroot = '${copy}'\n${code}\nzarr.consolidate_metadata(root)"
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
# The same on a square cube, where the swap changes no length and only the pixels can tell.
set(square "${OUTPUT_DIR}/square.zarr")
run("generate.py square" "${SOURCE_DIR}/tools/zarr-bench/generate.py" --synthetic --shape frequency=8,polarization=1,l=48,m=48
    --chunk l=16,m=16,frequency=8 --workers 1 --output "${square}")
run("zarr-to-fits.py square" "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${square}" "${OUTPUT_DIR}/square.fits" --block-mib 1)
file(COPY "${square}/" DESTINATION "${OUTPUT_DIR}/square_swapped.zarr")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 python -c "import json, zarr
root = '${OUTPUT_DIR}/square_swapped.zarr'
meta = json.load(open(root + '/SKY/zarr.json'))
meta['dimension_names'] = ['time', 'frequency', 'polarization', 'm', 'l']
json.dump(meta, open(root + '/SKY/zarr.json', 'w'))
zarr.consolidate_metadata(root)"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "swapping the square cube's names failed: ${result}")
endif()
run("verify.py square" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/square.fits" "${square}" --flag no --chunks 16,16,8)
refused("a square Zarr whose l and m are named the other way round" "${OUTPUT_DIR}/square.fits"
        "${OUTPUT_DIR}/square_swapped.zarr" --flag no --chunks 16,16,8)
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
changed(stokes "p = zarr.open_array(root + '/polarization', mode='r+')
p[0] = 'Q'")
refused("a Zarr calling Stokes I Q" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/stokes.zarr" --flag no --chunks 32,32,8)
changed(timeless "import shutil
shutil.rmtree(root + '/time')")
refused("a Zarr without a time coordinate" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/timeless.zarr" --flag no --chunks 32,32,8)
changed(later "t = zarr.open_array(root + '/time', mode='r+')
t[0] = t[0] + 1.0")
refused("a Zarr a day later than the FITS file" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/later.zarr" --flag no --chunks 32,32,8)
changed(rest "f = zarr.open_array(root + '/frequency', mode='r+')
f.attrs['rest_frequency'] = dict(f.attrs['rest_frequency'], data=1.4e9)")
refused("a Zarr with another rest frequency" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/rest.zarr" --flag no --chunks 32,32,8)
changed(observer "f = zarr.open_array(root + '/frequency', mode='r+')
reference = dict(f.attrs['reference_frequency'])
reference['attrs'] = dict(reference['attrs'], observer='bary')
f.attrs['reference_frequency'] = reference")
refused("a Zarr in another spectral frame" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/observer.zarr" --flag no --chunks 32,32,8)
changed(degenerate "l = zarr.open_array(root + '/l', mode='r+')
l[1] = l[0]")
refused("a Zarr whose first two l are the same" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/degenerate.zarr" --flag no
        --chunks 32,32,8)
# A node patched after the store was consolidated: carta-zarr reads the consolidated copy and refuses one
# that disagrees with the node's own zarr.json.
set(stale "${OUTPUT_DIR}/stale.zarr")
file(COPY "${cube}/" DESTINATION "${stale}")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed python -c "import json
path = '${stale}/SKY/zarr.json'
meta = json.load(open(path))
meta['attributes']['units'] = 'K'
json.dump(meta, open(path, 'w'))"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "patching stale failed: ${result}")
endif()
refused("a Zarr patched after it was consolidated" "${OUTPUT_DIR}/cube.fits" "${stale}" --flag no --chunks 32,32,8)
set(twice "${OUTPUT_DIR}/twice.zarr")
file(COPY "${cube}/" DESTINATION "${twice}")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed python -c "import json
path = '${twice}/zarr.json'
meta = json.load(open(path))
block = meta['consolidated_metadata']['metadata']
block['./SKY'] = block['SKY']
json.dump(meta, open(path, 'w'))"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "listing a node twice failed: ${result}")
endif()
refused("a Zarr whose consolidated metadata lists a node twice" "${OUTPUT_DIR}/cube.fits" "${twice}" --flag no
        --chunks 32,32,8)
changed(misnamed "meta = json.load(open(root + '/frequency/zarr.json'))
meta['dimension_names'] = ['m']
json.dump(meta, open(root + '/frequency/zarr.json', 'w'))")
refused("a Zarr whose frequency coordinate is named along m" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/misnamed.zarr"
        --flag no --chunks 32,32,8)
changed(gigahertz "f = zarr.open_array(root + '/frequency', mode='r+')
f.attrs['units'] = 'GHz'")
refused("a Zarr whose frequencies are said to be in GHz" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/gigahertz.zarr"
        --flag no --chunks 32,32,8)
changed(typed "a = zarr.open_array(root + '/SKY', mode='r+')
a.attrs['type'] = 'flag'")
refused("a Zarr whose SKY is typed flag" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/typed.zarr" --flag no --chunks 32,32,8)
changed(kelvin "a = zarr.open_array(root + '/SKY', mode='r+')
a.attrs['units'] = 'K'")
refused("a Zarr in another brightness unit" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/kelvin.zarr" --flag no --chunks 32,32,8)
changed(obsdate "a = zarr.open_array(root + '/SKY', mode='r+')
date = dict(a.attrs['obsdate'])
date['data'] = date['data'] + 1.0
a.attrs['obsdate'] = date")
refused("a Zarr observed a day later by its image's obsdate" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/obsdate.zarr" --flag no
        --chunks 32,32,8)
changed(beamless "a = zarr.open_array(root + '/SKY', mode='r+')
del a.attrs['beam_fit_params']")
refused("a Zarr without the FITS file's restoring beam" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/beamless.zarr" --flag no
        --chunks 32,32,8)
changed(widebeam "b = zarr.open_array(root + '/BEAM_FIT_PARAMS_SKY', mode='r+')
values = b[...]
values[0, 5, 0, 0] *= 2
b[...] = values")
refused("a Zarr whose beam is wider in one channel" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/widebeam.zarr" --flag no
        --chunks 32,32,8)
# The same pixels in shards whose chunks are stored transposed. The layout is what carta-zarr reads,
# in the image's axis order, so it verifies as that layout.
changed(transposed "import shutil
from zarr.codecs import TransposeCodec
old = zarr.open_array(root + '/SKY', mode='r')
values, attributes, names = old[...], old.attrs.asdict(), list(old.metadata.dimension_names)
shutil.rmtree(root + '/SKY')
sky = zarr.create_array(root + '/SKY', shape=values.shape, chunks=(1, 8, 1, 32, 16), shards=(1, 8, 1, 64, 64), dtype='float32',
                        filters=[TransposeCodec(order=(0, 1, 2, 4, 3))], dimension_names=names, fill_value=np.nan)
sky[...] = values
sky.attrs.update(attributes)")
run("verify.py with transposed chunks in shards" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/transposed.zarr" --flag no --chunks 32,16,8 --shards 64,64,8)
refused("a sharded Zarr verified as another layout" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/transposed.zarr" --flag no
        --chunks 16,32,8 --shards 64,64,8)
# The same pixels stored m before l: carta-zarr takes the axes by name, so this is the same image.
changed(reordered "import shutil
old = zarr.open_array(root + '/SKY', mode='r')
values, attributes = old[...], old.attrs.asdict()
shutil.rmtree(root + '/SKY')
sky = zarr.create_array(root + '/SKY', shape=values.shape[:3] + values.shape[:2:-1], chunks=(1, 8, 1, 32, 32), dtype='float32',
                        dimension_names=['time', 'frequency', 'polarization', 'm', 'l'], fill_value=np.nan)
sky[...] = values.transpose(0, 1, 2, 4, 3)
sky.attrs.update(attributes)")
run("verify.py with m stored before l" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/reordered.zarr" --flag no --chunks 32,32,8)
# SKY's dimension names kept in its attributes, where some XRADIO writers put them and carta-zarr looks.
changed(attributed "meta = json.load(open(root + '/SKY/zarr.json'))
meta['attributes']['dimension_names'] = meta.pop('dimension_names')
json.dump(meta, open(root + '/SKY/zarr.json', 'w'))")
run("verify.py with SKY's dimension names in its attributes" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/attributed.zarr" --flag no --chunks 32,32,8)
# A beam on as many planes as the image has, but not on its planes: half the channels, twice the
# polarizations.
changed(halfbeam "import shutil
old = zarr.open_array(root + '/BEAM_FIT_PARAMS_SKY', mode='r')
values, attributes, names = old[...], old.attrs.asdict(), list(old.metadata.dimension_names)
shutil.rmtree(root + '/BEAM_FIT_PARAMS_SKY')
shape = list(values.shape)
shape[names.index('frequency')] //= 2
shape[names.index('polarization')] *= 2
beam = zarr.create_array(root + '/BEAM_FIT_PARAMS_SKY', shape=shape, dtype=values.dtype, dimension_names=names)
beam[...] = values.reshape(shape)
beam.attrs.update(attributes)")
refused("a Zarr whose beams cover as many planes as it has, but other ones" "${OUTPUT_DIR}/cube.fits"
        "${OUTPUT_DIR}/halfbeam.zarr" --flag no --chunks 32,32,8)

# Pixels zarr-python decodes and carta-zarr cannot: numcodecs' zlib, which TensorStore does not register.
changed(zlibbed "import shutil
from numcodecs.zarr3 import Zlib
old = zarr.open_array(root + '/SKY', mode='r')
values, attributes, names = old[...], old.attrs.asdict(), list(old.metadata.dimension_names)
shutil.rmtree(root + '/SKY')
sky = zarr.create_array(root + '/SKY', shape=values.shape, chunks=old.chunks, dtype='float32', compressors=[Zlib(level=1)],
                        dimension_names=names, fill_value=np.nan)
sky[...] = values
sky.attrs.update(attributes)")
refused("a Zarr whose pixels carta-zarr cannot decode" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/zlibbed.zarr" --flag no
        --chunks 32,32,8)
# Frequencies written in GHz: carta-zarr gives the rest frequency in the axis's unit too, so one left in
# Hz beside them is a rest frequency a billion times too high, which RESTFRQ, in Hz, must expose.
changed(ghz "f = zarr.open_array(root + '/frequency', mode='r+')
f[...] = f[...] / 1e9
f.attrs['units'] = 'GHz'
reference = dict(f.attrs['reference_frequency'])
reference['data'] = reference['data'] / 1e9
f.attrs['reference_frequency'] = reference
rest = dict(f.attrs['rest_frequency'])
rest['data'] = rest['data'] / 1e9
f.attrs['rest_frequency'] = rest")
run("verify.py with the frequencies in GHz" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits"
    "${OUTPUT_DIR}/ghz.zarr" --flag no --chunks 32,32,8)
changed(ghz_rest_hz "f = zarr.open_array(root + '/frequency', mode='r+')
f[...] = f[...] / 1e9
f.attrs['units'] = 'GHz'
reference = dict(f.attrs['reference_frequency'])
reference['data'] = reference['data'] / 1e9
f.attrs['reference_frequency'] = reference")
refused("a Zarr in GHz with its rest frequency left in Hz" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/ghz_rest_hz.zarr" --flag no
        --chunks 32,32,8)
# The FITS side's linear transformation whole: CD alike to CDELT and PC, and sky axes mixed with the
# spectral one, which carta-zarr's direction cannot hold, refused.
function(fits_variant name code)
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with astropy==7.1.0 python -c "from astropy.io import fits
header = fits.getheader('${OUTPUT_DIR}/cube.fits')
data = fits.getdata('${OUTPUT_DIR}/cube.fits')
${code}
fits.writeto('${OUTPUT_DIR}/${name}.fits', data, header, overwrite=True)"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "writing ${name}.fits failed: ${result}")
    endif()
endfunction()
fits_variant(cd "for i in range(1, 5):
    header[f'CD{i}_{i}'] = header.pop(f'CDELT{i}')
for key in ('PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'):
    del header[key]")
fits_variant(coupled "header['PC1_4'] = 0.1")
run("verify.py against a FITS file in CD" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cd.fits" "${cube}" --flag no
    --chunks 32,32,8)
refused("a FITS file whose RA moves with frequency" "${OUTPUT_DIR}/coupled.fits" "${cube}" --flag no --chunks 32,32,8)
# A carta-zarr whose linear frequency axis were a channel off -- what carta-backend builds when it is
# given one -- with the table of frequencies right, which only a stand-in for the bench can make.
file(WRITE "${OUTPUT_DIR}/misfitted.py" "import json, subprocess, sys
result = subprocess.run(['${BENCH}'] + sys.argv[1:], capture_output=True)
out = result.stdout
if sys.argv[1] == 'probe' and '--describe' in sys.argv and not result.returncode:
    report = json.loads(out)
    report['image']['description']['spectral']['reference_pixel'] += 1
    out = json.dumps(report).encode()
sys.stdout.buffer.write(out)
sys.stderr.buffer.write(result.stderr)
sys.exit(result.returncode)
")
file(WRITE "${OUTPUT_DIR}/misfitted.sh" "#!/bin/sh\nexec \"${UV}\" run --quiet --no-project --python-preference only-managed python \"${OUTPUT_DIR}/misfitted.py\" \"$@\"\n")
file(CHMOD "${OUTPUT_DIR}/misfitted.sh" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
refused("a linear frequency axis a channel off" "${OUTPUT_DIR}/cube.fits" "${cube}" --flag no --chunks 32,32,8
        --bench "${OUTPUT_DIR}/misfitted.sh")
# obsdate as an ISO date, which carta-zarr reads, and as a numeric string, which it does not.
changed(isodate "a = zarr.open_array(root + '/SKY', mode='r+')
from astropy.time import Time
date = dict(a.attrs['obsdate'])
date['data'] = Time(date['data'], format='mjd', scale='utc').isot
a.attrs['obsdate'] = date")
run("verify.py with an ISO obsdate" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/isodate.zarr"
    --flag no --chunks 32,32,8)
changed(textdate "a = zarr.open_array(root + '/SKY', mode='r+')
date = dict(a.attrs['obsdate'])
date['data'] = str(date['data'])
a.attrs['obsdate'] = date")
refused("a Zarr whose obsdate is a number written as text" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/textdate.zarr" --flag no
        --chunks 32,32,8)
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
zarr.consolidate_metadata(root)
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
# The right flag in chunks other than the image's.
set(rechunked "${OUTPUT_DIR}/rechunked.zarr")
file(COPY "${OUTPUT_DIR}/flagged.zarr/" DESTINATION "${rechunked}")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
            python -c "import numpy as np, shutil, zarr
root = '${rechunked}'
old = zarr.open_array(root + '/FLAG_SKY', mode='r')
values, attributes, names = old[...], old.attrs.asdict(), list(old.metadata.dimension_names)
shutil.rmtree(root + '/FLAG_SKY')
flag = zarr.create_array(root + '/FLAG_SKY', shape=values.shape, chunks=(1, 8, 1, 16, 16), dtype=bool, dimension_names=names,
                         fill_value=False)
flag[...] = values
flag.attrs.update(attributes)
zarr.consolidate_metadata(root)"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "rechunking the flag failed: ${result}")
endif()
refused("a flag in chunks other than the image's" "${OUTPUT_DIR}/cube.fits" "${rechunked}" --flag yes --chunks 32,32,8)

# A flag beside SKY's shards, SKY's and the flag's each given as (transpose, inner l, inner m): a
# transpose "inner" reorders bytes within a chunk and changes no footprint, so a flag without one holds
# the same pixels in each chunk and verifies; one "outer", ahead of the shards, swaps the axes their
# chunk_shape is given in. zarr-python then reads both arrays as unsharded 64 x 64 chunks, whatever
# their inner chunks, so a flag cut 16 x 32 beside a SKY cut 32 x 16 must be told apart by the codecs.
function(sharded_flag name sky flag)
    set(copy "${OUTPUT_DIR}/${name}.zarr")
    file(COPY "${OUTPUT_DIR}/flagged.zarr/" DESTINATION "${copy}")
    execute_process(
        COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
                python -W ignore -c "import json, numpy as np, shutil, zarr
from zarr.codecs import TransposeCodec
root = '${copy}'
swap = (0, 1, 2, 4, 3)
for node, fill, (transpose, l, m) in (('SKY', np.nan, ${sky}), ('FLAG_SKY', False, ${flag})):
    old = zarr.open_array(root + '/' + node, mode='r')
    values, attributes, names = old[...], old.attrs.asdict(), list(old.metadata.dimension_names)
    shutil.rmtree(root + '/' + node)
    inner = (1, 8, 1, m, l) if transpose == 'outer' else (1, 8, 1, l, m)
    new = zarr.create_array(root + '/' + node, shape=values.shape, chunks=inner, shards=(1, 8, 1, 64, 64),
                            dtype=values.dtype, filters=[TransposeCodec(order=swap)] if transpose == 'inner' else None,
                            dimension_names=names, fill_value=fill)
    new.attrs.update(attributes)
    if transpose == 'outer':
        path = root + '/' + node + '/zarr.json'
        meta = json.load(open(path))
        meta['codecs'] = [{'name': 'transpose', 'configuration': {'order': list(swap)}}] + meta['codecs']
        json.dump(meta, open(path, 'w'))
    zarr.open_array(root + '/' + node, mode='r+')[...] = values
zarr.consolidate_metadata(root)"
        RESULT_VARIABLE result)
    if(result)
        message(FATAL_ERROR "sharding ${name} failed: ${result}")
    endif()
endfunction()
sharded_flag(flag_inner "('inner', 32, 16)" "('none', 32, 16)")
run("verify.py with a flag in SKY's shards, transposed only within SKY's chunks" "${SOURCE_DIR}/tools/testset/verify.py"
    "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flag_inner.zarr" --flag yes --chunks 32,16,8 --shards 64,64,8)
sharded_flag(flag_outer "('outer', 32, 16)" "('outer', 32, 16)")
run("verify.py with SKY and its flag sharded after a transpose" "${SOURCE_DIR}/tools/testset/verify.py"
    "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flag_outer.zarr" --flag yes --chunks 32,16,8 --shards 64,64,8)
sharded_flag(flag_across "('outer', 32, 16)" "('outer', 16, 32)")
refused("a flag sharded after a transpose in other chunks than SKY's" "${OUTPUT_DIR}/cube.fits"
        "${OUTPUT_DIR}/flag_across.zarr" --flag yes --chunks 32,16,8 --shards 64,64,8)

# A carta-zarr that applied the mask wrongly -- one masked NaN read back as a finite value -- must be
# caught by the masked comparison, which only a stand-in for the bench can make happen.
file(WRITE "${OUTPUT_DIR}/misapplied.py" "import subprocess, sys
import numpy as np
result = subprocess.run(['${BENCH}'] + sys.argv[1:], capture_output=True)
out = result.stdout
if sys.argv[1] == 'pixels' and '--unmasked' not in sys.argv and not result.returncode:
    pixels = np.frombuffer(out, dtype=np.float32).copy()
    pixels[np.flatnonzero(np.isnan(pixels))[0]] = 123.0
    out = pixels.tobytes()
sys.stdout.buffer.write(out)
sys.stderr.buffer.write(result.stderr)
sys.exit(result.returncode)
")
file(WRITE "${OUTPUT_DIR}/misapplied.sh" "#!/bin/sh\nexec \"${UV}\" run --quiet --no-project --python-preference only-managed --with numpy==2.3.1 python \"${OUTPUT_DIR}/misapplied.py\" \"$@\"\n")
file(CHMOD "${OUTPUT_DIR}/misapplied.sh" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
refused("a pixel mask applied where the FITS cube is not NaN" "${OUTPUT_DIR}/cube.fits" "${OUTPUT_DIR}/flagged.zarr" --flag yes
        --chunks 32,32,8 --bench "${OUTPUT_DIR}/misapplied.sh")

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

# stats.py takes the shape and the chunk as carta-zarr reads them: by name, so m stored before l is the
# same cube; and of the chunk a read decodes, which zarr-python's .chunks does not give past a transpose
# ahead of the shards.
function(stats_of dataset variable)
    execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/stats.py" "${OUTPUT_DIR}/${dataset}.zarr"
        RESULT_VARIABLE result OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(result)
        message(FATAL_ERROR "stats.py ${dataset} failed: ${result}\n${err}")
    endif()
    string(JSON shape GET "${out}" shape)
    string(JSON chunk GET "${out}" chunk)
    string(JSON nan GET "${out}" nan)
    string(JSON channels GET "${out}" nan_channels)
    string(REGEX REPLACE "[ \n]" "" summary "${shape}${chunk} ${nan} ${channels}")
    set(${variable} "${summary}" PARENT_SCOPE)
endfunction()
stats_of(cube plain)
stats_of(reordered reordered)
if(NOT reordered STREQUAL plain)
    message(FATAL_ERROR "stats.py read m stored before l as another cube: ${reordered}, not ${plain}")
endif()
stats_of(flag_outer outer)
if(NOT outer MATCHES "^\\[90,70,40\\]\\[32,16,8\\]")
    message(FATAL_ERROR "stats.py took the shards of a transposed layout for its chunks: ${outer}")
endif()

# zarr-to-fits.py writes the coordinates carta-zarr reads, whatever they are, so the FITS file it writes
# of a rotated, B1950, barycentric or otherwise unusual dataset verifies against that dataset; and it
# refuses one a FITS header cannot hold.
foreach(name rotated equinox observer pole tan gigahertz)
    run("zarr-to-fits.py ${name}" "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${OUTPUT_DIR}/${name}.zarr"
        "${OUTPUT_DIR}/${name}.fits" --block-mib 1)
    run("verify.py ${name} against its own FITS file" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/${name}.fits"
        "${OUTPUT_DIR}/${name}.zarr" --flag no --chunks 32,32,8)
endforeach()
# As many channels as the cigar has, evenly spaced: the increment taken from the first two would be off
# by their rounding 30,000 times over, and the cigar refused.
run("generate.py long" "${SOURCE_DIR}/tools/zarr-bench/generate.py" --synthetic --shape frequency=30000,polarization=1,l=8,m=8
    --chunk l=8,m=8,frequency=1000 --workers 1 --output "${OUTPUT_DIR}/long.zarr")
run("zarr-to-fits.py long" "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${OUTPUT_DIR}/long.zarr" "${OUTPUT_DIR}/long.fits")
run("verify.py long" "${SOURCE_DIR}/tools/testset/verify.py" "${OUTPUT_DIR}/long.fits" "${OUTPUT_DIR}/long.zarr" --flag no
    --chunks 8,8,1000)
changed(tai "t = zarr.open_array(root + '/time', mode='r+')
t.attrs['scale'] = 'tai'
a = zarr.open_array(root + '/SKY', mode='r+')
date = dict(a.attrs['obsdate'])
date['attrs'] = dict(date['attrs'], scale='tai')
a.attrs['obsdate'] = date")
foreach(case "middle|evenly spaced" "widebeam|plane to plane" "tai|wrong instant" "halfbeam|one on each")
    string(REPLACE "|" ";" case "${case}")
    list(GET case 0 name)
    list(GET case 1 reason)
    execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${OUTPUT_DIR}/${name}.zarr"
                            "${OUTPUT_DIR}/${name}.fits"
        RESULT_VARIABLE result OUTPUT_QUIET ERROR_VARIABLE err)
    if(NOT result OR NOT err MATCHES "${reason}")
        message(FATAL_ERROR "zarr-to-fits.py wrote ${name}, which a FITS header cannot hold: ${result}\n${err}")
    endif()
endforeach()

# zarr-to-fits.py refuses a dataset with a flag: FITS has only NaN to mark a pixel with, and a flagged
# dataset may hold finite values under its flag.
run("generate.py --flag" "${SOURCE_DIR}/tools/zarr-bench/generate.py" --synthetic --flag --shape frequency=8,polarization=1,l=20,m=20
    --chunk l=10,m=10,frequency=4 --workers 1 --output "${OUTPUT_DIR}/masked.zarr")
execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${OUTPUT_DIR}/masked.zarr"
                        "${OUTPUT_DIR}/masked.fits"
    RESULT_VARIABLE result OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT result OR NOT err MATCHES "carries a flag")
    message(FATAL_ERROR "zarr-to-fits.py wrote a flagged dataset: ${result}\n${err}")
endif()
# Nor one whose flag is declared by any data group, wherever it is kept.
set(nested "${OUTPUT_DIR}/nested.zarr")
file(COPY "${cube}/" DESTINATION "${nested}")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
            python -c "import json, numpy as np, zarr
root = '${nested}'
sky = zarr.open_array(root + '/SKY', mode='r')
zarr.create_array(root + '/masks/FLAG_SKY', shape=sky.shape, chunks=sky.chunks, dtype=bool,
                  dimension_names=list(sky.metadata.dimension_names), fill_value=False)
meta = json.load(open(root + '/zarr.json'))
meta['attributes']['data_groups']['other'] = {'sky': 'SKY', 'flag': 'masks/FLAG_SKY'}
json.dump(meta, open(root + '/zarr.json', 'w'))"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "nesting a flag failed: ${result}")
endif()
execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${nested}" "${OUTPUT_DIR}/nested.fits"
    RESULT_VARIABLE result OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT result OR NOT err MATCHES "carries a flag")
    message(FATAL_ERROR "zarr-to-fits.py wrote a dataset whose flag another group declares: ${result}\n${err}")
endif()
# Nor one with a flag typed as one and declared by nothing, kept behind a directory link.
file(COPY "${cube}/" DESTINATION "${OUTPUT_DIR}/linked.zarr")
execute_process(
    COMMAND "${UV}" run --quiet --no-project --python-preference only-managed --with zarr==3.2.1 --with numpy==2.3.1
            python -c "import zarr
sky = zarr.open_array('${cube}/SKY', mode='r')
flag = zarr.create_array('${OUTPUT_DIR}/linked_flag', shape=sky.shape, chunks=sky.chunks, dtype=bool,
                         dimension_names=list(sky.metadata.dimension_names), fill_value=False)
flag.attrs['type'] = 'flag'"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "writing the linked flag failed: ${result}")
endif()
file(CREATE_LINK "${OUTPUT_DIR}/linked_flag" "${OUTPUT_DIR}/linked.zarr/LINKED" SYMBOLIC)
execute_process(COMMAND "${UV}" run --quiet --script "${SOURCE_DIR}/tools/testset/zarr-to-fits.py" "${OUTPUT_DIR}/linked.zarr"
                        "${OUTPUT_DIR}/linked.fits"
    RESULT_VARIABLE result OUTPUT_QUIET ERROR_VARIABLE err)
if(NOT result OR NOT err MATCHES "carries a flag")
    message(FATAL_ERROR "zarr-to-fits.py wrote a dataset with an undeclared flag behind a link: ${result}\n${err}")
endif()
