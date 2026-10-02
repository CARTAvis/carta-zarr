# Writes the datasets tests/bench_generator_test.cc reads, with tools/zarr-bench/generate.py.
#
# Three synthetic cubes of one seed in three layouts -- plain zstd chunks; blosc in shards, without
# consolidated metadata; and a flag variable where the others write NaN -- and two rewrites of a
# committed fixture, cropped: one in a layout of its own, one in the fixture's. The test reads them back and checks that the layouts hold the same
# pixels, which is the property every comparison between layouts rests on.

if(NOT DEFINED UV OR NOT DEFINED SOURCE_DIR OR NOT DEFINED OUTPUT_DIR)
    message(FATAL_ERROR "UV, SOURCE_DIR and OUTPUT_DIR must be set")
endif()

set(generator "${SOURCE_DIR}/tools/zarr-bench/generate.py")
set(shape "frequency=12,polarization=2,l=150,m=130")
file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

function(generate name)
    execute_process(
        COMMAND "${UV}" run --quiet --script "${generator}" --output "${OUTPUT_DIR}/${name}" --workers 2 ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_QUIET)
    if(result)
        message(FATAL_ERROR "generate.py failed writing ${name}: ${result}")
    endif()
endfunction()

generate(plain --synthetic --shape ${shape} --chunk l=64,m=64,frequency=4 --codec zstd:3)
generate(sharded --synthetic --shape ${shape} --chunk l=32,m=128,frequency=2,polarization=2
         --shard l=128,m=128,frequency=8 --codec blosc:lz4:5:bitshuffle --no-consolidate --block-mib 1)
generate(flagged --synthetic --shape ${shape} --chunk l=64,m=64,frequency=4 --flag)
generate(rewritten --source "${SOURCE_DIR}/tests/data/images/zarr/xradio/pixels" --crop l=1:4
         --chunk l=3,m=2,frequency=2 --codec gzip:1)
generate(current --source "${SOURCE_DIR}/tests/data/images/zarr/xradio/pixels_wide" --crop frequency=1:3
         --layout-from-source)

# --force replaces whatever is at --output, so an --output that is the source, an alias of it, or a
# directory around it would delete the source before a byte of it was read. Each is refused, and the
# source is still there afterwards. A copy of a fixture stands in for the source, so that a
# regression costs the build tree rather than the repository.
set(overlap "${OUTPUT_DIR}/overlap")
set(victim "${overlap}/victim.zarr")
file(MAKE_DIRECTORY "${overlap}")
file(COPY "${SOURCE_DIR}/tests/data/images/zarr/xradio/pixels/" DESTINATION "${victim}")
file(CREATE_LINK "${victim}" "${overlap}/alias.zarr" SYMBOLIC)

function(refuse output)
    execute_process(
        COMMAND "${UV}" run --quiet --script "${generator}" --source "${victim}" --output "${output}" --force
                --workers 1 --layout-from-source
        RESULT_VARIABLE result
        OUTPUT_QUIET ERROR_QUIET)
    if(NOT result)
        message(FATAL_ERROR "generate.py accepted --output ${output} over its own --source")
    endif()
    if(NOT EXISTS "${victim}/zarr.json")
        message(FATAL_ERROR "generate.py deleted its --source writing --output ${output}")
    endif()
endfunction()

refuse("${victim}")
refuse("${overlap}/alias.zarr")
refuse("${overlap}")
refuse("${victim}/inside")
