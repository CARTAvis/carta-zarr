# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

if(NOT DEFINED CARTA_ZARR_BINARY_DIR OR NOT DEFINED CARTA_ZARR_SOURCE_DIR OR NOT DEFINED CARTA_ZARR_PACKAGE_TEST_DIR)
    message(FATAL_ERROR "Package consumer test variables are not set")
endif()

set(prefix "${CARTA_ZARR_PACKAGE_TEST_DIR}/prefix")
set(build_dir "${CARTA_ZARR_PACKAGE_TEST_DIR}/build")
file(REMOVE_RECURSE "${CARTA_ZARR_PACKAGE_TEST_DIR}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" --install "${CARTA_ZARR_BINARY_DIR}" --prefix "${prefix}"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "carta-zarr install failed: ${result}")
endif()
file(GLOB_RECURSE installed_notices "${prefix}/share/doc/*/THIRD_PARTY_NOTICES")
if(NOT installed_notices)
    message(FATAL_ERROR "The install carries no THIRD_PARTY_NOTICES for what the library links in")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}"
        -S "${CARTA_ZARR_SOURCE_DIR}/tests/consumer"
        -B "${build_dir}"
        -DCMAKE_PREFIX_PATH=${prefix}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_BUILD_RPATH=${prefix}/lib
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "External consumer configure failed: ${result}")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${build_dir}"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "External consumer build failed: ${result}")
endif()

execute_process(
    COMMAND "${build_dir}/carta_zarr_consumer"
    RESULT_VARIABLE result)
if(result)
    message(FATAL_ERROR "External consumer execution failed: ${result}")
endif()

foreach(example first_read statistics read_control animation)
    # Cover a synthetic store, XRADIO's own conversion and the single-channel
    # image whose width supplies the spectral increment.
    foreach(fixture pixels conformance conformance_continuum)
        execute_process(
            COMMAND "${build_dir}/examples/carta_zarr_${example}"
                    "${CARTA_ZARR_SOURCE_DIR}/tests/data/images/zarr/xradio/${fixture}"
            RESULT_VARIABLE result)
        if(result)
            message(FATAL_ERROR "Documented ${example} example failed on ${fixture}: ${result}")
        endif()
    endforeach()

    # Check the documented command-line contract, including the error path.
    execute_process(
        COMMAND "${build_dir}/examples/carta_zarr_${example}"
        RESULT_VARIABLE result)
    if(NOT result EQUAL 1)
        message(FATAL_ERROR "Documented ${example} example accepted missing arguments: ${result}")
    endif()
endforeach()
