# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

# The notices installed beside the library name every library linked into it.
#
# libcarta-zarr carries TensorStore and the libraries it brings in, statically, so a binary
# distribution of it is a distribution of theirs and has to carry their notices. The list of them is
# written by hand in cmake/ThirdPartyNotices.cmake; this reads the libraries the link line actually
# names out of the build tree and fails on any whose notice the generated file does not hold, so that
# a dependency TensorStore adds is not shipped without one.
#
# Takes BINARY_DIR, the build tree; DEPS_DIR, where the fetched libraries were built, which is the
# build tree's _deps unless FETCHCONTENT_BASE_DIR moved it; and NOTICES, the generated file.
cmake_minimum_required(VERSION 3.24)

if(NOT EXISTS "${NOTICES}")
    message(FATAL_ERROR "No third-party notices were generated at ${NOTICES}")
endif()
file(READ "${NOTICES}" notices)

# The link line, from whichever generator wrote the tree.
set(link "")
file(GLOB_RECURSE link_files "${BINARY_DIR}/CMakeFiles/carta_zarr.dir/link.txt")
if(link_files)
    list(GET link_files 0 link_file)
    file(READ "${link_file}" link)
elseif(EXISTS "${BINARY_DIR}/build.ninja")
    file(STRINGS "${BINARY_DIR}/build.ninja" ninja)
    set(in_library FALSE)
    foreach(line IN LISTS ninja)
        if(line MATCHES "^build [^:]*libcarta-zarr[^:]*: ")
            set(in_library TRUE)
        elseif(in_library AND line MATCHES "^  LINK_LIBRARIES = (.*)")
            set(link "${CMAKE_MATCH_1}")
            break()
        elseif(in_library AND NOT line MATCHES "^  ")
            set(in_library FALSE)
        endif()
    endforeach()
endif()
if(link STREQUAL "")
    message(FATAL_ERROR "Could not find how libcarta-zarr is linked in ${BINARY_DIR}")
endif()

# A library built in the tree is named relative to it, as _deps/<name>-build/; one built where
# FETCHCONTENT_BASE_DIR put it, as an offline build does, by its absolute path. Both are read alike.
string(REPLACE "${DEPS_DIR}/" "_deps/" link "${link}")
string(REGEX MATCHALL "_deps/[A-Za-z0-9_.+-]+-build/" built "${link}")
list(REMOVE_DUPLICATES built)
if(NOT built)
    message(FATAL_ERROR "The link line names no fetched library, which this check cannot be right about")
endif()
set(missing "")
foreach(entry IN LISTS built)
    string(REGEX REPLACE "^_deps/(.*)-build/$" "\\1" name "${entry}")
    string(FIND "${notices}" "== ${name} ==" at)
    if(at EQUAL -1)
        list(APPEND missing "${name}")
    endif()
endforeach()
if(missing)
    message(FATAL_ERROR "Linked into libcarta-zarr without a notice: ${missing}. "
                        "Add each to cmake/ThirdPartyNotices.cmake.")
endif()
list(LENGTH built count)
message(STATUS "Every one of the ${count} fetched libraries linked in has its notice")
