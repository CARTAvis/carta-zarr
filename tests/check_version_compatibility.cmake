# A version says which builds can stand in for one another, and the soname says the same thing.
#
# Before 1.0 the minor version is the one that breaks: 0.2 may change the ABI 0.1 had, so a consumer
# that asked for 0.1 must not be handed 0.2 by find_package, and a binary linked against 0.1 must not
# load 0.2 through the same soname. From 1.0 on it is the major. See docs/adr/0020.
#
# Takes CONFIG_VERSION, the generated CartaZarrConfigVersion.cmake; VERSION, the project's version;
# and SONAME, the library's soname file name.
cmake_minimum_required(VERSION 3.24)

if(NOT EXISTS "${CONFIG_VERSION}")
    message(FATAL_ERROR "No package version file was generated at ${CONFIG_VERSION}")
endif()

string(REPLACE "." ";" parts "${VERSION}")
list(GET parts 0 major)
list(GET parts 1 minor)
list(GET parts 2 patch)

# What the version file answers for one requested version. Each call starts clean, as find_package
# does, and the file is read again so nothing one answer set carries into the next.
function(answer requested result)
    string(REPLACE "." ";" asked "${requested}")
    list(LENGTH asked count)
    set(PACKAGE_FIND_VERSION "${requested}")
    set(PACKAGE_FIND_VERSION_COUNT ${count})
    set(PACKAGE_FIND_VERSION_MAJOR 0)
    set(PACKAGE_FIND_VERSION_MINOR 0)
    set(PACKAGE_FIND_VERSION_PATCH 0)
    list(GET asked 0 PACKAGE_FIND_VERSION_MAJOR)
    if(count GREATER 1)
        list(GET asked 1 PACKAGE_FIND_VERSION_MINOR)
    endif()
    if(count GREATER 2)
        list(GET asked 2 PACKAGE_FIND_VERSION_PATCH)
    endif()
    unset(PACKAGE_VERSION_COMPATIBLE)
    include("${CONFIG_VERSION}")
    set(${result} "${PACKAGE_VERSION_COMPATIBLE}" PARENT_SCOPE)
endfunction()

function(require requested expected why)
    answer("${requested}" compatible)
    if(NOT "${compatible}" STREQUAL "${expected}")
        message(FATAL_ERROR "find_package(CartaZarr ${requested}) against ${VERSION}: "
            "expected compatible=${expected}, got '${compatible}' -- ${why}")
    endif()
endfunction()

math(EXPR next_minor "${minor} + 1")
math(EXPR next_major "${major} + 1")

require("${major}.${minor}" TRUE "the version a consumer was written against is found")
require("${major}.${minor}.${patch}" TRUE "and so is its exact release")
require("${major}.${next_minor}" FALSE "a newer minor than this one is not this one")
require("${next_major}.0" FALSE "nor is a newer major")

if(major EQUAL 0)
    if(minor GREATER 0)
        math(EXPR previous_minor "${minor} - 1")
        require("0.${previous_minor}" FALSE "before 1.0 an older minor may have had another ABI")
        require("0" FALSE "and a request for 0 alone asks for 0.0, which is one of them")
    endif()
    set(expected_soname_version "${major}.${minor}")
else()
    if(minor GREATER 0)
        require("${major}.0" TRUE "from 1.0 on an older minor of the same major is kept")
    endif()
    require("${major}" TRUE "and a request for the major alone takes any release of it")
    set(expected_soname_version "${major}")
endif()

# libcarta-zarr.so.0.1 on Linux, libcarta-zarr.0.1.dylib on macOS: the version sits right after the
# name either way, and nothing more of it than says which ABI this is.
if(NOT SONAME MATCHES "^libcarta-zarr[.](so[.])?([0-9.]*[0-9])([.]dylib)?$")
    message(FATAL_ERROR "Cannot read a version out of the soname '${SONAME}'")
endif()
if(NOT CMAKE_MATCH_2 STREQUAL expected_soname_version)
    message(FATAL_ERROR "The soname '${SONAME}' carries version ${CMAKE_MATCH_2}; "
        "${VERSION} should carry ${expected_soname_version}")
endif()
