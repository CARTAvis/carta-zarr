# The notices of what libcarta-zarr carries.
#
# TensorStore and the libraries it brings in are linked into the library statically, so a binary
# distribution of libcarta-zarr is a distribution of each of them, and must carry their notices. These
# are they: the ten the link line names, and half, a header-only library TensorStore compiles in for
# float16. TensorStore declares some forty more that are downloaded and never compiled into this
# library; they are not here. tests/check_third_party_notices.cmake fails when the link line names a
# fetched library this list does not.
#
# Each entry is the name FetchContent fetched it under, the name it goes by, and its notice files,
# relative to its source and separated by commas. Where a project offers a choice of licence -- zstd
# is BSD or GPLv2, lz4's library BSD while its tools are GPL -- the file is the one for the library as
# it is used here.
set(CARTA_ZARR_THIRD_PARTY_NOTICES
    "tensorstore|TensorStore|LICENSE"
    "absl|Abseil|LICENSE"
    "riegeli|Riegeli|LICENSE"
    "re2|RE2|LICENSE"
    "zstd|Zstandard|LICENSE"
    "zlib|zlib|LICENSE"
    "blosc|c-blosc|LICENSE.txt,LICENSES/BITSHUFFLE.txt,LICENSES/FASTLZ.txt,LICENSES/LZ4.txt,LICENSES/SNAPPY.txt,LICENSES/STDINT.txt,LICENSES/ZLIB-NG.txt,LICENSES/ZLIB.txt"
    "snappy|Snappy|COPYING"
    "lz4|LZ4|lib/LICENSE"
    "nlohmann_json|JSON for Modern C++|LICENSE.MIT"
    "half|half|LICENSE.txt")

# Writes the notices of every entry fetched into `deps` to `output`. One that was not fetched -- one a
# build takes from the system instead -- is not this library's to carry, and is left out. One that was
# fetched without the notice it is listed with is an error: a distribution that leaves a notice out
# is not one anybody may make, and the list being wrong is the likely reason.
function(carta_zarr_write_third_party_notices deps output)
    string(CONCAT text "Third-party notices for carta-zarr\n\n"
                       "libcarta-zarr contains the following libraries, linked into it statically. Each is\n"
                       "distributed under the licence that follows its name.\n")
    foreach(entry IN LISTS CARTA_ZARR_THIRD_PARTY_NOTICES)
        string(REPLACE "|" ";" fields "${entry}")
        list(GET fields 0 name)
        list(GET fields 1 title)
        list(GET fields 2 files)
        set(source "${deps}/${name}-src")
        if(NOT IS_DIRECTORY "${source}")
            continue()
        endif()
        string(APPEND text "\n\n== ${name} ==\n${title}\n")
        string(REPLACE "," ";" files "${files}")
        foreach(file IN LISTS files)
            if(NOT EXISTS "${source}/${file}")
                message(FATAL_ERROR "${title} was fetched without ${file}, the notice "
                                    "cmake/ThirdPartyNotices.cmake lists for it")
            endif()
            file(READ "${source}/${file}" notice)
            string(APPEND text "\n--- ${file} ---\n\n${notice}")
        endforeach()
    endforeach()
    file(WRITE "${output}" "${text}")
endfunction()
