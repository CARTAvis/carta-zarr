cmake_minimum_required(VERSION 3.24)
project(carta_zarr_examples LANGUAGES CXX)

find_package(CartaZarr 0.2 CONFIG REQUIRED)

foreach(example first_read statistics read_control animation)
    add_executable(carta_zarr_${example} ${example}.cc)
    target_compile_features(carta_zarr_${example} PRIVATE cxx_std_17)
    target_link_libraries(carta_zarr_${example} PRIVATE CARTA::zarr)
endforeach()
