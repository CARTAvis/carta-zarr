cmake_minimum_required(VERSION 3.24)
project(carta_zarr_examples LANGUAGES CXX)

find_package(CartaZarr 0.1 CONFIG REQUIRED)

add_executable(carta_zarr_first_read first_read.cc)
target_compile_features(carta_zarr_first_read PRIVATE cxx_std_17)
target_link_libraries(carta_zarr_first_read PRIVATE CARTA::zarr)
