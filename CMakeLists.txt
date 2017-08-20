cmake_minimum_required(VERSION 3.7.1)
project(graph-clustering C CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON) # disable in source build

set(CMAKE_CXX_STANDARD 11) # enable C++11 standard

# -- Handle output

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -- Subdirectory that has CMakeLists

add_subdirectory(serial)
add_subdirectory(parallel)

