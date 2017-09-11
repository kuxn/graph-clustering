cmake_minimum_required(VERSION 3.5.2)
project(graph-clustering C CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON) # disable in source build

set(CMAKE_CXX_STANDARD 11) # enable C++11 standard

# -- Handle output

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -- Set common parts
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common/include)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/common/src COMMON_SOURCE_FILES)

# -- Subdirectory that has CMakeLists

message("/. CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message("/. CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("/. CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR})
message("/. CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})
message("/. graph-clustering_BINARY_DIR: " ${graph-clustering_BINARY_DIR})
message("/. graph-clustering_SOURCE_DIR: " ${graph-clustering_SOURCE_DIR})

# -- For add_test in the subdirectory, $ make test
enable_testing()

add_subdirectory(serial)
add_subdirectory(parallel)
add_subdirectory(tests)
