cmake_minimum_required(VERSION 3.5.2)
project(graph-clustering C CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON) # disable in source build

set(CMAKE_CXX_STANDARD 11) # enable C++11 standard

# -- Handle output

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -- Set common parts

set(COMMON_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/common/include)
set(SERIAL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/serial/include)
set(PARALLEL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/parallel/include)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/common/src COMMON_SOURCE_FILES)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/serial/src SERIAL_SOURCE_FILES)
list(REMOVE_ITEM SERIAL_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/serial/src/main.cc)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/parallel/src PARALLEL_SOURCE_FILES)
list(REMOVE_ITEM PARALLEL_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/parallel/src/main.cc)
include_directories(${COMMON_INCLUDE_DIR})

# -- For add_test in the subdirectory, $ make test

enable_testing()
function(add_gtest testName)
    add_executable(${testName}.gtest ${CMAKE_SOURCE_DIR}/tests/${testName}.cc)
    target_link_libraries(${testName}.gtest gtest_main ${ARGN})
    if (testName MATCHES "(.*)Serial(.*)")
        message("-- Added serial test" ${testName})
        add_test(${testName} ${CMAKE_BINARY_DIR}/bin/${testName}.gtest)
    else()
        message("-- Added parallel test" ${testName})
        add_test(ParallelTest ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_BINARY_DIR}/bin/${testName}.gtest)
    endif()
endfunction()

# -- Subdirectory that has CMakeLists.txt

add_subdirectory(serial)
add_subdirectory(parallel)
add_subdirectory(tests)
