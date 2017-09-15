# -- Include

include_directories(${PARALLEL_INCLUDE_DIR}) # add given directories to current list of directories

# -- Source

# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SOURCE_FILES)

# -- Packages
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

find_package(Boost 1.58 REQUIRED COMPONENTS mpi serialization timer program_options)
include_directories(${Boost_INCLUDE_DIRS})

# -- Flags

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# -- Libs

add_library(parallel_core ${COMMON_SOURCE_FILES} ${PARALLEL_SOURCE_FILES})

# -- Binary

add_executable(main_parallel ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc)
target_link_libraries(main_parallel ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} parallel_core)

# -- Add gtests to ctest

add_gtest(ParallelTest ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} parallel_core)

# -- Options

option(Median_ "Partition the graph by comparing the elements with the median of the eigenvector instead of Fielder method" OFF)  # display with $ cmake -LAH
if(Median_)
    message("-- Median method is enabled")
    target_compile_definitions(main_parallel PUBLIC -DMedian_)  # cmake .. -DMedian_=ON
endif()
