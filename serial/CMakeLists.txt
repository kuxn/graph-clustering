# -- Include

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include) # add given directories to current list of directories

# -- Source

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SOURCE_FILES)

# -- Packages

find_package(Boost 1.58 REQUIRED COMPONENTS timer program_options)
include_directories(${Boost_INCLUDE_DIRS})

# -- Flags

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# -- Binary

add_executable(main_serial ${COMMON_SOURCE_FILES} ${SOURCE_FILES})
target_link_libraries(main_serial ${Boost_LIBRARIES})

# -- Options

option(Median_ "Partition the graph by comparing the elements with the median of the eigenvector instead of Fielder method" OFF)  # display with $ cmake -LAH
if(Median_)
    message("-- Median method is enabled")
    target_compile_definitions(main_serial PUBLIC -DMedian_)  # cmake .. -DMedian_=ON
endif()

# -- Profiling

option(profile "Enable profiling" OFF)
if (profile)
    find_library(gperftools profiler)
    if (gperftools)
        message("-- Found Google performance tools")
        target_link_libraries(main_serial profiler)
    else()
        message("-- Google performance tools NOT Found")
    endif()
endif()

