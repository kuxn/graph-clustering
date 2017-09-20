# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
message("tests/. CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message("tests/. CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("tests/. CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR})
message("tests/. CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})
message("tests/. graph-clustering_BINARY_DIR: " ${graph-clustering_BINARY_DIR})
message("tests/. graph-clustering_SOURCE_DIR: " ${graph-clustering_SOURCE_DIR})
if(result)
  message("CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR})
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Only build gtest
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
message("Add subdirectory, CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
message("<gtest>/. gtest_SOURCE_DIR: " ${gtest_SOURCE_DIR})
message("<gtest>/. gtest_BINARY_DIR: " ${gtest_BINARY_DIR})
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

set(ENV{TEST_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/dotfiles)
message("set environment variable: "  $ENV{TEST_FILES})
