# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
message("<gtest>/. CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message("<gtest>/. CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("<gtest>/. CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR})
message("<gtest>/. CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})
message("<gtest>/. graph-clustering_BINARY_DIR: " ${graph-clustering_BINARY_DIR})
message("<gtest>/. graph-clustering_SOURCE_DIR: " ${graph-clustering_SOURCE_DIR})
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

# Now simply link against gtest or gtest_main as needed. Eg
# enable_testing() // If added here, `ctest -N` can only be invoked in the directory
function(add_gtest testName)
    message("adding test ", ${testName})
    add_executable(${testName}.gtest ${ARGN})
    target_link_libraries(${testName}.gtest gtest_main)
    add_test(${testName}_test ${CMAKE_BINARY_DIR}/bin/${testName}.gtest)
#add_test(NAME example_test COMMAND example.gtest)
endfunction()

add_gtest(example example.cpp)
add_gtest(example2 example2.cpp)
add_gtest(example3 example3.cpp)
