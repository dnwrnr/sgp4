CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)
PROJECT(SGP4)

if (POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "..." FORCE)
endif()

set(CMAKE_INSTALL_LIBDIR "/lib64")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enable AVX2 optimizations (without -march=native to avoid ISA mismatches)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -fopenmp")
endif()

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-align")
endif(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(libsgp4)

add_subdirectory(libsgp4)
add_subdirectory(sattrack)
add_subdirectory(runtest)
add_subdirectory(passpredict)

# Optional test programs not present in this repo snapshot are skipped.
# add_executable(residual_trace residual_trace.cc)
# target_link_libraries(residual_trace sgp4)
#
# add_executable(test_bstar test_bstar.cc)
# target_link_libraries(test_bstar sgp4)
#
# add_executable(testbench testbench.cc)
# target_link_libraries(testbench sgp4)
#
add_executable(full_test full_test.cc)
target_link_libraries(full_test sgp4)

add_executable(benchmark benchmark.cc)
target_link_libraries(benchmark sgp4)

# Optional test programs not present in this repo snapshot are skipped.
# add_executable(baseline_test baseline_test.cc)
# target_link_libraries(baseline_test sgp4)
#
# add_executable(calibrate_k calibrate_k.cc)
# target_link_libraries(calibrate_k sgp4)
#
# add_executable(oracle_test oracle_test.cc)
# target_link_libraries(oracle_test sgp4)
#
# add_executable(piecewise_test piecewise_test.cc)
# target_link_libraries(piecewise_test sgp4)
#
# add_executable(multiyear_bench multiyear_bench.cc)
# target_link_libraries(multiyear_bench sgp4)
#
# add_executable(trace_cpp trace_cpp.cc)
# target_link_libraries(trace_cpp sgp4)

file(COPY SGP4-VER.TLE DESTINATION ${PROJECT_BINARY_DIR})
