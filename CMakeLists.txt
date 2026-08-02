cmake_minimum_required(VERSION 3.14)

project(SGP4
    VERSION 2.0.0
    DESCRIPTION "SGP4 satellite orbit propagation library"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(BUILD_EXAMPLES "Build example programs (sattrack, runtest, passpredict)" ON)
option(BUILD_TESTS "Build unit tests" ON)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0)
FetchContent_MakeAvailable(googletest)

enable_testing()
include(GoogleTest)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/SGP4-VER.TLE"
    "${CMAKE_CURRENT_BINARY_DIR}/SGP4-VER.TLE"
    COPYONLY)

add_subdirectory(libsgp4)

if(BUILD_EXAMPLES)
    add_subdirectory(sattrack)
    add_subdirectory(runtest)
    add_subdirectory(passpredict)
    add_subdirectory(csvprop)
endif()

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format clang-format-21 clang-format-20 clang-format-19 clang-format-18)

if(CLANG_FORMAT_EXECUTABLE)
    file(GLOB_RECURSE SGP4_FORMAT_SOURCES
        "${PROJECT_SOURCE_DIR}/libsgp4/*.h"
        "${PROJECT_SOURCE_DIR}/libsgp4/*.cc"
        "${PROJECT_SOURCE_DIR}/tests/*.cc"
        "${PROJECT_SOURCE_DIR}/sattrack/*.cc"
        "${PROJECT_SOURCE_DIR}/runtest/*.cc"
        "${PROJECT_SOURCE_DIR}/passpredict/*.cc"
        "${PROJECT_SOURCE_DIR}/csvprop/*.cc")

    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -i ${SGP4_FORMAT_SOURCES}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Applying clang-format to all source files"
        USES_TERMINAL)

    add_custom_target(format-check
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -i --dry-run --Werror ${SGP4_FORMAT_SOURCES}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Checking clang-format compliance"
        USES_TERMINAL)
endif()
