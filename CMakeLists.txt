CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(SGP4)

if (POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
endif()

macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        include(CheckCXXCompilerFlag)
        CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
        if(COMPILER_SUPPORTS_CXX11)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        else()
            message(STATUS "Unable to determine C++11 option for compiler ${CMAKE_CXX_COMPILER}.")
        endif()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)

use_cxx11()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-align")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion")
endif(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

include_directories(libsgp4)

add_subdirectory(libsgp4)
add_subdirectory(sattrack)
add_subdirectory(runtest)
add_subdirectory(passpredict)

file(COPY SGP4-VER.TLE DESTINATION ${PROJECT_BINARY_DIR})
