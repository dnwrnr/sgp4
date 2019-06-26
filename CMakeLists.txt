CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(SGP4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_COMPILER_IS_GNUCC)
    add_compile_options(-Wall)
    add_compile_options(-Wextra)
    add_compile_options(-Wshadow)
    add_compile_options(-Wnon-virtual-dtor)
    add_compile_options(-pedantic)
    add_compile_options(-Wno-long-long)
    add_compile_options(-Wcast-align)
    add_compile_options(-Wsign-conversion)
endif()

include_directories(libsgp4)

add_subdirectory(libsgp4)
add_subdirectory(sattrack)
add_subdirectory(runtest)
add_subdirectory(passpredict)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/SGP4-VER.TLE DESTINATION ${PROJECT_BINARY_DIR})
