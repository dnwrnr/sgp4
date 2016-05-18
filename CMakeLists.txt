CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(SGP4)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

file(GLOB SRCS libsgp4/*.cpp)

include_directories(libsgp4)

add_library(sgp4 STATIC ${SRCS})
add_executable(sattrack sattrack/sattrack.cpp)
target_link_libraries(sattrack sgp4)
add_executable(runtest runtest/runtest.cpp)
target_link_libraries(runtest sgp4)
add_executable(passpredict passpredict/passpredict.cpp)
target_link_libraries(passpredict sgp4)

file(COPY SGP4-VER.TLE DESTINATION ${PROJECT_BINARY_DIR})
