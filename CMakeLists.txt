CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(SGP4)

file(GLOB SRCS libsgp4/*.cpp)

include_directories(libsgp4)

add_library(sgp4 STATIC ${SRCS})
add_executable(sattrack sattrack/sattrack.cpp)
target_link_libraries(sattrack sgp4 rt)
add_executable(runtest runtest/runtest.cpp)
target_link_libraries(runtest sgp4 rt)
add_executable(passpredict passpredict/passpredict.cpp)
target_link_libraries(passpredict sgp4 rt)

file(COPY SGP4-VER.TLE DESTINATION ${PROJECT_BINARY_DIR})
