set(SRCS
    sattrack.cpp)

add_executable(sattrack
    ${SRCS})
target_link_libraries(sattrack
    sgp4
    rt)
