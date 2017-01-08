set(SRCS
    runtest.cpp)

add_executable(runtest
    ${SRCS})
target_link_libraries(runtest
    sgp4
    rt)
