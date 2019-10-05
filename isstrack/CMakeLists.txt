set(SRCS
    isstrack.cc)

add_executable(isstrack
    ${SRCS})
target_link_libraries(isstrack
    sgp4
    rt)

