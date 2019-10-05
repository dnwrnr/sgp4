set(SRCS
    isstrack.cc)

add_executable(isstrack
    ${SRCS})
target_link_libraries(isstrack
    curl
    sgp4
    rt)

