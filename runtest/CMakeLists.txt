set(SRCS
    runtest.cc)

add_executable(runtest
    ${SRCS})
target_link_libraries(runtest
    sgp4)

set(SRCS_ISO8601
    runtest_iso8601.cc)

add_executable(runtest_iso8601
    ${SRCS_ISO8601})
target_link_libraries(runtest_iso8601
    sgp4)
