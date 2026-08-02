add_executable(sgp4_tests
    test_util.cc
    test_timespan.cc
    test_vector.cc
    test_datetime.cc
    test_tle.cc
    test_sgp4.cc
    test_observer.cc
    test_csv_loader.cc)

target_link_libraries(sgp4_tests PRIVATE sgp4 gtest_main)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(sgp4_tests PRIVATE stdc++fs)
endif()

gtest_discover_tests(sgp4_tests)
