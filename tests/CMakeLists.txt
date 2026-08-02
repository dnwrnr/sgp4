add_executable(sgp4_tests
    test_util.cpp
    test_timespan.cpp
    test_vector.cpp
    test_datetime.cpp
    test_tle.cpp
    test_sgp4.cpp
    test_observer.cpp
    test_csv_loader.cpp)

target_link_libraries(sgp4_tests PRIVATE sgp4 gtest_main)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(sgp4_tests PRIVATE stdc++fs)
endif()

gtest_discover_tests(sgp4_tests)
