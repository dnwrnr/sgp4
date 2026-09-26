add_executable(csvprop csvprop.cpp)
target_link_libraries(csvprop PRIVATE SGP4::sgp4)
get_filename_component(CSV_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data" ABSOLUTE)
target_compile_definitions(csvprop PRIVATE CSV_DATA_DIR="${CSV_DATA_DIR}")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(csvprop PRIVATE stdc++fs)
endif()
