set(SOURCE_LIST runtest.cpp)

add_executable(runtest ${SOURCE_LIST})

set_target_properties(runtest PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

target_compile_features(runtest PRIVATE cxx_std_17)
target_compile_options(runtest PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-pedantic;-Wno-long-long;-Wcast-align;-Wsign-conversion>
)
target_link_libraries(runtest PRIVATE libsgp4::sgp4)

add_test(NAME runtest COMMAND runtest)

if(WIN32 AND SGP4_BUILD_SHARED_LIBS)
  # copy sgp4 dll into the tests directory
  add_custom_command(TARGET runtest POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:libsgp4::sgp4>" "$<TARGET_FILE_DIR:runtest>")
endif()
file(COPY SGP4-VER.TLE DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
