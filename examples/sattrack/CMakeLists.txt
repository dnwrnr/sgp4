set(SOURCE_LIST sattrack.cpp)

add_executable(sattrack ${SOURCE_LIST})

set_target_properties(sattrack PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

target_compile_features(sattrack PRIVATE cxx_std_17)
target_compile_options(sattrack PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-pedantic;-Wno-long-long;-Wcast-align;-Wsign-conversion>
)
target_link_libraries(sattrack PRIVATE libsgp4::sgp4)

if(WIN32 AND SGP4_BUILD_SHARED_LIBS)
  # copy sgp4 dll into the examples/sattrack directory
  add_custom_command(TARGET sattrack POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:libsgp4::sgp4>" "$<TARGET_FILE_DIR:sattrack>")
endif()
