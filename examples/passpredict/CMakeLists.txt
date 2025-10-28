set(SOURCE_LIST passpredict.cpp)

add_executable(passpredict ${SOURCE_LIST})

set_target_properties(passpredict PROPERTIES
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

target_compile_features(passpredict PRIVATE cxx_std_17)
target_compile_options(passpredict PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-pedantic;-Wno-long-long;-Wcast-align;-Wsign-conversion>
)
target_link_libraries(passpredict PRIVATE libsgp4::sgp4)

if(WIN32 AND SGP4_BUILD_SHARED_LIBS)
  # copy sgp4 dll into the examples/passpredict directory
  add_custom_command(TARGET passpredict POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:libsgp4::sgp4>" "$<TARGET_FILE_DIR:passpredict>")
endif()
