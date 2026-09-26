add_executable(passpredict passpredict.cpp)
target_link_libraries(passpredict PRIVATE SGP4::sgp4)
