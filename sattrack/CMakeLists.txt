add_executable(sattrack sattrack.cpp)
target_link_libraries(sattrack PRIVATE SGP4::sgp4)
