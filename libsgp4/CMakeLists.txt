set(SGP4_SOURCES
    CoordGeodetic.cpp
    CoordTopocentric.cpp
    CsvTleLoader.cpp
    DateTime.cpp
    DecayedException.cpp
    Eci.cpp
    Globals.cpp
    Observer.cpp
    OrbitalElements.cpp
    SGP4.cpp
    SatelliteException.cpp
    SolarPosition.cpp
    TimeSpan.cpp
    Tle.cpp
    TleException.cpp
    Util.cpp
    Vector.cpp)

set(SGP4_HEADERS
    CoordGeodetic.h
    CoordTopocentric.h
    CsvTleLoader.h
    DateTime.h
    DecayedException.h
    Eci.h
    Globals.h
    Observer.h
    OrbitalElements.h
    SatelliteException.h
    SGP4.h
    SolarPosition.h
    TimeSpan.h
    TleException.h
    Tle.h
    Util.h
    Vector.h)

add_library(sgp4 STATIC ${SGP4_SOURCES})
add_library(sgp4s SHARED ${SGP4_SOURCES})

add_library(SGP4::sgp4 ALIAS sgp4)
add_library(SGP4::sgp4s ALIAS sgp4s)

target_include_directories(sgp4
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_include_directories(sgp4s
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set_target_properties(sgp4s PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    WINDOWS_EXPORT_ALL_SYMBOLS ON)

foreach(target sgp4 sgp4s)
    target_compile_features(${target} PUBLIC cxx_std_17)

    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wshadow -Wnon-virtual-dtor
            -pedantic -Wno-long-long -Wcast-align -Wsign-conversion
            -fno-math-errno)
    endif()
endforeach()

install(TARGETS sgp4 sgp4s
    EXPORT SGP4Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${SGP4_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libsgp4)

install(EXPORT SGP4Targets
    FILE SGP4Targets.cmake
    NAMESPACE SGP4::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SGP4)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/SGP4Config.cmake.in"
    "${PROJECT_BINARY_DIR}/SGP4Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SGP4)

write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/SGP4ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
    "${PROJECT_BINARY_DIR}/SGP4Config.cmake"
    "${PROJECT_BINARY_DIR}/SGP4ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SGP4)
