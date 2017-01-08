set(SRCS
    CoordGeodetic.cpp
    CoordTopocentric.cpp
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

add_library(sgp4 STATIC ${SRCS})
