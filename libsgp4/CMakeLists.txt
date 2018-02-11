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

  set(INCS
     CoordGeodetic.h
     CoordTopocentric.h
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
     Vector.h
     )

add_library(sgp4 STATIC ${SRCS} ${INCS})
add_library(sgp4s SHARED ${SRCS} ${INCS})
install( TARGETS sgp4s LIBRARY DESTINATION lib )
install( FILES ${INCS} DESTINATION include/SGP4 )
