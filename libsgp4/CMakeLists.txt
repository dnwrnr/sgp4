set(SRCS
    CoordGeodetic.cc
    CoordTopocentric.cc
    DateTime.cc
    DecayedException.cc
    Eci.cc
    Globals.cc
    Observer.cc
    OrbitalElements.cc
    SGP4.cc
    SatelliteException.cc
    SolarPosition.cc
    TimeSpan.cc
    Tle.cc
    TleException.cc
    Util.cc
    Vector.cc)

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
if(WIN32)
  install(TARGETS sgp4s
    RUNTIME DESTINATION lib)
else()
  install(TARGETS sgp4s
    LIBRARY DESTINATION lib)
endif()  
install( FILES ${INCS} DESTINATION include/SGP4 )
