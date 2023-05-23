set(SRCS
        Satellite.cc)
set(INCS
        Satellite.h)

add_library(Satellite STATIC ${SRCS} ${INCS})
add_library(Satellites SHARED ${SRCS} ${INCS})

target_link_libraries(Satellite sgp4 UE)