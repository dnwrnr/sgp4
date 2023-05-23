set(SRCS
        Topology.cc)
set(INCS
        Topology.h)

add_library(Topology STATIC ${SRCS} ${INCS})
add_library(Topologies SHARED ${SRCS} ${INCS})

target_link_libraries(Topology UE Satellite)