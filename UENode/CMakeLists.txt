set(SRCS
    UserEquipment.cc)
set(INCS
    UserEquipment.h)

add_library(UE STATIC ${SRCS} ${INCS})
add_library(UEs SHARED ${SRCS} ${INCS})

target_link_libraries(UE sgp4 Satellite)