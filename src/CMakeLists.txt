set(HEADER_LIST
  "${PROJECT_SOURCE_DIR}/include/SGP4/CoordGeodetic.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/CoordTopocentric.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/DateTime.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/DecayedException.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Eci.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Globals.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Observer.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/OrbitalElements.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/SGP4.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/SatelliteException.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/SolarPosition.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/TimeSpan.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Tle.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/TleException.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Util.h"
  "${PROJECT_SOURCE_DIR}/include/SGP4/Vector.h"
)
set(SOURCE_LIST
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
  Vector.cpp
)

if(SGP4_BUILD_SHARED_LIBS)
  add_library(${SGP4_DEFAULT_TARGET_NAME} SHARED ${SOURCE_LIST})
else()
  add_library(${SGP4_DEFAULT_TARGET_NAME} STATIC ${SOURCE_LIST})
  set_target_properties(${SGP4_DEFAULT_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ${SGP4_ENABLE_PIC})
endif()
set(SGP4_TARGETS ${SGP4_DEFAULT_TARGET_NAME})

add_library(${SGP4_NAMESPACE_NAME}::${SGP4_DEFAULT_TARGET_NAME} ALIAS ${SGP4_DEFAULT_TARGET_NAME})

if(SGP4_BUILD_STATIC_LIBS)
  add_library(${SGP4_STATIC_TARGET_NAME} STATIC ${SOURCE_LIST})
  list(APPEND SGP4_TARGETS ${SGP4_STATIC_TARGET_NAME})

  add_library(${SGP4_NAMESPACE_NAME}::${SGP4_STATIC_TARGET_NAME} ALIAS ${SGP4_STATIC_TARGET_NAME})
  set_target_properties(${SGP4_STATIC_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ${SGP4_ENABLE_PIC})
endif()

foreach(TARGET ${SGP4_TARGETS})
  set_target_properties(${TARGET} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    # TODO: possibly replace with __declspec(dllexport) function annotations
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )

  target_include_directories(${TARGET} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_compile_features(${TARGET} PUBLIC cxx_std_17)

  target_compile_options(${TARGET} PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wshadow;-Wnon-virtual-dtor;-pedantic;-Wno-long-long;-Wcast-align;-Wsign-conversion>
  )
endforeach()

source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST}
)
