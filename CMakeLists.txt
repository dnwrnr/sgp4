cmake_minimum_required(VERSION 3.16...4.0)

project(SGP4
  VERSION 1.1.0
  DESCRIPTION "Simplified perturbations models"
  HOMEPAGE_URL https://github.com/dnwrnr/sgp4
  LANGUAGES CXX
)

set(MAIN_PROJECT OFF)
# determine if SGP4 is built as a subproject (using add_subdirectory)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MAIN_PROJECT ON)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  option(BUILD_SHARED_LIBS "Build SGP4 as a shared library." ON)

  include(CTest)
endif()

option(SGP4_BUILD_SHARED_LIBS "Build sgp4 target as SHARED instead of STATIC."                                                     ${BUILD_SHARED_LIBS})
option(SGP4_BUILD_STATIC_LIBS "Build sgp4-static target along with sgp4."                                                          ${MAIN_PROJECT})
option(SGP4_INSTALL           "Install CMake targets during install step."                                                         ${MAIN_PROJECT})
option(SGP4_BUILD_TESTS       "Enable SGP4 tests."                                                                                 ${MAIN_PROJECT})
option(SGP4_BUILD_EXAMPLES    "Enable SGP4 examples."                                                                              ${MAIN_PROJECT})
option(SGP4_ENABLE_PIC        "Enable position-independent code (ignored for sgp4 target when SGP4_BUILD_SHARED_LIBS is enabled)." ON)

if(SGP4_BUILD_STATIC_LIBS AND NOT SGP4_BUILD_SHARED_LIBS)
  message(WARNING "SGP4_BUILD_STATIC_LIBS only makes sense if BUILD_SHARED_LIBS is set to ON.")
  message(WARNING "You might be building and installing two identical static libraries.")
  message(WARNING "NOTE: to override BUILD_SHARED_LIBS for SGP4 only, you may use SGP4_BUILD_SHARED_LIBS.")
endif()

if(NOT DEFINED SGP4_DEFAULT_TARGET_NAME)
  # Allow overriding the target name when using FetchContent / add_subdirectory.
  set(SGP4_DEFAULT_TARGET_NAME sgp4)
endif()
if(NOT DEFINED SGP4_STATIC_TARGET_NAME)
  # Allow overriding the target name when using FetchContent / add_subdirectory.
  set(SGP4_STATIC_TARGET_NAME sgp4-static)
endif()

set(SGP4_NAMESPACE_NAME libsgp4)

add_subdirectory(src)

if(SGP4_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(SGP4_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}"
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pkg-config.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  COMPATIBILITY AnyNewerVersion
)

if(SGP4_INSTALL)
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  )

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}"
  )

  install(TARGETS ${SGP4_DEFAULT_TARGET_NAME}
    EXPORT "${PROJECT_NAME}Targets"
    INCLUDES DESTINATION ${CMAKE_INSTALL_DATADIR}
    # .so on UNIX
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    # .dll on Windows
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    # .lib on Windows
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  if(SGP4_BUILD_STATIC_LIBS)
    install(TARGETS ${SGP4_STATIC_TARGET_NAME}
      EXPORT "${PROJECT_NAME}Targets"
      INCLUDES DESTINATION ${CMAKE_INSTALL_DATADIR}
      # .lib on Windows, .a on UNIX
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  endif()

  install(EXPORT "${PROJECT_NAME}Targets" NAMESPACE ${SGP4_NAMESPACE_NAME}:: DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}")

  install(DIRECTORY "include/" TYPE INCLUDE)
endif()
