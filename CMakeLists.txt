##
#  CMake for the AMMBER interface library
##

# =========================================================
# CMake setup
# =========================================================

# CMake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(
    FATAL_ERROR
    "In-source build prohibited. \n"
    "Use a separate build directory with:\n  cmake -B build"
  )
endif()

# Project declaration, language, and version
project(AMMBER LANGUAGES CXX VERSION 1.0.1)

# Export compile command for IDE lsp support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# =========================================================
# Configuration setup
# =========================================================

set(AMMBER_VERSION ${PROJECT_VERSION})
message(STATUS "Configuring AMMBER v${AMMBER_VERSION}")

# =========================================================
# Find dependencies
# =========================================================

message(STATUS "Configuring external libraries")

# Find PRISMS-PF
find_package(PRISMS-PF 4.0 REQUIRED HINTS ${PRISMS_PF_DIR} $ENV{PRISMS_PF_DIR})

#=========================================================
# Configure targets
# =========================================================

message(STATUS "Configuring build targets")

# Define interface library and alias
add_library(ammber INTERFACE)
add_library(AMMBER::ammber ALIAS ammber)

# Glob the vendor headers
file(
  GLOB_RECURSE VENDOR_HEADERS
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/vendor/*.hpp"
)

# Glob the AMMBER headers
file(
  GLOB_RECURSE AMMBER_HEADERS
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

# Add headers to target so they get installed
target_sources(
  ammber
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}/vendor
        ${CMAKE_CURRENT_SOURCE_DIR}/include
      FILES ${VENDOR_HEADERS} ${AMMBER_HEADERS}
)

# Link PRISMS-PF
target_link_libraries(
  ammber
  INTERFACE
    $<IF:$<AND:$<CONFIG:Debug>,$<STREQUAL:${PRISMS_PF_BUILD_TYPE},DebugRelease>>,
    PRISMS-PF::prisms_pf_debug,
    PRISMS-PF::prisms_pf>
)

# =========================================================
# Installation
# =========================================================

message(STATUS "Configuring installation targets")

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Export the target to a targets file
install(
  TARGETS ammber
  EXPORT AMMBERTargets
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Export the targets file so downstream can use find_package(AMMBER)
install(
  EXPORT AMMBERTargets
  FILE AMMBERTargets.cmake
  NAMESPACE AMMBER::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMMBER
)

# Generate and install the config file for find_package support
configure_package_config_file(
  cmake/AMMBERConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/AMMBERConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMMBER
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/AMMBERConfigVersion.cmake
  VERSION ${AMMBER_VERSION}
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/AMMBERConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/AMMBERConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AMMBER
)

# Export targets so you can link to build directory
export(
  TARGETS ammber
  NAMESPACE AMMBER::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/AMMBERTargets.cmake"
)
