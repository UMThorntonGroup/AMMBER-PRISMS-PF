##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25)

# Create a project for the application
project(AMMBER CXX)

# Find the PRISMS-PF package
#find_package(
#  PRISMS-PF
#  REQUIRED
#  HINTS
#    ${PRISMS_PF_DIR}
#    $ENV{PRISMS_PF_DIR}
#)

# Define the header-only library as an INTERFACE target
add_library(ammber INTERFACE)

# Add an alias so it can be used like an installed package locally
add_library(AMMBER::ammber ALIAS ammber)

# Include the headers in this repo
target_include_directories(ammber INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor>
    $<INSTALL_INTERFACE:include>
)

# Require C++20
target_compile_features(ammber INTERFACE cxx_std_20)

#target_link_libraries(
#  ammber
#  INTERFACE
#    $<IF:$<AND:$<CONFIG:Debug>,$<STREQUAL:${PRISMS_PF_BUILD_TYPE},DebugRelease>>,
#    PRISMS-PF::prisms_pf_debug,
#    PRISMS-PF::prisms_pf>
#)

# Install Instructions

# Define where the CMake config files will be installed
set(PROJECT_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Install the include directory
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY vendor/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the target to a targets file
install(TARGETS ammber
    EXPORT ${PROJECT_NAME}Targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Include standard layout paths (defines CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_LIBDIR, etc.)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

## Generate the version file (ensures find_package can check version constraints)
#write_basic_package_version_file(
#    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
#    VERSION ${PROJECT_VERSION}
#    COMPATIBILITY SameMajorVersion
#)

# Create the main Package Config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${PROJECT_CONFIG_INSTALL_DIR}
)

# Install all the generated CMake files
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${PROJECT_CONFIG_INSTALL_DIR}
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
#    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${PROJECT_CONFIG_INSTALL_DIR}
)
