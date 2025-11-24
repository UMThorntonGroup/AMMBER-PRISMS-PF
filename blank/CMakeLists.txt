##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.13.4)

# Create a project for the application
project(myapp CXX)

include_directories(${CMAKE_SOURCE_DIR}/../include)
include_directories(${CMAKE_SOURCE_DIR}/../vendor)
# Find the PRISMS-PF package
find_package(prisms_pf QUIET HINTS ${PRISMS_PF_DIR} $ENV{PRISMS_PF_DIR})
if(PRISMS_PF_FOUND)
    # Include setup script
    include(${PRISMS_PF_CORE_DIR}/cmake/setup_application.cmake)

    # Set up the application
    prisms_pf_autopilot(${PRISMS_PF_CORE_DIR})
else()
    message(
        FATAL_ERROR
        "PRISMS-PF package not found. Please make sure it is in path."
    )
endif()

#add_executable(generate_bc_amr
#  ${CMAKE_SOURCE_DIR}/generate_bc_amr.cc)
#  set_property(TARGET generate_bc_amr PROPERTY OUTPUT_NAME generate-bc-amr)
#  deal_ii_setup_target(generate_bc_amr RELEASE)