cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
#at least 3.14 for FetchContent_MakeAvailable
#at least 3.19 for a correctly working FindHDF5 cmake module
project(LADDS)

set(EXECUTABLE_NAME ladds)

set(LADDS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LADDS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")
# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(
            CMAKE_BUILD_TYPE
            "Release"
            CACHE
            STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE
    )
endif ()

# Offer option to use submodules / fetchcontent with ssh instead of https.
# This might be needed on machines with restricted outward connections.
option(GIT_SUBMODULES_SSH "Use SSH for git submodules instead of HTTPS" OFF)

# add custom modules directory
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

# modules that do NOT need the main target
include(version)
include(ladds_sanitizers)
include(ladds_clang-format)
include(ladds_doxygen)
include(ladds_mpi)

include(autopas)
include(propagator)
include(breakupModel)
include(ladds_hdf5)

add_subdirectory(src)

# needs to be in the root directory
enable_testing()
add_subdirectory(tests)

# modules that need the main target
include(ladds_compilerWarnings)

# fetch the default configuration file from the cfg folder
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/cfg)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfg/default_cfg.yaml ${CMAKE_CURRENT_BINARY_DIR}/cfg/default_cfg.yaml COPYONLY)
