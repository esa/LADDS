cmake_minimum_required(VERSION 3.14) #at least 3.8.2 bc target_compile_features(${EXECUTABLE_NAME} PUBLIC cxx_std_17)
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

include(clang-format)

include(autopas)

add_subdirectory(src)