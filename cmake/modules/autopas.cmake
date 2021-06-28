# autopas library
message(STATUS "Adding AutoPas.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
if (GIT_SUBMODULES_SSH)
    set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
endif ()

# Configure AutoPas to not build anything but the library
set(AUTOPAS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(AUTOPAS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(AUTOPAS_BUILD_TARGET_DOC OFF CACHE BOOL "" FORCE)
set(AUTOPAS_FORMATTING_TARGETS OFF CACHE BOOL "" FORCE)
set(AUTOPAS_OPENMP ${OPENMP} CACHE BOOL "" FORCE)
set(spdlog_ForceBundled ON CACHE BOOL "" FORCE)
set(Eigen3_ForceBundled ON CACHE BOOL "" FORCE)

FetchContent_Declare(
        autopasfetch
        GIT_REPOSITORY ${autopasRepoPath}
        GIT_TAG 7691cc8eb9315f28f0e0523373ae1fa144dbef5e
)

# Get autopas source and binary directories from CMake project
FetchContent_GetProperties(autopasfetch)

if (NOT autopasfetch_POPULATED)
    FetchContent_Populate(autopasfetch)

    add_subdirectory(${autopasfetch_SOURCE_DIR} ${autopasfetch_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
