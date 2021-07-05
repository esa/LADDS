# autopas library
message(STATUS "Adding AutoPas.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(autopasRepoPath https://github.com/AutoPas/AutoPas.git)
if (GIT_SUBMODULES_SSH)
    set(autopasRepoPath git@github.com:AutoPas/AutoPas.git)
endif ()

FetchContent_Declare(
        autopasfetch
        GIT_REPOSITORY ${autopasRepoPath}
        GIT_TAG 427939164532052d8155cd0c5bff7dfc7d4cb77c
)

# Get autopas source and binary directories from CMake project
FetchContent_GetProperties(autopasfetch)

if (NOT autopasfetch_POPULATED)
    FetchContent_Populate(autopasfetch)

    add_subdirectory(${autopasfetch_SOURCE_DIR} ${autopasfetch_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
