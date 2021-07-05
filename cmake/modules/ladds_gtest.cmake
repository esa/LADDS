message(STATUS "Adding gtest")

include(FetchContent)

# Select https (default) or ssh path.
set(gTestRepoPath https://github.com/google/googletest.git)
if (GIT_SUBMODULES_SSH)
    set(gTestRepoPath git@github.com:google/googletest.git)
endif ()

FetchContent_Declare(
        googletest
        GIT_REPOSITORY ${gTestRepoPath}
        GIT_TAG release-1.11.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# Disable actual installation of gtest
option(INSTALL_GTEST "" OFF)
# Hide variables that should not be touched by te user
mark_as_advanced(
    BUILD_GMOCK
    INSTALL_GTEST
)
FetchContent_MakeAvailable(googletest)

# also add the gtest cmake module for cmake functions
include(GoogleTest)