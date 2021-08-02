# Breakup Model library
message(STATUS "Adding Breakup Model.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(BreakupModelRepoPath https://github.com/schuhmaj/nasa-breakup-model-cpp.git)
if (GIT_SUBMODULES_SSH)
    set(BreakupModelRepoPath git@github.com:schuhmaj/nasa-breakup-model-cpp.git)
endif ()

# disable breakup stuff we do not need
set(BUILD_BREAKUP_MODEL_TESTS OFF CACHE INTERNAL "")
set(BUILD_BREAKUP_MODEL_SIM OFF CACHE INTERNAL "")

FetchContent_Declare(
        BreakupModelfetch
        GIT_REPOSITORY ${BreakupModelRepoPath}
        # branch tle-fix
        GIT_TAG c308f2c162a158ecb96336d0589bb5ee3bcb450d
)

# Get Breakup Model source and binary directories from CMake project
FetchContent_MakeAvailable(BreakupModelfetch)

# Disable warnings from the library target
target_compile_options(breakupModel_lib PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval breakupModel_lib INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(breakupModel_lib SYSTEM PUBLIC "${propval}")
