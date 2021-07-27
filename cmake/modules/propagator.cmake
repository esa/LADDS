# Propagator library
message(STATUS "Adding Propagator.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(PropagatorRepoPath https://github.com/Wombatwarrior/BA_space_debris.git)
if (GIT_SUBMODULES_SSH)
    set(PropagatorRepoPath git@github.com:Wombatwarrior/BA_space_debris.git)
endif ()

# set(BUILD_GMOCK OFF)
# set(INSTALL_GTEST OFF)

FetchContent_Declare(
        Propagatorfetch
        GIT_REPOSITORY ${PropagatorRepoPath}
        GIT_TAG origin/GenericDebris
)

# Get Propagator source and binary directories from CMake project
FetchContent_MakeAvailable(Propagatorfetch)

# Disable warnings from the library target
target_compile_options(debris_sim_lib PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval debris_sim_lib INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(debris_sim_lib SYSTEM PUBLIC "${propval}")
