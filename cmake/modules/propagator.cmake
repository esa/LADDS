# Propagator library
message(STATUS "Adding Propagator.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(PropagatorRepoPath https://github.com/Wombatwarrior/BA_space_debris.git)
if (GIT_SUBMODULES_SSH)
    set(PropagatorRepoPath git@github.com:Wombatwarrior/BA_space_debris.git)
endif ()

# disable stuff we do not need
set(DebrisSim_Thesis OFF CACHE INTERNAL "")
set(DebrisSim_Tests OFF CACHE INTERNAL "")
set(DebrisSim_Simulation OFF CACHE INTERNAL "")

FetchContent_Declare(
        Propagatorfetch
        GIT_REPOSITORY ${PropagatorRepoPath}
        # branch tmp (=GenericDebris + main)
        GIT_TAG 9ad61f25686d34aa4be6642998a157c0f7cb3a5c
)

# Get Propagator source and binary directories from CMake project
FetchContent_MakeAvailable(Propagatorfetch)

# Disable warnings from the library target
target_compile_options(debris_sim_lib PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval debris_sim_lib INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(debris_sim_lib SYSTEM PUBLIC "${propval}")
