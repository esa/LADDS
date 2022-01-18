# Propagator library
message(STATUS "Adding Propagator.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(PropagatorRepoPath https://github.com/FG-TUM/OrbitPropagator.git)
if (GIT_SUBMODULES_SSH)
    set(PropagatorRepoPath git@github.com:FG-TUM/OrbitPropagator.git)
endif ()

# disable stuff we do not need
set(DebrisSim_Thesis OFF CACHE INTERNAL "")
set(DebrisSim_Tests OFF CACHE INTERNAL "")
set(DebrisSim_Simulation OFF CACHE INTERNAL "")

FetchContent_Declare(
        Propagatorfetch
        GIT_REPOSITORY ${PropagatorRepoPath}
        # branch: main; latest feature: ActivityState
        GIT_TAG bc1504c86d0025378caf9d1b34508f30d64b5404
)

# Get Propagator source and binary directories from CMake project
FetchContent_MakeAvailable(Propagatorfetch)

# Disable warnings from the library target
target_compile_options(debris_sim_lib PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval debris_sim_lib INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(debris_sim_lib SYSTEM PUBLIC "${propval}")
