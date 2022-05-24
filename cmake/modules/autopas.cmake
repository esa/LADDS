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
        # merge after AutoPas::deleteParticle(particle &) 03.02.22
        # GIT_TAG 2b262a83a24311bdd9b50d8b4db468726a659758
        # Temporary: diffuse LB branch (bc it has lots of MPI fixes)
        #GIT_TAG 2db4f5b66980baaa495fa1d56cbcb980a50929e7
        GIT_TAG md-flex/diffuse-loadbalancing # FIXME: reset to stable hash
)
# Populate dependency
FetchContent_MakeAvailable(autopasfetch)

# Disable warnings from the library target
target_compile_options(autopas PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval autopas INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(autopas SYSTEM PUBLIC "${propval}")