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
        # AutoPas::deleteParticle(particle &) 24.01.22
        GIT_TAG exposeDeleteParticleRef
)
# Populate dependency
FetchContent_MakeAvailable(autopasfetch)

# Disable warnings from the library target
target_compile_options(autopas PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval autopas INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(autopas SYSTEM PUBLIC "${propval}")