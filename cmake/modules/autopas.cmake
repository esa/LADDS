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
        # TODO: freeze 17.09.24
        GIT_TAG master
)
# Populate dependency
FetchContent_MakeAvailable(autopasfetch)

# Disable warnings from the library target
target_compile_options(autopas PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval autopas INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(autopas SYSTEM PUBLIC "${propval}")