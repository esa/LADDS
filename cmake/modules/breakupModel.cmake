# Breakup Model library
message(STATUS "Adding Breakup Model.")

# Enable ExternalProject CMake module
include(FetchContent)

# Select https (default) or ssh path.
set(BreakupModelRepoPath https://github.com/esa/NASA-breakup-model-cpp.git)
if (GIT_SUBMODULES_SSH)
    set(BreakupModelRepoPath git@github.com:esa/NASA-breakup-model-cpp.git)
endif ()

# disable breakup stuff we do not need
set(BUILD_BREAKUP_MODEL_TESTS OFF CACHE INTERNAL "")
set(BUILD_BREAKUP_MODEL_SIM OFF CACHE INTERNAL "")

FetchContent_Declare(
        BreakupModelfetch
        GIT_REPOSITORY ${BreakupModelRepoPath}
        # initial commit after moving repo to ESA da4da134c1fb32799869b7eb15d323cd39c13da9 (9.12.21)
        GIT_TAG da4da134c1fb32799869b7eb15d323cd39c13da9
)

# Get Breakup Model source and binary directories from CMake project
FetchContent_MakeAvailable(BreakupModelfetch)

# fetch the satcat file from the breakup repository to our data directory
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/data)
configure_file(${breakupmodelfetch_SOURCE_DIR}/example-config/satcat.csv ${PROJECT_SOURCE_DIR}/data/satcat_breakupModel.csv COPYONLY)

# Disable warnings from the library target
target_compile_options(breakupModel_lib PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval breakupModel_lib INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(breakupModel_lib SYSTEM PUBLIC "${propval}")
