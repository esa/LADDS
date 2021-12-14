option(LADDS_HDF5 "Add HDF5 and HighFive to support HDF5 output" ON)

if (LADDS_HDF5)
    #Select https(default) or ssh path.
    set(highFiveRepoPath https://github.com/BlueBrain/HighFive.git)
    if (GIT_SUBMODULES_SSH)
        set(highFiveRepoPath git@github.com:BlueBrain/HighFive.git)
    endif ()

    FetchContent_Declare(
            highFive
            GIT_REPOSITORY ${highFiveRepoPath}
            GIT_TAG v2.3.1
    )

    set(HIGHFIVE_BUILD_DOCS OFF)
    set(HIGHFIVE_EXAMPLES OFF)
    set(HIGHFIVE_UNIT_TESTS OFF)
    set(HIGHFIVE_USE_BOOST OFF)
    # TODO: Maybe something to look into in the future
    # set(HIGHFIVE_PARALLEL_HDF5 ON)

    FetchContent_MakeAvailable(highFive)

    # TODO might be more ugly since this is an INTERFACE library
    # Disable warnings from the library target
#    target_compile_options(HighFive PRIVATE -w)
    # Disable warnings from included headers
#    get_target_property(propval HighFive INTERFACE_INCLUDE_DIRECTORIES)
#    target_include_directories(HighFive SYSTEM INTERFACE "${propval}")
endif()