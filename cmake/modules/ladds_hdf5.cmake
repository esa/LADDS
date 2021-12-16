option(LADDS_HDF5 "Add HDF5 and h5pp to support HDF5 output" ON)

if (LADDS_HDF5)

    find_package(HDF5 1.8 COMPONENTS C HL REQUIRED)

    #Select https(default) or ssh path.
    set(h5ppRepoPath https://github.com/DavidAce/h5pp.git)
    if (GIT_SUBMODULES_SSH)
        set(h5ppRepoPath git@github.com:DavidAce/h5pp.git)
    endif ()

    FetchContent_Declare(
            h5pp
            GIT_REPOSITORY ${h5ppRepoPath}
            GIT_TAG v1.9.0
    )


    FetchContent_MakeAvailable(h5pp)

    # TODO might be more ugly since this is an INTERFACE library
    # Disable warnings from the library target
#    target_compile_options(HighFive PRIVATE -w)
    # Disable warnings from included headers
#    get_target_property(propval HighFive INTERFACE_INCLUDE_DIRECTORIES)
#    target_include_directories(HighFive SYSTEM INTERFACE "${propval}")
endif()
