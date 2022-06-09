option(LADDS_HDF5 "Add HDF5 and h5pp to support HDF5 output" ON)

if (LADDS_HDF5)

    find_package(HDF5 1.8 COMPONENTS C HL REQUIRED)

    #Select https(default) or ssh path.
    set(h5ppRepoPath https://github.com/DavidAce/h5pp.git)
    if (GIT_SUBMODULES_SSH)
        set(h5ppRepoPath git@github.com:DavidAce/h5pp.git)
    endif ()

    # explicitly tell cmake that the checks that are about to happen (via check_cxx_source_compiles in h5pp) need c++17
    set(CMAKE_REQUIRED_FLAGS -std=c++17)

    FetchContent_Declare(
            h5pp
            GIT_REPOSITORY ${h5ppRepoPath}
            GIT_TAG v1.9.0
    )

    if (LADDS_MPI)
        option(H5PP_ENABLE_MPI "" ON)
    endif ()

    FetchContent_MakeAvailable(h5pp)

    # get includes from the HEADER target and set them to system includes
    get_target_property(interfaceIncludeDirs_h5pp h5pp::headers INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(h5pp SYSTEM INTERFACE "${interfaceIncludeDirs_h5pp}")

    # remove the _FORTIFY_SOURCE definition which leaks all over the project producing warnings in debug mode
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        get_target_property(interfaceCompileDefs_hdf5 hdf5::hdf5 INTERFACE_COMPILE_DEFINITIONS)
        string(REPLACE "_FORTIFY_SOURCE=2" "" interfaceCompileDefsFixed_hdf5 "${interfaceCompileDefs_hdf5}")
        set_target_properties(hdf5::hdf5 PROPERTIES INTERFACE_COMPILE_DEFINITIONS "${interfaceCompileDefsFixed_hdf5}")
    endif ()
endif ()
