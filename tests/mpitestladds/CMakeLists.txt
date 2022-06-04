option(LADDS_BUILD_TESTS "Build tests based on the Google Test Framework (will be downloaded if not available)" ON)

# do nothing if nothing should be done
if (NOT LADDS_BUILD_TESTS OR NOT LADDS_MPI)
    message(STATUS "LADDS_BUILD_TESTS or LADDS_MPI is disabled. Not building MPI tests.")
    return()
endif ()

# download gtest and includes the gtest modules
include(ladds_gtest)

set(MPITEST_EXECUTABLE_NAME ladds_tests_mpi)

# collect test sources
file(
        GLOB_RECURSE
        TEST_SRC
        "*.cpp"
        "*.h"
)
# append main executable sources since these are the ones we want to test
list(APPEND TEST_SRC ${EXECUTABLE_SRC})
# kick out main to avoid having to main functions
list(FILTER TEST_SRC EXCLUDE REGEX "main.cpp")

add_executable(${MPITEST_EXECUTABLE_NAME} ${TEST_SRC})

target_compile_features(${MPITEST_EXECUTABLE_NAME} PUBLIC cxx_std_17)

target_link_libraries(
        ${MPITEST_EXECUTABLE_NAME}
        gmock
        autopas
        autopasTools
        debris_sim_lib
        breakupModel_lib
        "$<$<BOOL:${LADDS_HDF5}>:hdf5::hdf5_hl;hdf5::hdf5;h5pp::h5pp>"
)

target_include_directories(
        ${MPITEST_EXECUTABLE_NAME}
        PUBLIC
        ${LADDS_SOURCE_DIR}/src
)

target_compile_definitions(${MPITEST_EXECUTABLE_NAME}
                           PUBLIC
                           SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${LADDS_MAX_LOGLVL}
                           DATADIR=\"${PROJECT_SOURCE_DIR}/data/\"
                           CFGDIR=\"${PROJECT_SOURCE_DIR}/cfg/\"
                           CFGDIRTEST=\"${PROJECT_SOURCE_DIR}/cfg/\"
                           $<$<BOOL:${LADDS_HDF5}>:LADDS_HDF5=1>
                           )

set(LADDS_TESTS_NUM_MPI_RANKS 8 CACHE STRING "Number of MPI Ranks to run LADDS MPI tests with.")
message(STATUS "XXXXXXXXXX MPIEXEC_EXECUTABLE       : ${MPIEXEC_EXECUTABLE}")
message(STATUS "XXXXXXXXXX MPIEXEC_NUMPROC_FLAG     : ${MPIEXEC_NUMPROC_FLAG}")
message(STATUS "XXXXXXXXXX LADDS_TESTS_NUM_MPI_RANKS: ${LADDS_TESTS_NUM_MPI_RANKS}")
set_property(TARGET ${MPITEST_EXECUTABLE_NAME} PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} --tag-output ${MPIEXEC_NUMPROC_FLAG} ${LADDS_TESTS_NUM_MPI_RANKS}')
#TODO: Something is not working here...
gtest_discover_tests(${MPITEST_EXECUTABLE_NAME}
                     TEST_PREFIX "MPITests/"
                     PROPERTIES
                     ENVIRONMENT "${TSAN_OPTIONS_STR}"
                     ENVIRONMENT "${LSAN_OPTIONS_STR}"
                     )
