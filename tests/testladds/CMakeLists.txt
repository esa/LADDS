option(LADDS_BUILD_TESTS "Build tests based on the Google Test Framework (will be downloaded if not available)" ON)

# do nothing if nothing should be done
if (NOT LADDS_BUILD_TESTS)
    message(STATUS "LADDS_BUILD_TESTS is disabled. Not building tests.")
    return()
endif ()

# download gtest and includes the gtest modules
include(ladds_gtest)

set(TEST_EXECUTABLE_NAME ladds_tests)

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

add_executable(${TEST_EXECUTABLE_NAME} ${TEST_SRC})

target_compile_features(${TEST_EXECUTABLE_NAME} PUBLIC cxx_std_17)

target_link_libraries(
    ${TEST_EXECUTABLE_NAME}
    gmock_main
    autopas
    autopasTools
    debris_sim_lib
    breakupModel_lib
    "$<$<BOOL:${LADDS_HDF5}>:hdf5::hdf5_hl;hdf5::hdf5;h5pp::h5pp>"
)

target_include_directories(
    ${TEST_EXECUTABLE_NAME}
    PUBLIC
    ${LADDS_SOURCE_DIR}/src
)

target_compile_definitions(${TEST_EXECUTABLE_NAME}
    PUBLIC
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${LADDS_MAX_LOGLVL}
    DATADIR=\"${PROJECT_SOURCE_DIR}/data/\"
    CFGDIR=\"${PROJECT_SOURCE_DIR}/cfg/\"
    CFGDIRTEST=\"${PROJECT_SOURCE_DIR}/cfg/\"
    $<$<BOOL:${LADDS_HDF5}>:LADDS_HDF5=1>
)

gtest_discover_tests(${TEST_EXECUTABLE_NAME}
    PROPERTIES
        ENVIRONMENT "${TSAN_OPTIONS_STR}"
        ENVIRONMENT "${LSAN_OPTIONS_STR}"
)

add_test(
        NAME validationtest.testValidation.yaml
        COMMAND
        ${PROJECT_SOURCE_DIR}/tests/validationtest.py
        ${CMAKE_BINARY_DIR}/src/ladds/ladds
        ${PROJECT_SOURCE_DIR}/cfg/testValidation.yaml
        output_209.vtu
        ${PROJECT_SOURCE_DIR}/tests/data/reference_validationTest.vtu
        CONFIGURATIONS validation
)
set_tests_properties(validationtest.testValidation.yaml
        PROPERTIES
        ENVIRONMENT "${TSAN_OPTIONS_STR};${LSAN_OPTIONS_STR}"
)
