# this module runs git to determine the version and state of the code.
# it then edits version.h so that the version is available in C++

find_package(Git)
# we need git and the folder where git stores its info
if (Git_FOUND AND EXISTS "${LADDS_SOURCE_DIR}/.git" )
    # branch name
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE LADDS_VERSION_BRANCH)
    string(STRIP "${LADDS_VERSION_BRANCH}" LADDS_VERSION_BRANCH)
    # commit hash
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            OUTPUT_VARIABLE LADDS_VERSION_HASH)
    string(STRIP "${LADDS_VERSION_HASH}" LADDS_VERSION_HASH)
    # check dirty
    execute_process(
            COMMAND ${GIT_EXECUTABLE} diff --quiet --exit-code
            RESULT_VARIABLE is_dirty_result )
    if (${is_dirty_result})
        set(LADDS_VERSION_IS_DIRTY "_dirty")
    else()
        set(LADDS_VERSION_IS_DIRTY "")
    endif ()
else()
    message(WARNING "Could not find git or ${LADDS_SOURCE_DIR}/.git folder! LADDS_VERSION will be missing information.")
endif ()

configure_file(src/ladds/Version.h.in ${LADDS_BINARY_DIR}/src/ladds/Version.h @ONLY)