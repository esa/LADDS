option(LADDS_DOXYGEN "Create \"make doc_doxygen\" target (requires Doxygen)" ON)

# Require dot, treat the other components as optional
if (LADDS_DOXYGEN)
    find_package(Doxygen
                 REQUIRED dot
                 OPTIONAL_COMPONENTS mscgen dia)
    if (DOXYGEN_FOUND)
        set(ladds_doc ladds_doc)
        # Find all files
        file(GLOB_RECURSE ${ladds_doc}
             "${PROJECT_SOURCE_DIR}/src/*.h"
             "${PROJECT_SOURCE_DIR}/src/*.cpp"
             )
        set(DOXYGEN_INPUT_DIRECTORY ${ladds_doc})
        set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
#        set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIRECTORY}/html/index.html)
        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_GENERATE_XML NO)
        set(DOXYGEN_GENERATE_LATEX NO)
        set(DOXYGEN_HAVE_DOT YES)
        doxygen_add_docs(doc_doxygen
                         ${DOXYGEN_INPUT_DIRECTORY}
                         COMMENT "Generating doxygen documentation")
    else ()
        message(WARNING "Doxygen not found. No documentation will be generated. Please install doxygen and dot.")
        set(LADDS_DOXYGEN OFF CACHE BOOL "" FORCE)
    endif ()
endif ()
