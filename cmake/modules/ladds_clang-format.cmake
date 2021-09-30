# additional target to perform clang-format run, requires clang-format
set(DESIRED_CLANG_FORMAT clang-format-9)

find_program(CLANG_FORMAT NAMES ${DESIRED_CLANG_FORMAT})

if (CLANG_FORMAT)
    message(STATUS "clang format found, added clangformat target")

    # get all project source files
    foreach (suffix IN ITEMS cpp h)
        file(
                GLOB_RECURSE
                CF_ALL_SOURCE_FILES_TMP
                "${PROJECT_SOURCE_DIR}/src/*.${suffix}"
                "${PROJECT_SOURCE_DIR}/tests/*.${suffix}"
        )
        list(APPEND CF_ALL_SOURCE_FILES ${CF_ALL_SOURCE_FILES_TMP})
    endforeach ()

    set(dummyfiles)
    foreach (_file ${CF_ALL_SOURCE_FILES})
        string(
                REPLACE
                "."
                "_"
                file_cf
                ${_file}
        )
        string(
                REPLACE
                ".."
                "."
                file_cf
                ${file_cf}
        )
        set(file_cf ".dummy/cf/${file_cf}_cf")
        add_custom_command(
                OUTPUT ${file_cf}
                COMMAND
                ${CLANG_FORMAT}
                -style=file
                -i
                ${_file}
                DEPENDS ${_file}
        )
        list(APPEND dummyfiles ${file_cf})
    endforeach ()
    add_custom_command(OUTPUT .dummy/cf/clang_dummy COMMAND true DEPENDS ${dummyfiles})
    add_custom_target(clangformat DEPENDS .dummy/cf/clang_dummy)
else ()
    message(
            STATUS
            "${DESIRED_CLANG_FORMAT} not found, not adding clang format target. Other Versions not supported!"
    )
endif ()
