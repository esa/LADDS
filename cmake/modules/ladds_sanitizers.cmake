option(LADDS_ENABLE_ADDRESS_SANITIZER "Adds clang's address sanitizer to CMAKE_CXX_FLAGS and CMAKE_LINKER_FLAGS" OFF)
option(LADDS_ENABLE_MEMORY_SANITIZER "Adds clang's memory sanitizer to CMAKE_CXX_FLAGS and CMAKE_LINKER_FLAGS" OFF)
option(LADDS_ENABLE_THREAD_SANITIZER "Adds clang's thread sanitizer to CMAKE_CXX_FLAGS and CMAKE_LINKER_FLAGS" OFF)

if (
        LADDS_ENABLE_ADDRESS_SANITIZER
        OR LADDS_ENABLE_MEMORY_SANITIZER
        OR LADDS_ENABLE_THREAD_SANITIZER
        )

    if ( NOT (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
        message(WARNING "Sanitizer flags are intended for Clang/GCC")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer")
else ()
    message(STATUS "LADDS: clang sanitizers disabled")
endif ()

if (LADDS_ENABLE_ADDRESS_SANITIZER)
    message(STATUS "LADDS: ADDRESS SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
else ()
    message(STATUS "LADDS: ADDRESS SANITIZER DISABLED")
endif ()

if (LADDS_ENABLE_MEMORY_SANITIZER)
    message(STATUS "LADDS: MEMORY SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=memory")
else ()
    message(STATUS "LADDS: MEMORY SANITIZER DISABLED")
endif ()

if (LADDS_ENABLE_THREAD_SANITIZER)
    message(STATUS "LADDS: THREAD SANITIZER ENABLED!!!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=thread")
else ()
    message(STATUS "LADDS: THREAD SANITIZER DISABLED")
endif ()

