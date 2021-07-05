file(
        GLOB_RECURSE
        EXECUTABLE_SRC
        "*.cpp"
        "*.h"
)
# Create cache entry for this var so it is accessible from everywhere
set(EXECUTABLE_SRC ${EXECUTABLE_SRC} CACHE INTERNAL "Source files for the main executable")

add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SRC})

target_compile_features(${EXECUTABLE_NAME} PUBLIC cxx_std_17)

target_include_directories(${EXECUTABLE_NAME} PUBLIC ${LADDS_SOURCE_DIR}/src/)

target_link_libraries(${EXECUTABLE_NAME}
                      autopas
                      )