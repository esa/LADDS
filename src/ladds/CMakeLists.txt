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
        debris_sim_lib
        breakupModel_lib
)

# note: in spdlog trace is the highest level
set(
        LADDS_MAX_LOGLVL
        "DEBUG"
        CACHE
        STRING "Choose the most verbose log level that shall be compiled."
)
# let ccmake and cmake-gui offer the default build type options
set_property(CACHE LADDS_MAX_LOGLVL PROPERTY STRINGS "TRACE;DEBUG;INFO;WARN;ERR;CRITICAL;OFF")

target_compile_definitions(${EXECUTABLE_NAME}
        PUBLIC
        SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_${LADDS_MAX_LOGLVL}
        DATADIR=\"${PROJECT_SOURCE_DIR}/data/\"
)
