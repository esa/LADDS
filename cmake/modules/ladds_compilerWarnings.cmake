# add warning flags depending on the compiler
target_compile_options(
        ${EXECUTABLE_NAME}
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:
        -Wsuggest-override
        -Wall
        -Weffc++
        -Wno-unused-variable
        -Wno-unused-function
        >
        $<$<CXX_COMPILER_ID:Clang>:
        -Wall
        -Weffc++
        -Wextra
        -Wno-unused-parameter # triggered by functions with disabled bodies or parameters that are needed for interfaces
        >
)