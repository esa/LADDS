set(TSAN_OPTIONS_STR "TSAN_OPTIONS=suppressions=${LADDS_SOURCE_DIR}/.tsanIgnoreList.txt ignore_noninstrumented_modules=1")
set(LSAN_OPTIONS_STR "LSAN_OPTIONS=suppressions=${LADDS_SOURCE_DIR}/.lsanIgnoreList.txt ignore_noninstrumented_modules=1")

add_subdirectory(testladds)
add_subdirectory(mpitestladds)

