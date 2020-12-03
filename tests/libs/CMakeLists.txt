file(GLOB TEST_FILES "test_*_g.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs/)
add_executable(kmtests "googletest_main.cpp" ${TEST_FILES})
target_include_directories(kmtests PUBLIC ${GINCLUDE} ${LZ_INCLUDE})
target_link_directories(kmtests PUBLIC ${GLIB_DIR} ${PROJECT_SOURCE_DIR}/bin/lib)
target_link_libraries(kmtests gtest gtest_main pthread kmtricks lz4 ZLIB::ZLIB ${OPENMP})
add_dependencies(kmtests googletest LZ4 kmtricks)

add_test(
    NAME kmtests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./kmtests"
)

add_test(
    NAME lib_test
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./run_test.sh"
)