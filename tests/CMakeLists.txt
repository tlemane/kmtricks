file(GLOB_RECURSE TEST_FILES "*_test.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
add_executable(${PROJECT_NAME}-tests ${TEST_FILES})
add_dependencies(${PROJECT_NAME}-tests ${DEPS})
target_include_directories(${PROJECT_NAME}-tests PUBLIC ${INCLUDES})
target_compile_definitions(${PROJECT_NAME}-tests PRIVATE DMAX_C=${MAX_C})
target_link_directories(${PROJECT_NAME}-tests PUBLIC ${GTEST_LIB})
target_link_libraries(${PROJECT_NAME}-tests gtest gtest_main pthread fmt)

add_executable(${PROJECT_NAME}-task-tests task_main.cpp)
add_dependencies(${PROJECT_NAME}-task-tests ${DEPS})
target_include_directories(${PROJECT_NAME}-task-tests PUBLIC ${INCLUDES})
target_compile_definitions(${PROJECT_NAME}-task-tests PRIVATE DMAX_C=${MAX_C})
target_link_directories(${PROJECT_NAME}-task-tests PUBLIC ${GTEST_LIB})
target_link_libraries(${PROJECT_NAME}-task-tests gtest gtest_main pthread fmt)

add_test(
    NAME kmdiff-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)

add_test(
    NAME kmdiff-task-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-task-tests --verbose"
)
