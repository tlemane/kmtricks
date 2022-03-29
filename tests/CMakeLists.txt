file(GLOB_RECURSE TEST_FILES "*_test.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
add_executable(${PROJECT_NAME}-tests ${TEST_FILES})
target_compile_definitions(${PROJECT_NAME}-tests PRIVATE DMAX_C=${MAX_C})
target_link_libraries(${PROJECT_NAME}-tests PRIVATE build_type_flags headers links deps)

add_executable(${PROJECT_NAME}-task-tests task_main.cpp)
target_compile_definitions(${PROJECT_NAME}-task-tests PRIVATE DMAX_C=${MAX_C})
target_link_libraries(${PROJECT_NAME}-task-tests PRIVATE build_type_flags headers links deps)

add_test(
    NAME kmtricks-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-tests --verbose"
)

add_test(
    NAME kmtricks-task-tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./${PROJECT_NAME}-task-tests --verbose"
)
