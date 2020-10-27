file(GLOB TestFiles RELATIVE ${PROJECT_SOURCE_DIR}/tests/libs *.cpp)
foreach( source_file ${TestFiles})
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    string ( REPLACE ".cpp" "" name ${source_file})
    add_executable(${name} ${source_file})
    add_dependencies(${name} Criterion LZ4)
    target_include_directories(${name} PUBLIC ${CRITERION_INCLUDE} ${LZ_INCLUDE})
    target_link_directories(${name} PUBLIC ${CRITERION_LIB_DIR} ${TESTTEST})
    target_link_libraries(${name} criterion z ${OPENMP} ${TESTTESTL} lz4)
    add_test(
            NAME ${name}
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./${name} ${CRITERION_VERBOSE}"
    )
endforeach( source_file ${TestFiles})

add_test(
        NAME lib_test
        COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./run_test.sh"
)