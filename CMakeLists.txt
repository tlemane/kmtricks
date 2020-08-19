project(kbf)

cmake_minimum_required(VERSION 3.10.0)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (debug 1)
endif()

set(CMAKE_CXX_FLAGS "-std=c++0x")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

if (debug)
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
    set(GATB_BUILD_LIB "lib/Debug")
else()
    set(GATB_DEBUG "")
    set(GATB_BUILD_LIB "lib")
endif()

find_package(OpenMP)

include(ExternalProject)
ExternalProject_Add(gatb
        PREFIX gatb
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped
        CMAKE_ARGS
            ${GATB_DEBUG}
        BUILD_ALWAYS 1
        INSTALL_COMMAND "")

ExternalProject_Add(FMT
        PREFIX FMT
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/fmt
        INSTALL_COMMAND "")

ExternalProject_Add(SDSL
        PREFIX sdsl
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/sdsl-lite
        INSTALL_COMMAND "")

if (TEST)
    ExternalProject_Add(Criterion
            PREFIX Criterion
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/Criterion
            INSTALL_COMMAND "")
endif()

add_definitions(${FLAGS})

if (NOT KMERTYPE)
    set(KMERTYPE 64)
endif()


if (NOT COUNTTYPE)
    set(COUNTTYPE 8)
endif()

add_definitions(-DKTYPE=${KMERTYPE})
add_definitions(-DCNTYPE=${COUNTTYPE})

include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped/src)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/fmt/include)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped/thirdparty)
include_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/include)
include_directories(${PROJECT_BINARY_DIR}/sdsl/src/SDSL-build/include)
if (debug)
    include_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/include/Debug)
endif()
include_directories(${PROJECT_SOURCE_DIR}/thirdparty)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/libs)

file(GLOB_RECURSE ProjectFiles src/*.cpp)
link_directories(${PROJECT_BINARY_DIR}/FMT/src/FMT-build)
link_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/${GATB_BUILD_LIB})
link_directories(${PROJECT_BINARY_DIR}/sdsl/src/SDSL-build/lib)
#
install(DIRECTORY DESTINATION ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set (LIBS gatbcore fmt z pthread)
set (DEPS gatb FMT)

add_executable(kmtricks src/kmtricks.cpp)
add_dependencies(kmtricks ${DEPS})
target_link_libraries(kmtricks ${LIBS})
#
add_executable(km_minim_repart src/km_minim_repart.cpp)
add_dependencies(km_minim_repart ${DEPS})
target_link_libraries(km_minim_repart ${LIBS})
#
add_executable(km_reads_to_superk src/km_reads_to_superk.cpp)
add_dependencies(km_reads_to_superk ${DEPS})
target_link_libraries(km_reads_to_superk ${LIBS})
#
add_executable(km_superk_to_kmer_counts src/km_superk_to_kmer_counts.cpp)
add_dependencies(km_superk_to_kmer_counts ${DEPS})
target_link_libraries(km_superk_to_kmer_counts ${LIBS})
#
add_executable(km_merge_within_partition src/km_merge_within_partition.cpp)
add_dependencies(km_merge_within_partition ${DEPS})
if(OpenMP_CXX_FOUND)
    target_link_libraries(km_merge_within_partition PUBLIC OpenMP::OpenMP_CXX ${LIBS})
else()
    target_link_libraries(km_merge_within_partition ${LIBS})
endif()

add_executable(km_output_convert src/km_output_convert.cpp)
add_dependencies(km_output_convert ${DEPS} SDSL)
target_link_libraries(km_output_convert ${LIBS} sdsl)

enable_testing()

add_executable(bitmatrix_ex libs/snippets/bitmatrix_ex.cpp)
set_target_properties(bitmatrix_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
if(OpenMP_CXX_FOUND)
    target_link_libraries(bitmatrix_ex PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(bitmatrix_ex PRIVATE -Wall)


add_executable(merger_ex libs/snippets/merger_ex.cpp)
set_target_properties(merger_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_link_libraries(merger_ex z)
target_compile_options(merger_ex PRIVATE -Wall)

add_executable(minim_ex libs/snippets/minim_repart_ex.cpp)
set_target_properties(minim_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_compile_options(minim_ex PRIVATE -Wall)

add_executable(sequences_ex libs/snippets/sequences_ex.cpp)
set_target_properties(sequences_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_compile_options(sequences_ex PRIVATE -Wall)

add_executable(skreader_ex libs/snippets/skreader_ex.cpp)
set_target_properties(skreader_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_compile_options(skreader_ex PRIVATE -Wall)

add_executable(repartition_ex libs/snippets/repartition_ex.cpp)
set_target_properties(repartition_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_compile_options(repartition_ex PRIVATE -Wall)


if (TEST)
    set(TEST_DEPS Criterion)
    set(TEST_LIBS ${PROJECT_BINARY_DIR}/Criterion/src/Criterion-build)
    set(TEST_INCS ${PROJECT_SOURCE_DIR}/thirdparty/Criterion/include ${PROJECT_SOURCE_DIR}/libs/kmtricks)

    link_directories(${TEST_LIBS})
    include_directories(${TEST_INCS})
    link_libraries(criterion)

    add_executable(test_code tests/libs/test_code.cpp)
    set_target_properties(test_code PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    add_dependencies(test_code ${TEST_DEPS})

    add_executable(test_sequences tests/libs/test_sequences.cpp)
    set_target_properties(test_sequences PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    add_dependencies(test_sequences ${TEST_DEPS})

    add_executable(test_skreader tests/libs/test_skreader.cpp)
    set_target_properties(test_skreader PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    add_dependencies(test_skreader ${TEST_DEPS})

    add_executable(test_merger tests/libs/test_merger.cpp)
    set_target_properties(test_merger PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    target_link_libraries(test_merger z)
    add_dependencies(test_merger ${TEST_DEPS})

    add_executable(test_bitmatrix tests/libs/test_bitmatrix.cpp)
    set_target_properties(test_bitmatrix PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    add_dependencies(test_bitmatrix ${TEST_DEPS})

    add_executable(test_repartition tests/libs/test_repartition.cpp)
    set_target_properties(test_repartition PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/libs)
    add_dependencies(test_repartition ${TEST_DEPS})

    if (TEST_VERBOSITY EQUAL 1)
        set(CRITERION_VERBOSE "--verbose")
    else()
        set(CRITERION_VERBOSE "")
    endif()

    add_test(
            NAME test_code
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_code ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME test_sequences
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_sequences ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME test_skreader
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_skreader ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME test_merger
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_merger ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME test_bitmatrix
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_bitmatrix ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME test_repartition
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/ ; ./test_repartition ${CRITERION_VERBOSE}"
    )

    add_test(
            NAME lib_test
            COMMAND sh -c  "cd ${PROJECT_SOURCE_DIR}/tests/libs/; ./run_test.sh"
    )

    add_test(
            NAME kmtricks_test
            COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/kmtricks/; ./run_test.sh"
    )


endif()


MESSAGE("KMER BIT REPRESENTATION = ${KMERTYPE}")
MESSAGE("COUNT BIT REPRESENTATION = ${COUNTTYPE}")
