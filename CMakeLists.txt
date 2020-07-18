project(kbf)

cmake_minimum_required(VERSION 3.10.0)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (debug 1)
endif()

if (debug)
    set(FLAGS -std=c++0x -O0 -g -p)
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
    set(GATB_BUILD_LIB "lib/Debug")
else()
    set(FLAGS -std=c++0x -O3)
    set(GATB_DEBUG "")
    set(GATB_BUILD_LIB "lib")
endif()

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
set (DEPS gatb FMT SDSL)

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
add_executable(km_superk_to_kmer_count src/km_superk_to_kmer_count.cpp)
add_dependencies(km_superk_to_kmer_count ${DEPS})
target_link_libraries(km_superk_to_kmer_count ${LIBS})
#
add_executable(km_merge_within_partition src/km_merge_within_partition.cpp)
add_dependencies(km_merge_within_partition ${DEPS})
target_link_libraries(km_merge_within_partition ${LIBS})

add_executable(km_output_convert src/km_output_convert.cpp)
add_dependencies(km_output_convert ${DEPS})
target_link_libraries(km_output_convert ${LIBS} sdsl)

enable_testing()

add_executable(bitmatrix_ex libs/snippets/bitmatrix_ex.cpp)
set_target_properties(bitmatrix_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)

add_executable(merger_ex libs/snippets/merger_ex.cpp)
set_target_properties(merger_ex PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/libs/snippets)
target_link_libraries(merger_ex z)

add_test(
        NAME lib_test
        COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/libs/; ./run_test.sh"
)

add_test(
        NAME kmtricks_test
        COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/kmtricks/; ./run_test.sh"
)

MESSAGE("KMER BIT REPRESENTATION = ${KMERTYPE}")
MESSAGE("COUNT BIT REPRESENTATION = ${COUNTTYPE}")
