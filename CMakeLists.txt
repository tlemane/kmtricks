project(kbf)

cmake_minimum_required(VERSION 3.10.0)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (debug 1)
endif()

if (debug)
    set(FLAGS -O0 -g -p)
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
    set(GATB_BUILD_LIB "lib/Debug")
else()
    set(FLAGS -O3)
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

ExternalProject_Add(sdsl
        PREFIX sdsl
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/sdsl-lite
        INSTALL_COMMAND "")

add_definitions(${FLAGS})
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped/src)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/fmt/include)
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped/thirdparty)
include_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/include)
if (debug)
    include_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/include/Debug)
endif()
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/libs)

file(GLOB_RECURSE ProjectFiles src/*.cpp)
link_directories(${PROJECT_BINARY_DIR}/FMT/src/FMT-build)
link_directories(${PROJECT_BINARY_DIR}/gatb/src/gatb-build/${GATB_BUILD_LIB})
#
install(DIRECTORY DESTINATION ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set (LIBS gatbcore fmt z pthread)
set (DEPS gatb FMT sdsl)

add_executable(kmtricks src/kmtricks.cpp)
add_dependencies(kmtricks ${DEPS})
target_link_libraries(kmtricks ${LIBS})
#
add_executable(km_part src/kmpart.cpp)
add_dependencies(km_part ${DEPS})
target_link_libraries(km_part ${LIBS})
#
add_executable(km_superk src/kmsuperk.cpp)
add_dependencies(km_superk ${DEPS})
target_link_libraries(km_superk ${LIBS})
#
add_executable(km_count src/kmcount.cpp)
add_dependencies(km_count ${DEPS})
target_link_libraries(km_count ${LIBS})
#
add_executable(km_merge src/kmmerge.cpp)
add_dependencies(km_merge ${DEPS})
target_link_libraries(km_merge ${LIBS})

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
