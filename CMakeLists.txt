cmake_minimum_required(VERSION 3.13.0)
project(kmtricks VERSION 0.0.0)

set(CONTACT "teo.lemane@inria.fr")
set(INFOS ${PROJECT_BINARY_DIR}/build_infos.txt)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

execute_process(COMMAND
  git describe --always
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  OUTPUT_STRIP_TRAILING_WHITESPACE)

include(CheckTypeSize)
set(INT128_FOUND 0)
CHECK_TYPE_SIZE("__uint128_t" INT128)
if (INT128 EQUAL 16)
    set(INT128_FOUND 1)
endif()

configure_file("./config.h.in" ${PROJECT_BINARY_DIR}/include/kmconfig.hpp)

find_package(ZLIB)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(OPENMP OpenMP::OpenMP_CXX)
endif()

set (debug 0)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (debug 1)
else ()
    set (CMAKE_BUILD_TYPE "Release")
endif()

# flags
set(SSE_FLAGS "-msse -msse2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS} -std=c++0x")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
    set(GATB_BUILD_LIB "lib/Debug")
else()
    set(GATB_DEBUG "")
    set(GATB_BUILD_LIB "lib/Release")
endif()

# build thirdparty libraries
add_subdirectory(thirdparty)

# type size
if (NOT KMER_NB_BIT)
    set(KMER_NB_BIT 32 64)
endif()
if (NOT COUNT_NB_BIT)
    set(COUNT_NB_BIT 8 16 32)
endif()

if (KMER_NB_BIT STREQUAL "ALL")
    if (INT128_FOUND EQUAL 1)
        set(KMER_NB_BIT 8 16 32 64 128)
    else()
        set(KMER_NB_BIT 8 16 32 64)
    endif()
endif()

if (COUNT_NB_BIT STREQUAL "ALL")
    set (COUNT_NB_BIT 8 16 32)
endif()

message("KMER NB BIT ${KMER_NB_BIT}")
message("COUNT NB BIT ${COUNT_NB_BIT}")

# thirdparty source directories
set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)
set(GATB_DIR ${THIRDPARTY_DIR}/gatb-core-stripped)
set(FMT_DIR ${THIRDPARTY_DIR}/fmt)
set(LZ_DIR ${THIRDPARTY_DIR}/lz4)

# thirdparty build directories
set(THIRDPARTY_BUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty)
set(GATB_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/gatb/src/gatb-build)
set(SDSL_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/sdsl/src/SDSL-build)
set(FMT_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/FMT/src/FMT-build)
set(LZ_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/lz4/src/LZ4-build)

# includes
set(GATB_INCLUDE ${GATB_DIR}/src ${GATB_DIR}/thirdparty ${GATB_BUILD_DIR}/include)
set(FMT_INCLUDE ${FMT_DIR}/include)
set(SDSL_INCLUDE ${SDSL_BUILD_DIR}/include)
set(KMTRICKS_INCLUDE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/libs ${THIRDPARTY_DIR})
set(LZ_INCLUDE ${LZ_DIR}/lib)
set(CONFIG_INCLUDE ${PROJECT_BINARY_DIR}/include)

include_directories(${ZLIB_INCLUDE_DIR} ${GATB_INCLUDE} ${FMT_INCLUDE} ${SDSL_INCLUDE} ${KMTRICKS_INCLUDE} ${CONFIG_INCLUDE})

# thirdparty lib directories
set(GATB_LIB_DIR ${GATB_BUILD_DIR}/${GATB_BUILD_LIB})
set(SDSL_LIB_DIR ${SDSL_BUILD_DIR}/lib)
set(FMT_LIB_DIR ${FMT_BUILD_DIR})
set(LZ_LIB_DIR ${LZ_BUILD_DIR})

# links
link_directories(${GATB_LIB_DIR} ${SDSL_LIB_DIR} ${FMT_LIB_DIR} ${LZ_LIB_DIR})

# build kmtricks library
add_subdirectory(libs/kmtricks) 

# build kmtricks
add_subdirectory(src)

# build snippets
add_subdirectory(libs/snippets)

# tests
enable_testing()
set(TEST_VERBOSITY 0)
if (TEST)
    set(GDIR ${THIRDPARTY_DIR}/googletest)
    set(GBUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty/googletest/src/googletest-build)
    set(GINCLUDE ${GDIR}/googletest/include)
    set(GLIB_DIR ${GBUILD_DIR}/lib)
    add_subdirectory(tests/libs) # kmtricks lib tests
    add_subdirectory(tests/kmtricks) # kmtricks pipeline test
else()
    set(TEST 0)
endif()

if (HOWDE)
    include(ExternalProject)
    ExternalProject_Add(howdesbt
        PREFIX howdesbt
        SOURCE_DIR ${THIRDPARTY_DIR}/HowDeSBT
        CMAKE_ARGS -DSDSLINC=${SDSL_INCLUDE}
        CMAKE_ARGS -DSDSLLIB=${SDSL_LIB_DIR}
        CMAKE_ARGS -DKMTRICKS_LIB=${PROJECT_SOURCE_DIR}/bin/lib/
        CMAKE_ARGS -DKMTRICKS_INC=${PROJECT_SOURCE_DIR}/libs/
        CMAKE_ARGS -DOUTPUT=${PROJECT_SOURCE_DIR}/bin/
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        INSTALL_COMMAND "")
    add_dependencies(howdesbt SDSL)
else()
    set(HOWDE 0)
endif()

file(REMOVE ${INFOS})
file(APPEND ${INFOS} "System: ${CMAKE_SYSTEM_NAME}\n")
file(APPEND ${INFOS} "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\n")
file(APPEND ${INFOS} "OpenMp: ${OpenMP_CXX_FOUND}\n")
file(APPEND ${INFOS} "OpenMpVersion: ${OpenMP_CXX_VERSION}\n")
file(APPEND ${INFOS} "HowDe: ${HOWDE}\n")
file(APPEND ${INFOS} "KMER_NB_BIT: ${KMER_NB_BIT}\n")
file(APPEND ${INFOS} "COUNT_NB_BIT: ${COUNT_NB_BIT}\n")
file(APPEND ${INFOS} "Debug: ${debug}\n")
file(APPEND ${INFOS} "TEST: ${TEST}\n")
file(APPEND ${INFOS} "TEST_VERBOSITY: ${TEST_VERBOSITY}\n")
file(APPEND ${INFOS} "GIT SHA1: ${GIT_SHA1}\n")
file(APPEND ${INFOS} "Binaries: ${PROJECT_SOURCE_DIR}/bin\n")
file(APPEND ${INFOS} "kmtricks static library: ${PROJECT_SOURCE_DIR}/bin/lib\n")
file(APPEND ${INFOS} "kmtricks ex snippets: ${PROJECT_SOURCE_DIR}/bin/snippets\n")
if (HOWDE)
    file(APPEND ${INFOS} "howdesbt binary: ${PROJECT_SOURCE_DIR}/bin/howdesbt\n")
endif()
file(APPEND ${INFOS}  "Run test: ctest CTestTestfile.cmake --verbose\n")

message("${TARGETS}")

function(get_all_targets var)
    set(targets)
    get_all_targets_recursive(targets ${CMAKE_CURRENT_SOURCE_DIR})
    set(${var} ${targets} PARENT_SCOPE)
endfunction()

macro(get_all_targets_recursive targets dir)
    get_property(subdirectories DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    foreach(subdir ${subdirectories})
        get_all_targets_recursive(${targets} ${subdir})
    endforeach()

    get_property(current_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    list(APPEND ${targets} ${current_targets})
endmacro()

get_all_targets(all_targets)
add_custom_target( end ALL
        ${CMAKE_COMMAND} -E cmake_echo_color --cyan
        "Binaries: ${PROJECT_SOURCE_DIR}/bin"
        "kmtricks static library: ${PROJECT_SOURCE_DIR}/bin/lib"
        "kmtricks ex snippets: ${PROJECT_SOURCE_DIR}/bin/snippets"
        COMMENT "" )
add_dependencies( end ${all_targets} )


if (CMAKE_BUILD_TYPE STREQUAL "Release" AND MAKE_RELEASE EQUAL 1)
    include(CPackConfig)
endif()
