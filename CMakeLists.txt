cmake_minimum_required(VERSION 3.13.0)
project(kmtricks VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(WITH_MODULES "Build all modules." OFF)
option(WITH_HOWDE "Build km_howdesbt." OFF)
option(WITH_SOCKS "Build socks interface." OFF)
option(COMPILE_TESTS "Compile tests." OFF)
option(MAKE_PACKAGE "Build package." OFF)
option(CONDA_BUILD "Build inside conda env." OFF)
option(STATIC "Static build (requires static zlib)." OFF)
option(NATIVE "Build with -march=native" ON)

set(DEV_MODE OFF)

include(CheckFilesystem)

if (STATIC AND NOT APPLE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

find_package(ZLIB REQUIRED)

if (NOT KMER_LIST)
  set(KMER_LIST "32 64 96 128")
endif()

set(KMER_LIST_GATB ${KMER_LIST})
string(REPLACE " " ";" KMER_LIST2 ${KMER_LIST})
list(LENGTH KMER_LIST2 LLEN)
set(KMER_N ${LLEN})
string(REPLACE " " "," KMER_LIST ${KMER_LIST})

if (NOT MAX_K)
  set(MAX_K 32)
endif()

if (NOT MAX_C)
  set(MAX_C 4294967295)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

set(IS_CONDA_BUILD No)
if (CONDA_BUILD)
  set(IS_CONDA_BUILD Yes)
endif()

set(PROJECT_DESCRIPTION "kmtricks - k-mer matrices and Bloom filters construction.")
set(CONTACT "teo.lemane@inria.fr")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if (NATIVE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -g -DNDEBUG")

if (DEV_MODE)
    add_compile_definitions(KM_DEV_MODE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_subdirectory(thirdparty)

set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)

include(GitUtils)
get_sha1(${PROJECT_SOURCE_DIR} GIT_SHA1)
get_sha1(${THIRDPARTY_DIR}/bcli BCLI_SHA1)
get_sha1(${THIRDPARTY_DIR}/fmt FMT_SHA1)
get_sha1(${THIRDPARTY_DIR}/kff-cpp-api KFF_SHA1)
get_sha1(${THIRDPARTY_DIR}/lz4 LZ4_SHA1)
get_sha1(${THIRDPARTY_DIR}/spdlog SPDLOG_SHA1)
get_sha1(${PROJECT_SOURCE_DIR}/gatb-core-stripped GATB_SHA1)
get_sha1(${THIRDPARTY_DIR}/sdsl-lite SDSL_SHA1)
get_sha1(${THIRDPARTY_DIR}/xxHash XXHASH_SHA1)
get_sha1(${THIRDPARTY_DIR}/googletest GTEST_SHA1)
get_sha1(${THIRDPARTY_DIR}/CRoaring CROAR_SHA1)
get_sha1(${THIRDPARTY_DIR}/indicators IND_SHA1)
get_sha1(${THIRDPARTY_DIR}/robin-hood-hashing ROBIN_SHA1)
get_sha1(${THIRDPARTY_DIR}/TurboPFor-Integer-Compression TURBOP_SHA1)
get_sha1(${THIRDPARTY_DIR}/cfrcat CFR_SHA1)

configure_file("${PROJECT_SOURCE_DIR}/include/kmtricks/config.h.in" ${PROJECT_BINARY_DIR}/include/kmtricks/config.hpp)

set(THIRDPARTY_BUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty)

set(SPD_INCLUDE ${THIRDPARTY_DIR}/spdlog/include)
set(FMT_INCLUDE ${THIRDPARTY_DIR}/fmt/include)
set(BCLI_INCLUDE ${THIRDPARTY_DIR}/bcli/include)
set(LZ4_INCLUDE ${THIRDPARTY_DIR}/lz4/lib)
set(TURBOP_INCLUDE ${THIRDPARTY_DIR}/TurboPFor-Integer-Compression)
set(GATB_INCLUDE ${THIRDPARTY_DIR}/gatb-core-stripped/src)
set(GATB_INCLUDE ${GATB_INCLUDE} ${THIRDPARTY_DIR}/gatb-core-stripped/thirdparty)
set(GATB_INCLUDE ${GATB_INCLUDE} ${THIRDPARTY_BUILD_DIR}/GATB/src/GATB-build/include)
set(KFF_INCLUDE ${THIRDPARTY_BUILD_DIR}/KFF/src/KFF-build)
set(GTEST_INCLUDE ${THIRDPARTY_DIR}/googletest/googletest/include)
set(SDSL_INCLUDE ${THIRDPARTY_BUILD_DIR}/SDSL/src/SDSL-build/include)
set(CROAR_INCLUDE ${THIRDPARTY_DIR}/CRoaring/include)
set(HOWDE_INCLUDE ${PROJECT_SOURCE_DIR}/km_howdesbt)
set(INDICATORS_INCLUDE ${THIRDPARTY_DIR}/indicators/include)
set(XXHASH_INCLUDE ${THIRDPARTY_DIR}/xxHash/)
set(ROBIN_INCLUDE ${THIRDPARTY_DIR}/robin-hood-hashing/src/include)
set(CFR_INCLUDE ${THIRDPARTY_DIR}/cfrcat/include)
set(FMT_LIB ${THIRDPARTY_BUILD_DIR}/FMT/src/FMT-build/)
set(KFF_LIB ${THIRDPARTY_BUILD_DIR}/KFF/src/KFF-build)
set(LZ4_LIB ${THIRDPARTY_BUILD_DIR}/LZ4/src/LZ4-build)
set(GTEST_LIB ${THIRDPARTY_BUILD_DIR}/GTEST/src/GTEST-build/lib)
set(SDSL_LIB ${THIRDPARTY_BUILD_DIR}/SDSL/src/SDSL-build/lib)
set(CROAR_LIB ${THIRDPARTY_BUILD_DIR}/CROAR/src/CROAR-build/src)
set(XXHASH_LIB ${THIRDPARTY_BUILD_DIR}/XXHASH/src/XXHASH-build)
set(TURBOP_LIB ${THIRDPARTY_BUILD_DIR}/TURBOP)

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  set(GATB_LIB_DIR ${THIRDPARTY_BUILD_DIR}/GATB/src/GATB-build/lib/Release)
else()
  set(GATB_LIB_DIR ${THIRDPARTY_BUILD_DIR}/GATB/src/GATB-build/lib/Debug)
endif()

set(INCLUDES ${SPD_INCLUDE} ${FMT_INCLUDE} ${BCLI_INCLUDE} ${GTEST_INCLUDE} ${KFF_INCLUDE}
             ${ZLIB_INCLUDE_DIRS} ${LZ4_INCLUDE} ${GATB_INCLUDE} ${HOWDE_INCLUDE} ${INDICATORS_INCLUDE}
             ${SDSL_INCLUDE} ${CROAR_INCLUDE} ${XXHASH_INCLUDE} ${ROBIN_INCLUDE} ${TURBOP_INCLUDE}
             ${CFR_INCLUDE})

set(KM_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)
set(DEPS SPDLOG FMT LZ4 KFF GATB SDSL XXHASH GATB TURBOP)

include_directories(SYSTEM ${INCLUDES})
include_directories(${KM_INCLUDES})

link_directories(${FMT_LIB} ${KMTRICKS_LIB} ${KFF_LIB} ${LZ4_LIB} ${GATB_LIB_DIR} ${SDSL_LIB}
                 ${CROAR_LIB} ${XXHASH_LIB} ${TURBOP_LIB})

link_libraries(fmt lz4 pthread dl kff gatbcore ZLIB::ZLIB xxhash sdsl ic)

if (WITH_HOWDE)
  add_subdirectory(km_howdesbt)
endif()

add_subdirectory(src)

if (COMPILE_TESTS)
  message(STATUS "COMPILE_TESTS=ON - Add target test.")
  enable_testing()
  add_subdirectory(tests)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release" AND MAKE_PACKAGE)
  message(STATUS "MAKE_PACKAGE=ON - Add target package.")
  include(CPackConfig)
endif()

add_custom_target(end ALL
        ${CMAKE_COMMAND} -E cmake_echo_color --cyan
        "Binaries: ${PROJECT_SOURCE_DIR}/bin"
        COMMENT "" )


if (WITH_SOCKS)
  add_dependencies(end ${PROJECT_NAME}-socks)
endif()

if (WITH_HOWDE)
  add_dependencies(end howdesbt)
endif()

if (COMPILE_TESTS)
  add_dependencies(end ${PROJECT_NAME} ${PROJECT_NAME}-tests ${PROJECT_NAME}-task-tests)
else()
  add_dependencies(end ${PROJECT_NAME})
endif()

