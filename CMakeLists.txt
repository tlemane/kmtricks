cmake_minimum_required(VERSION 3.13.0)
project(kmtricks VERSION 1.4.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(WITH_MODULES "Build all modules." OFF)
option(WITH_HOWDE "Build km_howdesbt." OFF)
option(WITH_SOCKS "Build socks interface." OFF)
option(WITH_PLUGIN "Build plugins" OFF)
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

if (DEV_MODE)
    add_compile_definitions(KM_DEV_MODE)
endif()

add_compile_definitions(NONCANONICAL)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_library(build_type_flags INTERFACE)
target_compile_options(build_type_flags INTERFACE
  $<$<CONFIG:Debug>:-O0 -g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:Profile>:-O3 -ggdb3 -DNDEBUG -fno-inline>
  $<$<CONFIG:Coverage>:-O0 -g --coverage>
)
if (NATIVE)
  target_compile_options(build_type_flags INTERFACE "-march=native")
endif()

add_library(headers INTERFACE)
target_include_directories(headers INTERFACE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include
  ${CMAKE_SOURCE_DIR}/km_howdesbt
)

add_library(links INTERFACE)
target_link_libraries(links INTERFACE pthread dl)

add_library(deps INTERFACE)

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

if (WITH_HOWDE)
  add_subdirectory(km_howdesbt)
endif()

add_subdirectory(src)

if (NOT STATIC AND WITH_PLUGIN)
  add_subdirectory(plugins)
endif()

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

