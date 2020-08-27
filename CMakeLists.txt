cmake_minimum_required(VERSION 3.13.0)
project(kmtricks VERSION 0.0.0)

set(CONTACT "teo.lemane@inria.fr")

configure_file("./config.h.in" ${PROJECT_BINARY_DIR}/include/kmconfig.hpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(OPENMP OpenMP::OpenMP_CXX)
endif()

set (debug 0)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (debug 1)
endif()

# flags
set(CMAKE_CXX_FLAGS "-std=c++0x")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
    set(GATB_BUILD_LIB "lib/Debug")
else()
    set(GATB_DEBUG "")
    set(GATB_BUILD_LIB "lib")
endif()

add_subdirectory(thirdparty) # Build thirdparty libraries
add_subdirectory(libs/kmtricks) # Build kmtricks library

# type size
if (NOT KMERTYPE)
    set(KMERTYPE 64)
endif()

if (NOT COUNTTYPE)
    set(COUNTTYPE 8)
endif()

add_definitions(-DKTYPE=${KMERTYPE})
add_definitions(-DCNTYPE=${COUNTTYPE})

# thirdparty source directories
set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/thirdparty)
set(GATB_DIR ${THIRDPARTY_DIR}/gatb-core-stripped)
set(FMT_DIR ${THIRDPARTY_DIR}/fmt)

# thirdparty build directories
set(THIRDPARTY_BUILD_DIR ${PROJECT_BINARY_DIR}/thirdparty)
set(GATB_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/gatb/src/gatb-build)
set(SDSL_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/sdsl/src/SDSL-build)
set(FMT_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/FMT/src/FMT-build)

# includes
set(GATB_INCLUDE ${GATB_DIR}/src ${GATB_DIR}/thirdparty ${GATB_BUILD_DIR}/include)
set(FMT_INCLUDE ${FMT_DIR}/include)
set(SDSL_INCLUDE ${SDSL_BUILD_DIR}/include)
set(KMTRICKS_INCLUDE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/libs ${THIRDPARTY_DIR})
set(CONFIG_INCLUDE ${PROJECT_BINARY_DIR}/include)

include_directories(${GATB_INCLUDE} ${FMT_INCLUDE} ${SDSL_INCLUDE} ${KMTRICKS_INCLUDE} ${CONFIG_INCLUDE})

# thirdparty lib directories
set(GATB_LIB_DIR ${GATB_BUILD_DIR}/${GATB_BUILD_LIB})
set(SDSL_LIB_DIR ${SDSL_BUILD_DIR}/lib)
set(FMT_LIB_DIR ${FMT_BUILD_DIR})

# links
link_directories(${GATB_LIB_DIR} ${SDSL_LIB_DIR} ${FMT_LIB_DIR})

# build kmtricks
add_subdirectory(src)

# build snippets
add_subdirectory(libs/snippets)

# tests
enable_testing()
set(TEST_VERBOSITY 0)
if (TEST)
    if (TEST_VERBOSITY EQUAL 1)
        set(CRITERION_VERBOSE "--verbose")
    else()
        set(CRITERION_VERBOSE "")
    endif()

    set(CRITERION_DIR ${THIRDPARTY_DIR}/Criterion)
    set(CRITERION_BUILD_DIR ${THIRDPARTY_BUILD_DIR}/Criterion/src/Criterion-build)

    set(CRITERION_INCLUDE ${CRITERION_DIR}/include)
    set(CRITERION_LIB_DIR ${CRITERION_BUILD_DIR})

    add_subdirectory(tests/libs) # kmtricks lib tests
    add_subdirectory(tests/kmtricks) # kmtricks pipeline test
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
            INSTALL_COMMAND "")
    add_dependencies(howdesbt SDSL)
endif()

set(INFOS ${PROJECT_BINARY_DIR}/build_infos.txt)
file(REMOVE ${INFOS})
file(APPEND ${INFOS} "kmtricks build summary\n\n")
file(APPEND ${INFOS} "KMERTYPE=${KMERTYPE}\n")
file(APPEND ${INFOS} "COUNTTYPE=${COUNTTYPE}\n")
file(APPEND ${INFOS} "Debug=${debug}\n")
file(APPEND ${INFOS} "TEST=${TEST}, TEST_VERBOSITY=${TEST_VERBOSITY}\n")
file(APPEND ${INFOS} "Binaries : ${PROJECT_SOURCE_DIR}/bin\n")
file(APPEND ${INFOS} "kmtricks static library : ${PROJECT_SOURCE_DIR}/bin/lib\n")
file(APPEND ${INFOS} "kmtricks ex snippets : ${PROJECT_SOURCE_DIR}/bin/snippets\n")
if (HOWDE)
    file(APPEND ${INFOS} "howdesbt binary : ${PROJECT_SOURCE_DIR}/bin/howdesbt\n")
endif()
file(APPEND ${INFOS}  "Run test : ctest CTestTestfile.cmake --verbose\n")
