include(ExternalProject)

ExternalProject_Add(FMT
        PREFIX FMT
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/fmt
        INSTALL_COMMAND "")

ExternalProject_Add(SPDLOG
        PREFIX SPDLOG
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/spdlog
        INSTALL_COMMAND "")

ExternalProject_Add(LZ4
        PREFIX LZ4
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/lz4
        SOURCE_SUBDIR build/cmake
        CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
        INSTALL_COMMAND "")

ExternalProject_Add(KFF
        PREFIX KFF
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kff-cpp-api
        INSTALL_COMMAND "")

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
endif()

ExternalProject_Add(GATB
        PREFIX GATB
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped
        CMAKE_ARGS ${GATB_DEBUG} -DKSIZE_LIST=${KMER_LIST_GATB}
        BUILD_ALWAYS 1
        INSTALL_COMMAND "")

ExternalProject_Add(SDSL
        PREFIX SDSL
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/sdsl-lite
        INSTALL_COMMAND "")

ExternalProject_Add(TURBOP
        PREFIX TURBOP
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/TurboPFor-Integer-Compression
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${CMAKE_COMMAND} -E env
                make
        COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/libic.a" ${PROJECT_BINARY_DIR}/thirdparty/TURBOP
        INSTALL_COMMAND make clean)

if (COMPILE_TESTS)
        ExternalProject_Add(GTEST
                PREFIX GTEST
                SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
                INSTALL_COMMAND "")
endif()

if (WITH_HOWDE)
        ExternalProject_add(CROAR
                PREFIX CROAR
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/CRoaring
                CMAKE_ARGS -DROARING_BUILD_STATIC=ON
                INSTALL_COMMAND "")
endif()

ExternalProject_Add(XXHASH
        PREFIX XXHASH
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/xxHash
        SOURCE_SUBDIR "cmake_unofficial"
        CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
        INSTALL_COMMAND "")
