include(ExternalProject)
ExternalProject_Add(gatb
        PREFIX gatb
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped
        CMAKE_ARGS ${GATB_DEBUG}
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

ExternalProject_Add(LZ4
        PREFIX lz4
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/lz4
        SOURCE_SUBDIR build/cmake
        CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF
        INSTALL_COMMAND "")

add_custom_target(copy_lz4
        COMMAND cp ${CMAKE_BINARY_DIR}/thirdparty/lz4/src/LZ4-build/liblz4.a ${CMAKE_SOURCE_DIR}/bin/lib/
        DEPENDS LZ4
)

if (TEST)
    ExternalProject_Add(googletest
        PREFIX googletest
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
        INSTALL_COMMAND "")
endif()