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
        INSTALL_COMMAND "")

if (TEST)
    ExternalProject_Add(Criterion
            PREFIX Criterion
            SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/Criterion
            INSTALL_COMMAND "")
endif()