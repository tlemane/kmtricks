include(ExternalProject)

set(THIRD_DIR ${CMAKE_SOURCE_DIR}/thirdparty)
set(THIRD_BINDIR ${CMAKE_BINARY_DIR}/thirdparty)

target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/bcli/include)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/indicators/include)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/cfrcat/include)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/robin-hood-hashing/src/include)

ExternalProject_Add(FMT
  PREFIX FMT
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/fmt
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/fmt/include)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/FMT/src/FMT-build)
target_link_libraries(links INTERFACE fmt)
add_dependencies(deps FMT)

ExternalProject_Add(SPDLOG
  PREFIX SPDLOG
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/spdlog
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/spdlog/include)
add_dependencies(deps SPDLOG)

ExternalProject_Add(LZ4
  PREFIX LZ4
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/lz4
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/lz4/lib)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/LZ4/src/LZ4-build)
target_link_libraries(links INTERFACE lz4)
add_dependencies(deps LZ4)

ExternalProject_Add(KFF
  PREFIX KFF
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/kff-cpp-api
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/KFF/src/KFF-build)
target_link_libraries(links INTERFACE kff)
add_dependencies(deps KFF)

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(GATB_DEBUG "-DCMAKE_BUILD_TYPE=Debug")
endif()

ExternalProject_Add(GATB
  PREFIX GATB
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core-stripped
  CMAKE_ARGS ${GATB_DEBUG} -DKSIZE_LIST=${KMER_LIST_GATB}
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE
  ${THIRD_DIR}/gatb-core-stripped/src
  ${THIRD_DIR}/gatb-core-stripped/thirdparty
  ${THIRD_BINDIR}/GATB/src/GATB-build/include)
target_link_directories(links INTERFACE
  $<$<CONFIG:Release>:${THIRD_BINDIR}/GATB/src/GATB-build/lib/Release>
  $<$<CONFIG:Debug>:${THIRD_BINDIR}/GATB/src/GATB-build/lib/Debug>
)
target_link_libraries(links INTERFACE gatbcore)
add_dependencies(deps GATB)

ExternalProject_Add(SDSL
  PREFIX SDSL
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/sdsl-lite
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_BINDIR}/SDSL/src/SDSL-build/include)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/SDSL/src/SDSL-build/lib)
target_link_libraries(links INTERFACE sdsl)
add_dependencies(deps SDSL)

ExternalProject_Add(TURBOP
  PREFIX TURBOP
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/TurboPFor-Integer-Compression
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E env
  make libic.a
  COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/libic.a" ${PROJECT_BINARY_DIR}/thirdparty/TURBOP
  INSTALL_COMMAND make clean
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/TurboPFor-Integer-Compression/include)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/TURBOP/)
target_link_libraries(links INTERFACE ic)
add_dependencies(deps TURBOP)

if (COMPILE_TESTS)
  ExternalProject_Add(GTEST
    PREFIX GTEST
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/googletest
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
  )
  target_include_directories(headers SYSTEM INTERFACE ${THIRD_DIR}/googletest/googletest/include)
  target_link_directories(links INTERFACE ${THIRD_BINDIR}/GTEST/src/GTEST-build/lib)
  target_link_libraries(links INTERFACE gtest gtest_main pthread)
  add_dependencies(deps GTEST)

endif()

if (WITH_HOWDE)
  ExternalProject_add(CROAR
    PREFIX CROAR
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/CRoaring
    CMAKE_ARGS -DROARING_BUILD_STATIC=ON
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
  )
  target_include_directories(headers INTERFACE ${THIRD_DIR}/CRoaring/include)
  target_link_directories(links INTERFACE ${THIRD_BINDIR}/CROAR/src/CROAR-build/src)
  add_dependencies(deps CROAR)
endif()

ExternalProject_Add(XXHASH
  PREFIX XXHASH
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/xxHash
  SOURCE_SUBDIR "cmake_unofficial"
  CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)
target_include_directories(headers INTERFACE ${THIRD_DIR}/xxHash)
target_link_directories(links INTERFACE ${THIRD_BINDIR}/XXHASH/src/XXHASH-build)
target_link_libraries(links INTERFACE xxhash)
add_dependencies(deps XXHASH)


target_link_libraries(links INTERFACE ZLIB::ZLIB)

add_library(bitpack INTERFACE)
target_include_directories(headers INTERFACE ${THIRD_DIR}/span-lite/include)
target_include_directories(headers INTERFACE ${THIRD_DIR}/bitpacker/include)
