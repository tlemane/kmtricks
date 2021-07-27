include(CheckCXXSourceRuns)

try_compile(HAS_FS "${CMAKE_BINARY_DIR}/TMP" "${CMAKE_SOURCE_DIR}/cmake/has_filesystem.cpp"
            CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON)

if (NOT HAS_FS)
  try_compile(HAS_FS "${CMAKE_BINARY_DIR}/TMP" "${CMAKE_SOURCE_DIR}/cmake/has_filesystem.cpp"
              CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON
              LINK_LIBRARIES c++fs)
  if (HAS_FS)
    MESSAGE(STATUS "Found std::filesystem (c++fs)")
    link_libraries(c++fs)
  else()
    try_compile(HAS_FS "${CMAKE_BINARY_DIR}/TMP" "${CMAKE_SOURCE_DIR}/cmake/has_filesystem.cpp"
              CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON
              LINK_LIBRARIES stdc++fs)
    if (HAS_FS)
      MESSAGE(STATUS "Found std::filesystem (stdc++fs)")
      link_libraries(stdc++fs)
    else()
      message(FATAL_ERROR "std::filesystem not found")
    endif()
  endif()
else()
  MESSAGE(STATUS "Found std::filesystem (std=c++17)")
endif()

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0.0)
  link_libraries(stdc++fs)
endif()