add_executable(${PROJECT_NAME} kmtricks.cpp cli.cpp utils.cpp)
target_compile_definitions(${PROJECT_NAME} PRIVATE DMAX_C=${MAX_C})
target_link_libraries(${PROJECT_NAME} PRIVATE build_type_flags headers links deps)

if (WITH_PLUGIN)
  target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_PLUGIN)
endif()

if (WITH_HOWDE)
  target_link_libraries(${PROJECT_NAME} PRIVATE howdesbt roaring)
  target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_HOWDE)
endif()

if (WITH_MODULES)
  target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_KM_MODULES)
endif()

if (NOT APPLE)
  target_link_options(${PROJECT_NAME} PUBLIC -export-dynamic)
  if (STATIC)
    target_link_options(${PROJECT_NAME} PUBLIC -static)
  endif()
endif()

target_link_libraries(${PROJECT_NAME} PRIVATE pthread)

if (WITH_SOCKS AND WITH_HOWDE)
  add_executable(${PROJECT_NAME}-socks kmtricks-socks.cpp utils.cpp socks_utils.cpp)
  target_compile_definitions(${PROJECT_NAME}-socks PRIVATE DMAX_K=${MAX_K} DMAX_C=${MAX_C})
  if (STATIC AND NOT APPLE)
    target_link_options(${PROJECT_NAME}-socks PUBLIC -static)
  endif()
  target_link_libraries(${PROJECT_NAME}-socks PRIVATE build_type_flags headers howdesbt roaring pthread links deps)
  target_compile_definitions(${PROJECT_NAME}-socks PRIVATE WITH_HOWDE)
endif()


