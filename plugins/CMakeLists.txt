file(GLOB_RECURSE PLUGIN_SOURCES "*.cpp")

add_custom_target(plugins)
foreach(PLUGIN ${PLUGIN_SOURCES})
  get_filename_component(PLUGIN_LIB ${PLUGIN} NAME_WLE)
  add_library(${PLUGIN_LIB} SHARED ${PLUGIN})
  target_link_libraries(${PLUGIN_LIB} PRIVATE build_type_flags headers)
  target_compile_definitions(${PLUGIN_LIB} PRIVATE DMAX_C=${MAX_C})
  add_dependencies(plugins ${PLUGIN_LIB})
endforeach()

