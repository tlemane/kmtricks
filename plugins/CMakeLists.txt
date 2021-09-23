file(GLOB_RECURSE PLUGIN_SOURCES "*.cpp")

foreach(PLUGIN ${PLUGIN_SOURCES})
  get_filename_component(PLUGIN_LIB ${PLUGIN} NAME_WLE)
  add_library(${PLUGIN_LIB} SHARED ${PLUGIN})
  target_include_directories(${PLUGIN_LIB} PRIVATE ${INCLUDES})
  target_compile_definitions(${PLUGIN_LIB} PRIVATE DMAX_C=${MAX_C})
endforeach()

