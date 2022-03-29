file(GLOB CPP_FILES RELATIVE ${PROJECT_SOURCE_DIR}/src "*.cc")

add_library(howdesbt ${CPP_FILES})
target_link_libraries(howdesbt PRIVATE build_type_flags headers links)
add_dependencies(howdesbt SDSL CROAR LZ4)
