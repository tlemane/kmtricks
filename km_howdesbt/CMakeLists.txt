file(GLOB CPP_FILES RELATIVE ${PROJECT_SOURCE_DIR}/src "*.cc")

add_library(howdesbt ${CPP_FILES})
add_dependencies(howdesbt SDSL CROAR LZ4)
