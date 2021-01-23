add_library(kmtricks kmlib.cpp)
target_include_directories(kmtricks PUBLIC ${LZ_INCLUDE})
file(GLOB_RECURSE LibHeaders *.hpp)
set_target_properties(kmtricks PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/lib"
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/lib"
        PUBLIC_HEADER "${LibHeaders}"
)
