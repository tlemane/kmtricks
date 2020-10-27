add_library(kmtricks STATIC kmlib.cpp)
if(OpenMP_CXX_FOUND)
        target_link_libraries(kmtricks PUBLIC OpenMP::OpenMP_CXX z lz4)
else()
target_link_libraries(kmtricks PUBLIC z lz4)
endif()
file(GLOB_RECURSE LibHeaders *.hpp)
set_target_properties(kmtricks PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/lib"
        PUBLIC_HEADER "${LibHeaders}"
)
