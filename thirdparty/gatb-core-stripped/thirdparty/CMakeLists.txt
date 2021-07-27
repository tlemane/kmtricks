# include other smaller libraries (json, Boophf)

add_custom_target (thirdparty_copy ALL)
add_custom_command (TARGET thirdparty_copy COMMAND ${CMAKE_COMMAND} -E
    make_directory ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json)
add_custom_command (TARGET thirdparty_copy  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/thirdparty/json/json.hpp   ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json/)

add_custom_command (TARGET thirdparty_copy COMMAND ${CMAKE_COMMAND} -E
    make_directory ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/sabuhash)
add_custom_command (TARGET thirdparty_copy  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/thirdparty/sabuhash/sabuhash.hpp   ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/sabuhash/)

#include(ExternalProject)
#
#ExternalProject_Add(LZ4
#    PREFIX LZ4
#    GIT_REPOSITORY "https://github.com/lz4/lz4"
#    GIT_TAG origin/dev
#    SOURCE_SUBDIR build/cmake
#    CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
#    INSTALL_COMMAND ""
#)

################################################################################
#  INSTALL 
################################################################################
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/json DESTINATION include)
    INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/sabuhash DESTINATION include)
ENDIF ()
