# -*- mode: cmake; -*-
#
#  Figure out the version of the used compiler
#  Variables set by this module
#  CMAKE_CXX_COMPILER_MAJOR  major version of compiler
#  CMAKE_CXX_COMPILER_MINR   minor version of compiler
#  CMAKE_CXX_COMPILER_PATCH  patch level (e.g. gcc 4.1.0)
#

# only available in Cmake 2.8.9, 
# extract version from command line if not available
if(NOT CMAKE_CXX_COMPILER_VERSION) 
  if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)

    string(STRIP ${CMAKE_CXX_COMPILER_VERSION} CMAKE_CXX_COMPILER_VERSION)
 endif( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")

#TODO Write same for clang
endif(NOT CMAKE_CXX_COMPILER_VERSION)


