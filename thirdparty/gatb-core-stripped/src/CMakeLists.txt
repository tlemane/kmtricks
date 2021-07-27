################################################################################
#  KMER SIZE DEFINITION 
################################################################################

# We define the default values for kmer sizes.
if (NOT KSIZE_DEFAULT_LIST) # only if the tool cmake didn't override it
    #list (APPEND KSIZE_DEFAULT_LIST  32   64)
    list (APPEND KSIZE_DEFAULT_LIST  32   64   96  128)
    #list (APPEND KSIZE_DEFAULT_LIST  32   64   96  128  160  192  224  256)
endif()

# We use a default list of kmer sizes if no one is provided by the user (via cmake -dKSIZE_LIST command)
if (NOT KSIZE_LIST)
    set (gatb-core-klist ${KSIZE_DEFAULT_LIST})
else()
    string(REPLACE " " ";" gatb-core-klist ${KSIZE_LIST})
endif()

set (gatb-core-klist ${gatb-core-klist} PARENT_SCOPE)

# We set all KSIZE_x variables (used in config.in files)
set (idx 0)
FOREACH (ksize ${gatb-core-klist})
    MATH(EXPR idx "${idx}+1")
    # We set the default value to KSIZE_i
    SET (KSIZE_${idx} ${ksize}) 
ENDFOREACH()

set (idx 0)
FOREACH (ksize ${gatb-core-klist})
    list (APPEND KSIZE_STRING_TYPE_TMP "boost::mpl::int_<${ksize}>")
ENDFOREACH()

string (REPLACE ";" "," KSIZE_STRING_COMMA "${gatb-core-klist}")
string (REPLACE ";" " " KSIZE_STRING_SPACE "${gatb-core-klist}")
string (REPLACE ";" "," KSIZE_STRING_TYPE  "${KSIZE_STRING_TYPE_TMP}")

# We dump some information
MESSAGE ("-- OPTIMIZED KMER SIZES INTERVALS ARE " ${KSIZE_STRING_SPACE} " <-- max supported kmer size without recompilation")


################################################################################
#  CONFIGURATION FILE
################################################################################

# We build the file holding system information
configure_file (
    ${PROJECT_SOURCE_DIR}/src/gatb/system/api/build_info.hpp.in
    ${PROJECT_BINARY_DIR}/include/gatb/system/api/build_info.hpp
)
# File holding configuration macros
configure_file (
    ${PROJECT_SOURCE_DIR}/src/gatb/system/api/config.hpp.in
    ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp
)

# We define the list of files that need template specializations for each possible kmer size
file (GLOB_RECURSE  TemplateFiles  ${PROJECT_SOURCE_DIR}/src/gatb/template/*.cpp.in)

FOREACH (KSIZE ${gatb-core-klist})

    FOREACH (template ${TemplateFiles})

        get_filename_component (template_out ${template} NAME)
        
        configure_file (
            ${template}
            ${PROJECT_BINARY_DIR}/src/gatb/template/${template_out}_${KSIZE}.cpp
        )

    ENDFOREACH ()

ENDFOREACH ()

################################################################################
#  LIBRARY GENERATION 
################################################################################

file (GLOB_RECURSE  ComponentFiles  *  ${PROJECT_BINARY_DIR}/src/gatb/*.cpp)

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

include_directories (${gatb-core-includes} ${PROJECT_BINARY_DIR/thirdparty/LZ4/src/LZ4/include})
link_directories(${PROJECT_BINARY_DIR}/thirdparty/LZ4/src/LZ4-build/)
add_library    (gatbcore-static  STATIC  ${ComponentFiles} )
#add_library   (gatbcore-dynamic SHARED  ${ComponentFiles} )

set_target_properties  (gatbcore-static   PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)

#set_target_properties (gatbcore-dynamic  PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)

################################################################################
#  INSTALLATION 
################################################################################

# We install the libraries
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    install (TARGETS gatbcore-static DESTINATION lib)
    install (DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp" PATTERN "*.pri" PATTERN "*.h")
    install (FILES ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp DESTINATION include/gatb/system/api/)
ENDIF ()
