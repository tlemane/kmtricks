################################################################################
#  MACROS 
################################################################################
FUNCTION (LOOKUP_PATH name items prefix result)

    SET (TMP "/notfound" )

    FOREACH (path ${items})
        IF (EXISTS "${prefix}/${path}")
            SET (TMP "${prefix}/${path}" )
        ENDIF()
    ENDFOREACH(path)

    GET_FILENAME_COMPONENT (TMP "${TMP}" REALPATH) 

    IF (NOT EXISTS ${TMP})
        MESSAGE (FATAL_ERROR  "-- UNABLE TO FIND A DIRECTORY FOR ${name}...")
    ELSE()
        MESSAGE ("-- ${name} is here '${TMP}'")
    ENDIF()
    
    SET (${result} ${TMP} PARENT_SCOPE)
    
ENDFUNCTION()

################################################################################
# GATB-CORE MANAGEMENT
################################################################################

# we don't want to install all GATB-CORE items
SET (GATB_CORE_INSTALL_EXCLUDE "1")

# We look for the gatb-core directory
# We look for the gatb-core directory
# WARNING: GATB-core has to be in one of the following HARDCODED locations in your project!!
LOOKUP_PATH ("gatb-core" "gatb-core/gatb-core;../thirdparty/gatb-core;../../thirdparty/gatb-core/gatb-core;thirdparty/gatb-core;thirdparty/gatb-core/gatb-core" ${PROJECT_SOURCE_DIR} GATB_CORE_PATH)

# we depend on gatb-core; here, we define where to find all the required material
add_subdirectory (${GATB_CORE_PATH}  "${CMAKE_CURRENT_BINARY_DIR}/ext/gatb-core")

# We copy the source of gatb-core to the third party directory and the cmake directory
SET (CPACK_SOURCE_INSTALLED_DIRECTORIES 
    "${CMAKE_CURRENT_SOURCE_DIR}"   "."
    "${GATB_CORE_PATH}"             "thirdparty/gatb-core"
    "${CMAKE_MODULE_PATH}"          "cmake"
)

# For the source archive, we exclude some unwanted directories.
SET (CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES}  
    "gatb-core/doc"  "gatb-core/examples" "gatb-core/test" "gatb-core/scripts"
)

# We set the cmake helpers directory
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${gatb-core-cmake})

