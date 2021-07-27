################################################################################
# DELIVERY
################################################################################

# We get the 'package' and 'package_source' targets
include (CPack)

# We get the user name
IF (NOT CPACK_USER_NAME)
    SET (CPACK_USER_NAME  $ENV{USER})
ENDIF (NOT CPACK_USER_NAME)

# We get the date
string (TIMESTAMP CPACK_DATE "%Y-%m-%d/%H:%M:%S")

# We may have to set (if not defined) the CPACK_GFORGE_PROJECT_NAME
IF (NOT CPACK_GFORGE_PROJECT_NAME)
    SET (CPACK_GFORGE_PROJECT_NAME  ${PROJECT_NAME})
ENDIF (NOT CPACK_GFORGE_PROJECT_NAME)

# We define the name of the bin archive
SET (CPACK_URI_BIN "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CPACK_SYSTEM_NAME}.tar.gz")
SET (CPACK_URI_BIN_INFO "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CPACK_SYSTEM_NAME}.info.txt")

# We set the text holding all the information about the delivery.
SET (CPACK_INFO_BIN ${CPACK_GFORGE_PROJECT_NAME} bin ${CPACK_PACKAGE_VERSION} on ${CPACK_DATE} for ${CPACK_SYSTEM_NAME} by ${CPACK_USER_NAME})

# We define the Inria Forge place where to place a copy of the bin archive
SET (CPACK_SERVER_ADDRESS   "${CPACK_USER_NAME}@scm.gforge.inria.fr")
SET (CPACK_SERVER_DIR_BIN   "/home/groups/${CPACK_GFORGE_PROJECT_NAME}/htdocs/versions/bin/")


################################################################################
# MAIN TARGET
################################################################################

# We add a custom target for delivery binaries
add_custom_target (delivery 

    COMMAND echo "================================================================"
    COMMAND echo ""
    COMMAND echo " Starting official delivery of GATB-CORE library ${CPACK_PACKAGE_VERSION}"
    COMMAND echo ""
    COMMAND echo "================================================================"
    COMMAND echo "Checking Inria Forge repository..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery_check_repo.sh 
    COMMAND echo "Compiling library..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery_compile.sh ${SILENT_MODE}
    COMMAND echo "Creating release tag on git repository..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/git_tag_manager.sh -M ${gatb-core_VERSION_MAJOR} -m ${gatb-core_VERSION_MINOR} -p ${gatb-core_VERSION_PATCH} -t \"'new release: ${CPACK_INFO_BIN}'\" 
    COMMAND echo "Uploading binary on Inria Forge..."
    COMMAND scp -q ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN} ${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery_dump_system.sh ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO} ${CMAKE_VERSION} ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} ${CMAKE_CXX_FLAGS} ${LIBRARY_COMPILE_DEFINITIONS}
    COMMAND scp -q ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO} ${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}
    COMMAND echo "Creating release on github..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/github_release_manager.sh -l ${GH_LOGIN} -t ${GH_TOKEN} -o ${GH_OWNER} -r ${GH_REPO} -d "v${CPACK_PACKAGE_VERSION}" -c create -m \"'new release: ${CPACK_INFO_BIN}'\" 
    COMMAND echo "Uploading binary on github..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/github_release_manager.sh -l ${GH_LOGIN} -t ${GH_TOKEN} -o ${GH_OWNER} -r ${GH_REPO} -d "v${CPACK_PACKAGE_VERSION}" -c upload ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/github_release_manager.sh -l ${GH_LOGIN} -t ${GH_TOKEN} -o ${GH_OWNER} -r ${GH_REPO} -d "v${CPACK_PACKAGE_VERSION}" -c upload ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO}
)

# We add a custom target for delivery binaries
add_custom_target (upload 

    COMMAND echo "Compiling library..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery_compile.sh true
    COMMAND echo "Uploading binary on Inria Forge..."
    COMMAND scp -q ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN} ${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery_dump_system.sh ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO} ${CMAKE_VERSION} ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} ${CMAKE_CXX_FLAGS} ${LIBRARY_COMPILE_DEFINITIONS}
    COMMAND scp -q ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO} ${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}
    COMMAND echo "Uploading binary on github..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/github_release_manager.sh -l ${GH_LOGIN} -t ${GH_TOKEN} -o ${GH_OWNER} -r ${GH_REPO} -d "v${CPACK_PACKAGE_VERSION}" -c upload ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/github_release_manager.sh -l ${GH_LOGIN} -t ${GH_TOKEN} -o ${GH_OWNER} -r ${GH_REPO} -d "v${CPACK_PACKAGE_VERSION}" -c upload ${CMAKE_CURRENT_SOURCE_DIR}/build/${CPACK_URI_BIN_INFO}
)


################################################################################
# TARGET 'help'
################################################################################

# We add a custom target for delivery sources
add_custom_target (delivery_help
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY TARGETS"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "delivery: build a targz for binaries, tag Inria repository, create a github release and upload to github"
    COMMAND echo "upload:   only build a targz for binaries and upload to github"
    COMMAND echo ""
)

