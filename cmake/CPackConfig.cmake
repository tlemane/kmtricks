SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "kmtricks")
SET (CPACK_PACKAGE_VENDOR               "T. Lemane, R. Chikhi, P. Peterlongo")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${CMAKE_PROJECT_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${CMAKE_PROJECT_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${CMAKE_PROJECT_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${CMAKE_PROJECT_VERSION}")  
SET (CPACK_PACKAGE_FILE_NAME            "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME}")
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")
SET (CPACK_SET_DESTDIR true)
SET (CPACK_INSTALL_PREFIX /)

INSTALL (TARGETS kmtricks
        LIBRARY
        DESTINATION lib
        PUBLIC_HEADER
        DESTINATION include/kmtricks
        )

foreach( target_name ${all_targets} )
    if (${target_name} MATCHES "-bin$")
        INSTALL (TARGETS ${target_name}
                DESTINATION bin 
                )
    elseif(${target_name} MATCHES "_ex$")
        INSTALL (TARGETS ${target_name}
                DESTINATION bin/snippets 
                )
    endif()
endforeach()

INSTALL (FILES ${CMAKE_BINARY_DIR}/thirdparty/lz4/src/LZ4-build/liblz4.a DESTINATION ext/lib)
INSTALL (FILES ${CMAKE_SOURCE_DIR}/kmtricks.py DESTINATION .)
INSTALL (FILES ${CMAKE_SOURCE_DIR}/README.md DESTINATION .)
INSTALL (FILES ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION .)

include (CPack)
