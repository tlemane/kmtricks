################################################################################
# CPPUNIT
################################################################################

FIND_PATH (CPPUNIT_INCLUDE_DIR cppunit/extensions/HelperMacros.h
  $ENV{HOME}/include
  $ENV{HOME}/.linuxbrew/include
  /local/include
  /usr/include
  /usr/local/include
  /usr/local/Cellar/cppunit/1.12.1/include
  NO_DEFAULT_PATH
)

FIND_LIBRARY (CPPUNIT_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
    /usr/lib
    /usr/local/lib
)

# A little hack here... We check whether the static library is reachable too.
IF (EXISTS "${CPPUNIT_INCLUDE_DIR}/cppunit/extensions/HelperMacros.h")
    if (NOT EXISTS "${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a")
        message ("-- CppUnit: found headers (in ${CPPUNIT_INCLUDE_DIR}) but not the static library (${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a)")
        SET (CPPUNIT_NO_STATIC_LIB_FOUND 1)   
    endif()
endif()

IF (CPPUNIT_INCLUDE_DIR)
    IF (CPPUNIT_LIBRARY)
        SET (CPPUNIT_FOUND "YES")
        SET (CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY})
        IF (NOT CPPUNIT_NO_STATIC_LIB_FOUND)
            SET (CPPUNIT_LIBRARY_STATIC ${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a)
        ENDIF()
    ENDIF (CPPUNIT_LIBRARY)
ENDIF (CPPUNIT_INCLUDE_DIR)

IF (DEFINED CPPUNIT_FOUND)
    message("-- CppUnit FOUND (${CPPUNIT_INCLUDE_DIR})")
ELSE()
    message("-- CppUnit NOT FOUND")
ENDIF()

