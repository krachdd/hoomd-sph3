#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "HOOMD::_hoomd" for configuration "Release"
set_property(TARGET HOOMD::_hoomd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HOOMD::_hoomd PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so"
  IMPORTED_SONAME_RELEASE "_hoomd.cpython-310-x86_64-linux-gnu.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS HOOMD::_hoomd )
list(APPEND _IMPORT_CHECK_FILES_FOR_HOOMD::_hoomd "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so" )

# Import target "HOOMD::quickhull" for configuration "Release"
set_property(TARGET HOOMD::quickhull APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HOOMD::quickhull PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so"
  IMPORTED_SONAME_RELEASE "libquickhull.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS HOOMD::quickhull )
list(APPEND _IMPORT_CHECK_FILES_FOR_HOOMD::quickhull "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so" )

# Import target "HOOMD::_nsearch" for configuration "Release"
set_property(TARGET HOOMD::_nsearch APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HOOMD::_nsearch PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/nsearch/_nsearch.cpython-310-x86_64-linux-gnu.so"
  IMPORTED_SONAME_RELEASE "_nsearch.cpython-310-x86_64-linux-gnu.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS HOOMD::_nsearch )
list(APPEND _IMPORT_CHECK_FILES_FOR_HOOMD::_nsearch "${_IMPORT_PREFIX}/lib/python3.10/site-packages/hoomd/nsearch/_nsearch.cpython-310-x86_64-linux-gnu.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
