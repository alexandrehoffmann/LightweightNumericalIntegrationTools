#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "LNIT::LNIT" for configuration "Release"
set_property(TARGET LNIT::LNIT APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(LNIT::LNIT PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libLNIT.a"
  )

list(APPEND _cmake_import_check_targets LNIT::LNIT )
list(APPEND _cmake_import_check_files_for_LNIT::LNIT "${_IMPORT_PREFIX}/lib/libLNIT.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
