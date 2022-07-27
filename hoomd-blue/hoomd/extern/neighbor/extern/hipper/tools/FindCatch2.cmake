if(Catch2_FOUND)
    return()
endif()

find_path(Catch2_INCLUDE_DIR
          NAMES catch2/catch.hpp
          PATHS ${Catch2_ROOT} ENV Catch2_ROOT
          PATH_SUFFIXES single_include
          NO_DEFAULT_PATH
          )
find_path(Catch2_INCLUDE_DIR
          NAMES catch2/catch.hpp
          PATH_SUFFIXES single_include
          )
mark_as_advanced(Catch2_INCLUDE_DIR)

# process package
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Catch2 REQUIRED_VARS Catch2_INCLUDE_DIR)
set(Catch2_INCLUDE_DIRS ${Catch2_INCLUDE_DIR})

if(Catch2_FOUND AND NOT TARGET Catch2::Catch2)
    add_library(Catch2::Catch2 INTERFACE IMPORTED)
    set_target_properties(Catch2::Catch2 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Catch2_INCLUDE_DIRS}")
endif()
