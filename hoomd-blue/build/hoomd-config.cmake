########################################################
# HOOMD CMake configuration for externally built plugins


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was hoomd-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

# build options
set(SINGLE_PRECISION "OFF")
set(ENABLE_HPMC_MIXED_PRECISION "ON")
set(HOOMD_GPU_PLATFORM "CUDA")

set(BUILD_MD "OFF")
set(BUILD_HPMC "OFF")
set(BUILD_METAL "OFF")
set(BUILD_DEM "")
set(BUILD_MPCD "OFF")

set(ENABLE_HIP "OFF")
set(HIP_PLATFORM "")
set(ENABLE_NVTOOLS "OFF")
set(ENABLE_ROCTRACER "")
set(ENABLE_MPI "ON")
set(ENABLE_MPI_CUDA "")
set(ENABLE_TBB "OFF")
set(ENABLE_LLVM "OFF")
set(ALWAYS_USE_MANAGED_MEMORY "OFF")

# CUDA architectures
set(CMAKE_CUDA_ARCHITECTURES "")

# C++ standard
set(CMAKE_CXX_STANDARD "17")
set(CMAKE_CUDA_STANDARD "14")

# installation locations
set(HOOMD_INSTALL_PREFIX "${PACKAGE_PREFIX_DIR}")
set(PYTHON_SITE_INSTALL_DIR "lib/python3.10/site-packages/hoomd")

# configure python
set(PYBIND11_PYTHON_VERSION 3)
find_package(pybind11 2.2 CONFIG REQUIRED)
find_package_message(pybind11 "Found pybind11: ${pybind11_DIR} ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION})" "[${pybind11_DIR}][${pybind11_INCLUDE_DIR}]")

find_package(Eigen3 3.2 CONFIG REQUIRED)
find_package_message(EIGEN3 "Found eigen: ${Eigen3_DIR} ${EIGEN3_INCLUDE_DIR} (version ${Eigen3_VERSION})" "[${Eigen3_DIR}][${EIGEN3_INCLUDE_DIR}]")

# find optional dependencies
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

if (ENABLE_HIP)
    include(HOOMDHIPSetup)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ")
endif()

if (ENABLE_MPI)
    find_dependency(MPI REQUIRED)

    find_package(cereal CONFIG)
    if (cereal_FOUND)
        find_package_message(cereal "Found cereal: ${cereal_DIR}" "[${cereal_DIR}]")

        if (NOT TARGET cereal::cereal AND TARGET cereal)
            message(STATUS "Found cereal target, adding cereal::cereal alias.")
            add_library(cereal::cereal ALIAS cereal)
        endif()
    else()
        # work around missing ceralConfig.cmake (common on Ubuntu 20.04)
        find_path(cereal_INCLUDE_DIR NAMES cereal/cereal.hpp
            PATHS ${CMAKE_INSTALL_PREFIX}/include)
        add_library(cereal::cereal INTERFACE IMPORTED)
        set_target_properties(cereal::cereal PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${cereal_INCLUDE_DIR}")
        find_package_message(cereal "Could not find cereal by config file, falling back to ${cereal_INCLUDE_DIR}" "[${cereal_INCLUDE_DIR}]")
    endif()

    # Work around broken cereal::cereal target (common on Ubuntu 22.04)
    get_target_property(_cereal_include cereal::cereal INTERFACE_INCLUDE_DIRECTORIES)
    if (_cereal_include STREQUAL "/include")
        find_path(cereal_INCLUDE_DIR NAMES cereal/cereal.hpp
            PATHS ${CMAKE_INSTALL_PREFIX}/include)
        set_target_properties(cereal::cereal PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${cereal_INCLUDE_DIR}")
        find_package_message(cereal "Fixing broken cereal::cereal target with ${cereal_INCLUDE_DIR}" "[${cereal_INCLUDE_DIR}]")
    endif()
endif()

if (ENABLE_TBB)
    find_dependency(TBB 4.3 REQUIRED)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/hoomd-targets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/hoomd-macros.cmake")

check_required_components(HOOMD)
