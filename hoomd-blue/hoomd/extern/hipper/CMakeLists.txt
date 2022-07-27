# Copyright (c) 2020, Michael P. Howard
# This file is released under the Modified BSD License.

# Maintainer: mphoward

cmake_minimum_required(VERSION 3.1)

project(hipper VERSION 0.2.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/tools")

# check if hipper is being used directly or via add_subdirectory
set(HIPPER_MASTER_PROJECT OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(HIPPER_MASTER_PROJECT ON)
endif()
option(HIPPER_INSTALL "Install hipper files?" ${HIPPER_MASTER_PROJECT})
option(HIPPER_TEST "Build hipper tests?" ${HIPPER_MASTER_PROJECT})

# header-only targets
include(GNUInstallDirs)
add_library(hipper INTERFACE)
add_library(hipper::hipper ALIAS hipper)
target_include_directories(hipper INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# add tests
if(HIPPER_TEST)
    enable_testing()
    add_subdirectory(tests)
endif()

# install headers
if(HIPPER_INSTALL)
    include(CMakePackageConfigHelpers)
    install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/hipper DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    configure_package_config_file(tools/hipperConfig.cmake.in
                                  "${CMAKE_CURRENT_BINARY_DIR}/hipperConfig.cmake"
                                  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/hipper)

    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/hipperConfigVersion.cmake
                                     VERSION ${hipper_VERSION}
                                     COMPATIBILITY SameMajorVersion)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hipperConfig.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/hipperConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/hipper)

    install(TARGETS hipper EXPORT hipperTargets)

    install(EXPORT hipperTargets
            NAMESPACE hipper::
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/hipper)
endif()
