# Install script for directory: /home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/USADR/ac130084/anaconda3/envs/sph3")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/include/hoomd/filter" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/export_filters.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterAll.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterCustom.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilter.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterIntersection.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterNull.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterRigid.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterSetDifference.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterTags.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterType.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/ParticleFilterUnion.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/filter" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/__init__.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/filter_.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/all_.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/null.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/rigid.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/set_.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/tags.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/type_.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/filter/custom.py"
    )
endif()

