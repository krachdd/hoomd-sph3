# Install script for directory: /home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so"
         RPATH "$ORIGIN/..:$ORIGIN")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph" TYPE SHARED_LIBRARY FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so"
         OLD_RPATH "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/nsearch:/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd:/usr/lib/x86_64-linux-gnu/openmpi/lib:"
         NEW_RPATH "$ORIGIN/..:$ORIGIN")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph/_sph.cpython-310-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/sph" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/__init__.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/integrate.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/force.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/kernel.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/eos.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/constrain.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/compute.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/sph/half_step_hook.py"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/sph/data/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/sph/methods/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/sph/sphmodel/cmake_install.cmake")

endif()

