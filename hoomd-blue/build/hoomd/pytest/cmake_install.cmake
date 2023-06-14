# Install script for directory: /home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/pytest" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/__init__.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_attr_tuner.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_balance.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_box.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_box_resize.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_collections.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_communicator.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_custom_tuner.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_custom_updater.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_custom_writer.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_dcd.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_device.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_filter_updater.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_mesh.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_trigger.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_type_parameter_dict.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_typeparam.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_operation.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_remove_drift.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_syncedlist.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_local_snapshot.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_logging.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_filter.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_parameter_dict.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/dummy.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_snapshot.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_state.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_simulation.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_table.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_tune_solve.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_variant.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_sorter.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/test_operations.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/pytest" TYPE PROGRAM FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest/pytest-openmpi.sh")
endif()

