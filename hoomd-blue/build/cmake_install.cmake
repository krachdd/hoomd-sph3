# Install script for directory: /home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd" TYPE FILE FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd-config-version.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd/hoomd-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd/hoomd-targets.cmake"
         "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/CMakeFiles/Export/lib/cmake/hoomd/hoomd-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd/hoomd-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd/hoomd-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd" TYPE FILE FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/CMakeFiles/Export/lib/cmake/hoomd/hoomd-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd" TYPE FILE FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/CMakeFiles/Export/lib/cmake/hoomd/hoomd-targets-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/hoomd" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/CMake/hoomd/FindTBB.cmake"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/CMake/hoomd/FindCUDALibs.cmake"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/CMake/hoomd/HOOMDHIPSetup.cmake"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/CMake/hoomd/hoomd-macros.cmake"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd-config.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/CMake/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
