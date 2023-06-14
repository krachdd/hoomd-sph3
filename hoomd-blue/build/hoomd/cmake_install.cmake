# Install script for directory: /home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/include" TYPE FILE FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/include/HOOMDVersion.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd" TYPE FILE FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/version_config.py")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so"
         RPATH "$ORIGIN")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd" TYPE SHARED_LIBRARY FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so"
         OLD_RPATH "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd:/usr/lib/x86_64-linux-gnu/openmpi/lib:"
         NEW_RPATH "$ORIGIN")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/_hoomd.cpython-310-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd" TYPE SHARED_LIBRARY FILES "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/libquickhull.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/libquickhull.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/box.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/communicator.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/_compile.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/conftest.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/device.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/__init__.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/error.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/operation.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/operations.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest_plugin_validate.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/util.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/variant.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/simulation.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/state.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/trigger.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/snapshot.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/logging.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/version.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/wall.py"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/pytest.ini"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python3.10/site-packages/hoomd/include/hoomd" TYPE FILE FILES
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/AABB.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/AABBTree.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Action.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Analyzer.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ArrayView.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Autotuned.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Autotuner.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BondedGroupData.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BondedGroupData.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BoxDim.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BoxResizeUpdater.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BoxResizeUpdaterGPU.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/BoxResizeUpdaterGPU.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/UpdaterRemoveDrift.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CachedAllocator.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CellListGPU.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CellListGPU.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CellList.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CellListStencil.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ClockSource.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CommunicatorGPU.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/CommunicatorGPU.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Communicator.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Compute.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/DomainDecomposition.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ExecutionConfiguration.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Filesystem.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ForceCompute.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ForceConstraint.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GlobalArray.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUArray.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUFlags.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUPartition.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUPolymorph.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUPolymorph.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GPUVector.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GSD.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GSDDumpWriter.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GSDReader.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/GSDShapeSpecWriter.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/HalfStepHook.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/HOOMDMath.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/HOOMDMPI.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Index1D.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Initializers.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Integrator.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Integrator.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/LoadBalancerGPU.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/LoadBalancerGPU.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/LoadBalancer.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/managed_allocator.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ManagedArray.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Messenger.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/MPIConfiguration.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ParticleData.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ParticleData.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ParticleGroup.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ParticleGroup.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/ParticleFilterUpdater.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/PythonLocalDataAccess.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/PythonUpdater.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/PythonAnalyzer.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/RandomNumbers.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/RNGIdentifiers.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SFCPackTunerGPU.cuh"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SFCPackTunerGPU.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SFCPackTuner.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SharedSignal.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SnapshotSystemData.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/SystemDefinition.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/System.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Trigger.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Tuner.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/TextureTools.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Updater.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/Variant.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/VectorMath.h"
    "/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/hoomd/WarpTools.cuh"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/extern/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/custom/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/data/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/filter/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/write/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/pytest/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/tune/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/update/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/nsearch/cmake_install.cmake")
  include("/home/USADR/ac130084/HOOMD-SPH/hoomd-sph3/hoomd-blue/build/hoomd/sph/cmake_install.cmake")

endif()

