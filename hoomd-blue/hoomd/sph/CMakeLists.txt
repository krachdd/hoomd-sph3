# Maintainer: David Krach

# add_subdirectory(extern)

# configure the version information file
# configure_file (HOOMDVersion.h.inc ${HOOMD_BINARY_DIR}/hoomd/include/HOOMDVersion.h)

# install(FILES ${HOOMD_BINARY_DIR}/hoomd/include/HOOMDVersion.h
#         DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include
#         )

# translate cmake true/false to Python
# if (BUILD_MD)
#     set(_md_built "True")
# else()
#     set(_md_built "False")
# endif()

# if (BUILD_MD)
#     set(_md_built "True")
# else()
#     set(_md_built "False")
# endif()

# if (BUILD_HPMC)
#     set(_hpmc_built "True")
# else()
#     set(_hpmc_built "False")
# endif()

# if (BUILD_METAL)
#     set(_metal_built "True")
# else()
#     set(_metal_built "False")
# endif()

# if (BUILD_MPCD)
#     set(_mpcd_built "True")
# else()
#     set(_mpcd_built "False")
# endif()

# if (ENABLE_LLVM)
#     set(_llvm_enabled "True")
# else()
#     set(_llvm_enabled "False")
# endif()

# configure_file (version_config.py.in ${HOOMD_BINARY_DIR}/hoomd/version_config.py)
# install(FILES ${HOOMD_BINARY_DIR}/hoomd/version_config.py
#         DESTINATION ${PYTHON_SITE_INSTALL_DIR}
#         )

##############################
## Source setup

# Set sources
set(_sph_sources module.cc)

# MODIFY THIS IF NEW CPP SOURCES ARE ADDED
set(_sph_cpu_sources SPHIntegratorTwoStep.cc
                     SPHIntegrationMethodTwoStep.cc
                     VelocityVerlet.cc
                     # SuspendedObjectIntegrator.cc
                     # RigidBodyIntegrator.cc
#                      StateEquationsExport.cc
                     SPHBaseClass.cc
                     SinglePhaseFlow.cc
#                      TwoPhaseFlow.cc
                     SmoothingKernelExport.cc
                     CustomForceCompute.cc
                   )

# ignore conversion warnings in external files
# if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     set_source_files_properties(extern/kiss_fft.cc PROPERTIES COMPILE_FLAGS "-Wno-conversion -Wno-float-conversion")
#     set_source_files_properties(extern/vmdsock.cc PROPERTIES COMPILE_FLAGS "-Wno-conversion")
#     set_source_files_properties(extern/imd.cc PROPERTIES COMPILE_FLAGS "-Wno-conversion")
# endif()

set(_sph_headers    #GeneralGPUDeviceFunctions.cuh  
#                     SinglePhaseFlowGPU.h  
                    SPHIntegrationMethodTwoStep.h     
#                     TwoPhaseFlowGPU.cuh
#                     GeneralGPUFunctions.cuh
                    SinglePhaseFlow.h     
                    SPHIntegratorTwoStep.h            
#                     TwoPhaseFlowGPU.h
                    # RigidBodyIntegratorGPU.cuh     
                    SmoothingKernel.h     
                    StateEquations.h                  
#                     TwoPhaseFlow.h
                    # RigidBodyIntegratorGPU.h       
                    SolidFluidTypeBit.h   
                    # SuspendedObjectIntegratorGPU.cuh  
#                     VelocityVerletGPU.cuh
                    # RigidBodyIntegrator.h          
#                     SPHBaseClass.cuh      
                    # SuspendedObjectIntegratorGPU.h    
#                     VelocityVerletGPU.h
#                     SinglePhaseFlowGPU.cuh         
                    SPHBaseClass.h        
#                     SuspendedObjectIntegrator.h       
                    VelocityVerlet.h
                    CustomForceCompute.h
    )

if (NOT ENABLE_HIP)
list(APPEND _sph_cpu_sources StateEquations.cc
                             SmoothingKernel.cc)
endif()

if (ENABLE_HIP)
list(APPEND _sph_cpu_sources VelocityVerletGPU.cc
                             # SuspendedObjectIntegratorGPU.cc
                             # RigidBodyIntegratorGPU.cc
                             SinglePhaseFlowGPU.cc
#                              TwoPhaseFlowGPU.cc
                           )
endif()

# set(_sph_cu_sources VelocityVerletGPU.cu
#                     #SuspendedObjectIntegratorGPU.cu
#                     # RigidBodyIntegratorGPU.cu
#                     SPHBaseClass.cu
#                     SinglePhaseFlowGPU.cu
#                     SmoothingKernel.cu
#                     StateEquations.cu
# #                     TwoPhaseFlowGPU.cu
#                     GeneralGPUFunctions.cu
#                     GeneralGPUDeviceFunctions.cu
#                     )

# if (ENABLE_HIP)
# set(_cuda_sources ${_sph_cu_sources})
# endif (ENABLE_HIP)

# #########################
# ## Build the module
# pybind11_add_module(_sph SHARED module.cc ${_sph_cpu_sources} ${_cuda_sources} ${_sph_headers} NO_EXTRAS)

# # alias into the HOOMD namespace so that plugins and symlinked components both work
# add_library(HOOMD::_sph ALIAS _sph)

# # Work around support for the delete operator with pybind11 and older versions of clang
# # https://github.com/pybind/pybind11/issues/1604
# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#     if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.11)
#         target_compile_options(_sph PUBLIC $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<STREQUAL:${HIP_PLATFORM},nvcc>>:-Xcompiler=>;-fsized-deallocation)
#     else()
#         target_compile_options(_sph PUBLIC -fsized-deallocation)
#     endif()
# endif()

# # add  hull as its own library so that it's symbols can be public
# # add_library (quickhull SHARED extern/quickhull/QuickHull.cpp)

# # ignore conversion warnings in quickull
# if(CMAKE_COMPILER_IS_GNUCXX)
#     if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.0.0)
#         target_compile_options(quickhull PRIVATE "-Wno-conversion;-Wno-pessimizing-move")
#     else()
#         target_compile_options(quickhull PRIVATE "-Wno-conversion")
#     endif()
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#     target_compile_options(quickhull PRIVATE "-Wno-conversion;-Wno-pessimizing-move")
# endif()

# # link the library to its dependencies
# target_link_libraries(_sph PUBLIC pybind11::pybind11 quickhull Eigen3::Eigen)
# if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
#     target_link_libraries(_sph PUBLIC execinfo) # on FreeBSD backtrace() is in libexecinfo
# endif()

# # specify required include directories
# target_include_directories(_sph PUBLIC
#                                   $<BUILD_INTERFACE:${HOOMD_SOURCE_DIR}>
#                                   $<INSTALL_INTERFACE:${PYTHON_SITE_INSTALL_DIR}/include>)

# target_include_directories(_sph PUBLIC
#                                   $<BUILD_INTERFACE:${HOOMD_BINARY_DIR}/hoomd/include>)

# # specify required definitions
# target_compile_definitions(_sph PUBLIC _REENTRANT EIGEN_MPL2_ONLY)

# if (SINGLE_PRECISION)
#     target_compile_definitions(_sph PUBLIC SINGLE_PRECISION)
# endif(SINGLE_PRECISION)

# # Libraries and compile definitions for CUDA enabled builds
# # if (ENABLE_HIP)
# #     if (HIP_PLATFORM STREQUAL "hcc" OR HIP_PLATFORM STREQUAL "hip-clang")
# #         target_link_libraries(_sph PUBLIC HIP::hipfft HIP::hiprt)
# #     elseif(HIP_PLATFORM STREQUAL "nvcc")
# #         target_link_libraries(_sph PUBLIC CUDA::cudart CUDA::cufft)
# #     endif()
# #     target_compile_definitions(_sph PUBLIC ENABLE_HIP CUDA_ARCH=${_cuda_min_arch})

# #     if(ALWAYS_USE_MANAGED_MEMORY)
# #         target_compile_definitions(_sph PUBLIC ALWAYS_USE_MANAGED_MEMORY)
# #     endif()

# #     if (ENABLE_NVTOOLS)
# #         target_link_libraries(_sph PUBLIC CUDA::nvToolsExt)
# #         target_compile_definitions(_sph PUBLIC ENABLE_NVTOOLS)
# #     endif()

# #     if (CUSOLVER_AVAILABLE)
# #         target_compile_definitions(_sph PUBLIC CUSOLVER_AVAILABLE)
# #     endif()

# #     target_link_libraries(_sph PUBLIC HIP::hip)

# #     if (ENABLE_ROCTRACER)
# #         target_link_libraries(_sph PUBLIC HIP::roctracer)
# #         target_compile_definitions(_sph PUBLIC ENABLE_ROCTRACER)
# #     endif()
# # endif()


# # Libraries and compile definitions for TBB enabled builds
# if (ENABLE_TBB)
#     find_package_config_first(TBB 4.3)

#     if (TBB_FOUND)
#         get_target_property(_tbb_library TBB::tbb IMPORTED_LOCATION_RELEASE)
#         get_target_property(_tbb_include_dir TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
#         find_package_message(tbb "Found TBB: ${TBB_DIR} ${_tbb_library} ${_tbb_include_dir}" "[${_tbb_library}][${_tbb_include_dir}]")
#     endif()

#     target_compile_definitions(_sph PUBLIC ENABLE_TBB)
#     target_link_libraries(_sph PUBLIC TBB::tbb)
# endif()

# # Libraries and compile definitions for MPI enabled builds
# if (ENABLE_MPI)
#     target_compile_definitions(_sph PUBLIC ENABLE_MPI)
#     target_link_libraries(_sph PUBLIC MPI::MPI_CXX cereal::cereal)

#     if (ENABLE_MPI_CUDA)
#           target_compile_definitions(_sph PUBLIC ENABLE_MPI_CUDA)
#     endif(ENABLE_MPI_CUDA)
# endif()

# if (ENABLE_HPMC_MIXED_PRECISION)
#     target_compile_definitions(_sph PUBLIC ENABLE_HPMC_MIXED_PRECISION)
# endif()

# if (APPLE)
# set_target_properties(_sph PROPERTIES INSTALL_RPATH "@loader_path")
# else()
# set_target_properties(_sph PROPERTIES INSTALL_RPATH "\$ORIGIN")
# endif()

# fix_cudart_rpath(_sph)

# # install the library
# install(TARGETS _sph quickhull EXPORT HOOMDTargets
#         LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph
#         )

# ################ Python only modules
# # copy python modules to the build directory to make it a working python package
# set(files   __init__.py
#             integrate.py
#             force.py
#             kernel.py
#             eos.py
#             models.py
#     )

# install(FILES ${files}
#         DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph
#        )

# # version_config.py is generated by configure_file, ignore here
# copy_files_to_build("${files}" "sph" "*.py" "version_config.py")

# # install headers in installation target
# install(FILES ${_sph_headers}
#         DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include/hoomd/sph
#        )

# subdirectories that are not components
# add_subdirectory(custom)
add_subdirectory(data)
# add_subdirectory(filter)
# add_subdirectory(write)
# add_subdirectory(pytest)
# add_subdirectory(tune)
# add_subdirectory(update)

# if (BUILD_TESTING)
#     # add_subdirectory(test-py)
#     add_subdirectory(test)
# endif()

##################################################
## Build components

# if (BUILD_MD)
#     if (ENABLE_MPI)
#         # add the distributed FFT library
#         add_subdirectory(${HOOMD_SOURCE_DIR}/hoomd/extern/dfftlib)
#     endif()
#     add_subdirectory(md)
# endif()

# if (NOT SINGLE_PRECISION AND BUILD_HPMC)
#     add_subdirectory(hpmc)
# endif()

# if (BUILD_METAL AND BUILD_MD)
#     add_subdirectory(metal)
# endif()

# if (NOT ENABLE_HIP OR HIP_PLATFORM STREQUAL "nvcc")
#     if (BUILD_MPCD)
#         add_subdirectory(mpcd)
#     endif()
# endif()

# foreach(entry ${PLUGINS})
#     if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${entry}/CMakeLists.txt)
#         add_subdirectory(${entry})
#     endif()
# endforeach()

#option(BUILD_SPH "Build the sph package" on)
#if (BUILD_SPH)
#    add_subdirectory(sph)
#endif()


# if (ENABLE_HIP)
# set(_cuda_sources ${_sph_cu_sources} )
# endif (ENABLE_HIP)



pybind11_add_module(_sph SHARED ${_sph_sources} #${_cuda_sources} 
    ${_sph_headers} NO_EXTRAS)
# alias into the HOOMD namespace so that plugins and symlinked components both work
add_library(HOOMD::_sph ALIAS _sph)

if(APPLE)
set_target_properties(_sph PROPERTIES INSTALL_RPATH "@loader_path/..;@loader_path")
else()
set_target_properties(_sph PROPERTIES INSTALL_RPATH "\$ORIGIN/..;\$ORIGIN")
endif()

# link the library to its dependencies
if (CUSOLVER_AVAILABLE)
    # CUDA 8.0 requires that we link in gomp
    target_link_libraries(_sph PUBLIC _hoomd CUDA::cusolver CUDA::cusparse gomp)
else()
    target_link_libraries(_sph PUBLIC _hoomd)
endif()
if (ENABLE_HIP)
    target_link_libraries(_sph PRIVATE neighbor)
endif()

fix_cudart_rpath(_sph)

# install the library
install(TARGETS _sph EXPORT HOOMDTargets
        LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
set(files   __init__.py
            integrate.py
            force.py
            kernel.py
            eos.py
            models.py
            constrain.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph
       )

copy_files_to_build("${files}" "sph" "*.py")

# install headers in installation target
install(FILES ${_sph_headers}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include/hoomd/sph
       )
