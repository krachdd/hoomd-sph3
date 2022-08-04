# Maintainer: David Krach

##############################
## Source setup

# Set sources
set(_sph_sources module.cc)

# MODIFY THIS IF NEW CPP SOURCES ARE ADDED
set(_sph_cpu_sources 
                     SPHIntegratorTwoStep.cc
                     SPHIntegrationMethodTwoStep.cc
                     VelocityVerlet.cc
                     # SuspendedObjectIntegrator.cc
                     # RigidBodyIntegrator.cc
#                      StateEquationsExport.cc
                     SPHBaseClass.cc
#                      SinglePhaseFlow.cc
#                      TwoPhaseFlow.cc
#                      SmoothingKernelExport.cc
                     CustomForceCompute.cc
                   )

set(_sph_headers    #GeneralGPUDeviceFunctions.cuh  
#                     SinglePhaseFlowGPU.h  
                    SPHIntegrationMethodTwoStep.h     
#                     TwoPhaseFlowGPU.cuh
#                     GeneralGPUFunctions.cuh
#                     SinglePhaseFlow.h     
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

# if (NOT ENABLE_HIP)
# list(APPEND _sph_cpu_sources 
# #                              StateEquations.cc
# #                              SmoothingKernel.cc
# )
# endif()

# if (ENABLE_HIP)
# list(APPEND _sph_cpu_sources VelocityVerletGPU.cc
#                              # SuspendedObjectIntegratorGPU.cc
#                              # RigidBodyIntegratorGPU.cc
#                              SinglePhaseFlowGPU.cc
# #                              TwoPhaseFlowGPU.cc
#                            )
# endif()

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

add_subdirectory(data)


pybind11_add_module(_sph SHARED ${_sph_cpu_sources} #${_cuda_sources} 
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
set(files   
            __init__.py
            integrate.py
            force.py
#             kernel.py
#             eos.py
#             models.py
#             constrain.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/sph
       )

copy_files_to_build("${files}" "sph" "*.py")

# install headers in installation target
install(FILES ${_sph_headers}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include/hoomd/sph
       )
