#-----------------------------------------------------------
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de
# ----------------------------------------------------------

set(COMPONENT_NAME sph)

set(_${COMPONENT_NAME}_sources
    module-sph.cc
    SPHIntegratorTwoStep.cc
    SPHIntegrationMethodTwoStep.cc
    VelocityVerlet.cc
    VelocityVerletBasic.cc
    KickDriftKickTV.cc
#     SuspendedObjectIntegrator.cc
#     RigidBodyIntegrator.cc
#     StateEquationsExport.cc
    StateEquations.cc
    SPHBaseClass.cc
    SinglePhaseFlow.cc
    SinglePhaseFlowTV.cc
#     TwoPhaseFlow.cc
#     SmoothingKernelExport.cc
    SmoothingKernel.cc
    CustomForceCompute.cc
    ComputeSPFBasicProperties.cc
    HalfStepHook.cc
    ConstantForceCompute.cc
    )

set(_${COMPONENT_NAME}_headers
#     GeneralGPUDeviceFunctions.cuh  
#     SinglePhaseFlowGPU.h  
    SPHIntegrationMethodTwoStep.h     
#     TwoPhaseFlowGPU.cuh
#     GeneralGPUFunctions.cuh
    SinglePhaseFlow.h     
    SinglePhaseFlowTV.h     
    SPHIntegratorTwoStep.h            
#     TwoPhaseFlowGPU.h
#     RigidBodyIntegratorGPU.cuh     
    SmoothingKernel.h     
    StateEquations.h                  
#     TwoPhaseFlow.h
#     RigidBodyIntegratorGPU.h       
#     SolidFluidTypeBit.h   
#     SuspendedObjectIntegratorGPU.cuh  
#     VelocityVerletGPU.cuh
#     RigidBodyIntegrator.h          
#     SPHBaseClass.cuh      
#     SuspendedObjectIntegratorGPU.h    
#     VelocityVerletGPU.h
#     SinglePhaseFlowGPU.cuh         
    SPHBaseClass.h        
#     SuspendedObjectIntegrator.h       
    VelocityVerlet.h
    VelocityVerletBasic.h
    KickDriftKickTV.h
    CustomForceCompute.h
    ComputeSPFBasicProperties.h
    ComputeMechanicalProp.h
    EvaluationMethodDefinition.h
    ConstantForceCompute.h
    )

set(_${COMPONENT_NAME}_cu_sources
#     VelocityVerletGPU.cu
#     SuspendedObjectIntegratorGPU.cu
#     RigidBodyIntegratorGPU.cu
#     SPHBaseClass.cu
#     SinglePhaseFlowGPU.cu
#     SmoothingKernel.cu
#     StateEquations.cu
#     TwoPhaseFlowGPU.cu
#     GeneralGPUFunctions.cu
#     GeneralGPUDeviceFunctions.cu
    )

if (ENABLE_HIP)
list (APPEND _${COMPONENT_NAME}_sources
#     VelocityVerletGPU.cc
#     SuspendedObjectIntegratorGPU.cc
#     RigidBodyIntegratorGPU.cc
#     SinglePhaseFlowGPU.cc
#     TwoPhaseFlowGPU.cc
    )
endif (ENABLE_HIP)

if (ENABLE_HIP)
set(_cuda_sources ${_md_cu_sources} ${DFFT_CU_SOURCES})
set_source_files_properties(${_md_cu_sources} PROPERTIES LANGUAGE ${HOOMD_DEVICE_LANGUAGE})
endif (ENABLE_HIP)

pybind11_add_module(_${COMPONENT_NAME} SHARED ${_${COMPONENT_NAME}_sources} ${_cuda_sources} NO_EXTRAS)
# alias into the HOOMD namespace so that plugins and symlinked components both work
add_library(HOOMD::_${COMPONENT_NAME} ALIAS _${COMPONENT_NAME})

if (APPLE)
set_target_properties(_${COMPONENT_NAME} PROPERTIES INSTALL_RPATH "@loader_path/..;@loader_path")
else()
set_target_properties(_${COMPONENT_NAME} PROPERTIES INSTALL_RPATH "\$ORIGIN/..;\$ORIGIN")
endif()

# link the library to its dependencies
target_link_libraries(_${COMPONENT_NAME} PUBLIC HOOMD::_hoomd)
target_link_libraries(_${COMPONENT_NAME} PUBLIC HOOMD::_nsearch)

# install the library
install(TARGETS _${COMPONENT_NAME}
        LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${COMPONENT_NAME}
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
set(files
    __init__.py
    integrate.py
    force.py
    kernel.py
    eos.py
    constrain.py
    compute.py
    half_step_hook.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${COMPONENT_NAME}
       )

copy_files_to_build("${files}" "${COMPONENT_NAME}" "*.py")

# if (BUILD_TESTING)
#     # add_subdirectory(test-py)
# endif()

add_subdirectory(data)

add_subdirectory(methods)
add_subdirectory(sphmodel)