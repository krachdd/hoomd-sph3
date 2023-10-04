

set(_nsearch_sources module-nsearch.cc
#                    ActiveForceCompute.cc
#                    ActiveRotationalDiffusionUpdater.cc
#                    AlchemostatTwoStep.cc
#                    AnisoPotentialPairALJ2.cc
#                    AnisoPotentialPairALJ3.cc
#                    AnisoPotentialPairDipole.cc
#                    AnisoPotentialPairGB.cc
#                    AlchemyData.cc
#                    BondTablePotential.cc
#                    CommunicatorGrid.cc
#                    ComputeThermo.cc
#                    ComputeThermoHMA.cc
#                    CosineSqAngleForceCompute.cc
#                    CustomForceCompute.cc
#                    EvaluatorWalls.cc
#                    FIREEnergyMinimizer.cc
#                    ForceComposite.cc
#                    ForceDistanceConstraint.cc
#                    HarmonicAngleForceCompute.cc
#                    HarmonicDihedralForceCompute.cc
#                    HarmonicImproperForceCompute.cc
#                    IntegrationMethodTwoStep.cc
#                    IntegratorTwoStep.cc
#                    ManifoldZCylinder.cc
#                    ManifoldDiamond.cc
#                    ManifoldEllipsoid.cc
#                    ManifoldGyroid.cc
#                    ManifoldXYPlane.cc
#                    ManifoldPrimitive.cc
#                    ManifoldSphere.cc
#                    MolecularForceCompute.cc
#                    MuellerPlatheFlow.cc
                   NeighborListBinned.cc
                   NeighborList.cc
                   NeighborListStencil.cc
                   NeighborListTree.cc
#                    OPLSDihedralForceCompute.cc
#                    PPPMForceCompute.cc
#                    TableAngleForceCompute.cc
#                    TableDihedralForceCompute.cc
#                    TwoStepBD.cc
#                    TwoStepBerendsen.cc
#                    TwoStepLangevinBase.cc
#                    TwoStepLangevin.cc
#                    TwoStepNPTMTK.cc
#                    TwoStepNVE.cc
#                    TwoStepNVTMTK.cc
#                    TwoStepNVTAlchemy.cc
#                    WallData.cc
#                    ZeroMomentumUpdater.cc
                   )

set(_nsearch_headers 
#                 ActiveForceComputeGPU.h
#                 ActiveForceCompute.h
#                 ActiveForceConstraintCompute.h
#                 ActiveForceConstraintComputeGPU.h
#                 ActiveForceConstraintComputeGPU.cuh
#                 ActiveRotationalDiffusionUpdater.h
#                 AlchemyData.h
#                 AnisoPotentialPairGPU.cuh
#                 AnisoPotentialPairGPU.h
#                 AnisoPotentialPair.h
#                 BondTablePotentialGPU.h
#                 BondTablePotential.h
#                 CommunicatorGridGPU.h
#                 CommunicatorGrid.h
#                 ComputeThermoGPU.cuh
#                 ComputeThermoGPU.h
#                 ComputeThermoHMAGPU.cuh
#                 ComputeThermoHMAGPU.h
#                 ComputeThermo.h
#                 ComputeThermoHMA.h
#                 ComputeThermoTypes.h
#                 ComputeThermoHMATypes.h
#                 CosineSqAngleForceComputeGPU.h
#                 CosineSqAngleForceCompute.h
#                 CustomForceCompute.h
#                 EvaluatorBondFENE.h
#                 EvaluatorBondHarmonic.h
#                 EvaluatorBondTether.h
#                 EvaluatorSpecialPairLJ.h
#                 EvaluatorSpecialPairCoulomb.h
#                 EvaluatorExternalElectricField.h
#                 EvaluatorExternalPeriodic.h
#                 EvaluatorPairALJ.h
#                 EvaluatorPairBuckingham.h
#                 EvaluatorPairDipole.h
#                 EvaluatorPairDLVO.h
#                 EvaluatorPairDPDThermoLJ.h
#                 EvaluatorPairDPDThermoDPD.h
#                 EvaluatorPairEwald.h
#                 EvaluatorPairForceShiftedLJ.h
#                 EvaluatorPairGauss.h
#                 EvaluatorPairGB.h
#                 EvaluatorPairLJ.h
#                 EvaluatorPairLJ1208.h
#                 EvaluatorPairLJ0804.h
#                 EvaluatorPairMie.h
#                 EvaluatorPairExpandedMie.h
#                 EvaluatorPairMoliere.h
#                 EvaluatorPairMorse.h
#                 EvaluatorPairOPP.h
#                 EvaluatorPairFourier.h
#                 EvaluatorPairReactionField.h
#                 EvaluatorPairExpandedLJ.h
#                 EvaluatorPairTable.h
#                 EvaluatorPairTWF.h
#                 EvaluatorPairYukawa.h
#                 EvaluatorPairZBL.h
#                 EvaluatorTersoff.h
#                 EvaluatorWalls.h
#                 FIREEnergyMinimizerGPU.h
#                 FIREEnergyMinimizer.h
#                 ForceCompositeGPU.h
#                 ForceComposite.h
#                 ForceDistanceConstraintGPU.h
#                 ForceDistanceConstraint.h
#                 HarmonicAngleForceComputeGPU.h
#                 HarmonicAngleForceCompute.h
#                 HarmonicDihedralForceComputeGPU.h
#                 HarmonicDihedralForceCompute.h
#                 HarmonicImproperForceComputeGPU.h
#                 HarmonicImproperForceCompute.h
#                 IntegrationMethodTwoStep.h
#                 IntegratorTwoStep.h
#                 ManifoldZCylinder.h
#                 ManifoldDiamond.h
#                 ManifoldEllipsoid.h
#                 ManifoldGyroid.h
#                 ManifoldXYPlane.h
#                 ManifoldPrimitive.h
#                 ManifoldSphere.h
#                 MolecularForceCompute.cuh
#                 MolecularForceCompute.h
#                 MuellerPlatheFlowEnum.h
#                 MuellerPlatheFlow.h
#                 MuellerPlatheFlowGPU.h
                NeighborListBinned.h
                NeighborListGPUBinned.h
                NeighborListGPU.h
                NeighborListGPUStencil.h
                NeighborListGPUTree.h
                NeighborList.h
                NeighborListStencil.h
                NeighborListTree.h
#                 OPLSDihedralForceComputeGPU.h
#                 OPLSDihedralForceCompute.h
#                 PotentialBondGPU.h
#                 PotentialBondGPU.cuh
#                 PotentialBond.h
#                 PotentialExternalGPU.h
#                 PotentialExternalGPU.cuh
#                 PotentialExternal.h
#                 PotentialPairAlchemical.h
#                 PotentialPairAlchemicalNormalized.h
#                 PotentialPairDPDThermoGPU.h
#                 PotentialPairDPDThermoGPU.cuh
#                 PotentialPairDPDThermo.h
#                 PotentialPairGPU.h
#                 PotentialPairGPU.cuh
#                 PotentialPair.h
#                 PotentialSpecialPairGPU.h
#                 PotentialSpecialPair.h
#                 PotentialTersoffGPU.h
#                 PotentialTersoff.h
#                 PPPMForceComputeGPU.h
#                 PPPMForceCompute.h
#                 TableAngleForceComputeGPU.h
#                 TableAngleForceCompute.h
#                 TableDihedralForceComputeGPU.h
#                 TableDihedralForceCompute.h
#                 TwoStepBDGPU.h
#                 TwoStepRATTLEBDGPU.h
#                 TwoStepRATTLEBDGPU.cuh
#                 TwoStepBD.h
#                 TwoStepRATTLEBD.h
#                 TwoStepBerendsenGPU.h
#                 TwoStepBerendsen.h
#                 TwoStepLangevinBase.h
#                 TwoStepLangevinGPU.h
#                 TwoStepRATTLELangevinGPU.h
#                 TwoStepRATTLELangevinGPU.cuh
#                 TwoStepLangevin.h
#                 TwoStepRATTLELangevin.h
#                 TwoStepNPTMTKGPU.h
#                 TwoStepNPTMTK.h
#                 TwoStepRATTLENVEGPU.h
#                 TwoStepRATTLENVEGPU.cuh
#                 TwoStepRATTLENVE.h
#                 TwoStepNVEGPU.h
#                 TwoStepNVE.h
#                 TwoStepNVTMTKGPU.h
#                 TwoStepNVTMTK.h
#                 AlchemostatTwoStep.h
#                 TwoStepNVTAlchemy.h
#                 WallData.h
#                 ZeroMomentumUpdater.h
                )

if (ENABLE_HIP)
list(APPEND _nsearch_sources 
#                             ActiveForceComputeGPU.cc
#                            AnisoPotentialPairALJ2GPU.cc
#                            AnisoPotentialPairALJ3GPU.cc
#                            AnisoPotentialPairDipoleGPU.cc
#                            AnisoPotentialPairGBGPU.cc
#                            BondTablePotentialGPU.cc
#                            CommunicatorGridGPU.cc
#                            ComputeThermoGPU.cc
#                            ComputeThermoHMAGPU.cc
#                            FIREEnergyMinimizerGPU.cc
#                            ForceCompositeGPU.cc
#                            ForceDistanceConstraintGPU.cc
#                            HarmonicAngleForceComputeGPU.cc
#                            HarmonicDihedralForceComputeGPU.cc
#                            HarmonicImproperForceComputeGPU.cc
#                            MolecularForceCompute.cu
                           NeighborListGPU.cc
                           NeighborListGPUBinned.cc
                           NeighborListGPUStencil.cc
                           NeighborListGPUTree.cc
#                            OPLSDihedralForceComputeGPU.cc
#                            PPPMForceComputeGPU.cc
#                            TableAngleForceComputeGPU.cc
#                            TableDihedralForceComputeGPU.cc
#                            TwoStepBDGPU.cc
#                            TwoStepBerendsenGPU.cc
#                            TwoStepLangevinGPU.cc
#                            TwoStepNPTMTKGPU.cc
#                            TwoStepNVEGPU.cc
#                            TwoStepNVTMTKGPU.cc
#                            MuellerPlatheFlowGPU.cc
#                            CosineSqAngleForceComputeGPU.cc
                           )
endif()

set(_nsearch_cu_sources 
#                     ActiveForceComputeGPU.cu
#                       AnisoPotentialPairALJ2GPUKernel.cu
#                       AnisoPotentialPairALJ3GPUKernel.cu
#                       AnisoPotentialPairDipoleGPUKernel.cu
#                       AnisoPotentialPairGBGPUKernel.cu
#                       ComputeThermoGPU.cu
#                       ComputeThermoHMAGPU.cu
#                       BondTablePotentialGPU.cu
#                       CommunicatorGridGPU.cu
#                       FIREEnergyMinimizerGPU.cu
#                       ForceCompositeGPU.cu
#                       ForceDistanceConstraintGPU.cu
#                       HarmonicAngleForceGPU.cu
#                       HarmonicDihedralForceGPU.cu
#                       HarmonicImproperForceGPU.cu
#                       MolecularForceCompute.cu
                      NeighborListGPUBinned.cu
                      NeighborListGPU.cu
                      NeighborListGPUStencil.cu
                      NeighborListGPUTree.cu
#                       OPLSDihedralForceGPU.cu
#                       PPPMForceComputeGPU.cu
#                       TableAngleForceGPU.cu
#                       TableDihedralForceGPU.cu
#                       TwoStepBDGPU.cu
#                       TwoStepBerendsenGPU.cu
#                       TwoStepLangevinGPU.cu
#                       TwoStepRATTLELangevinGPU.cu
#                       TwoStepNPTMTKGPU.cu
#                       TwoStepNVEGPU.cu
#                       TwoStepRATTLENVEGPU.cu
#                       TwoStepNVTMTKGPU.cu
#                       MuellerPlatheFlowGPU.cu
#                       CosineSqAngleForceGPU.cu
                      )

if (ENABLE_HIP)
set(_cuda_sources ${_md_cu_sources} ${DFFT_CU_SOURCES})
set_source_files_properties(${_md_cu_sources} PROPERTIES LANGUAGE ${HOOMD_DEVICE_LANGUAGE})
endif (ENABLE_HIP)

# generate pybind11 export cc files
set(_manifolds Cylinder
               Diamond
               Ellipsoid
               Gyroid
               Plane
               Primitive
               Sphere)

set(_rattle_methods BD Langevin NVE)

pybind11_add_module(_nsearch SHARED ${_nsearch_sources} ${_cuda_sources} ${_nsearch_headers} NO_EXTRAS)
# alias into the HOOMD namespace so that plugins and symlinked components both work
add_library(HOOMD::_nsearch ALIAS _nsearch)

if(APPLE)
set_target_properties(_nsearch PROPERTIES INSTALL_RPATH "@loader_path/..;@loader_path")
else()
set_target_properties(_nsearch PROPERTIES INSTALL_RPATH "\$ORIGIN/..;\$ORIGIN")
endif()

# link the library to its dependencies
if (CUSOLVER_AVAILABLE)
    # CUDA 8.0 requires that we link in gomp
    target_link_libraries(_nsearch PUBLIC _hoomd CUDA::cusolver CUDA::cusparse gomp)
else()
    target_link_libraries(_nsearch PUBLIC _hoomd)
endif()
if (ENABLE_HIP)
    target_link_libraries(_nsearch PRIVATE neighbor)
endif()

fix_cudart_rpath(_nsearch)

# install the library
install(TARGETS _nsearch EXPORT HOOMDTargets
        LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/nsearch
        )

################ Python only modules
# copy python modules to the build directory to make it a working python package
set(files __init__.py
#           angle.py
#           bond.py
#           compute.py
#           constrain.py
#           dihedral.py
#           force.py
          half_step_hook.py
#           improper.py
#           integrate.py
#           manifold.py
#           many_body.py
          nlist.py
#           update.py
#           special_pair.py
    )

install(FILES ${files}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/nsearch
       )

copy_files_to_build("${files}" "nsearch" "*.py")

# install headers in installation target
install(FILES ${_nsearch_headers}
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}/include/hoomd/nsearch
       )


