// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Maintainer: David Krach

// Include the defined classes that are to be exported to python
// #include <hoomd/HOOMDMath.h>

// #include "StateEquations.h"
// #include "SmoothingKernel.h"
// #include "SPHBaseClass.h"
// #include "SPHIntegratorTwoStep.h"
// #include "SPHIntegrationMethodTwoStep.h"
// #include "VelocityVerlet.h"
// // #include "SuspendedObjectIntegrator.h"
// // #include "RigidBodyIntegrator.h"
// #include "SinglePhaseFlow.h"
// // #include "TwoPhaseFlow.h"

// // include GPU classes
// #ifdef ENABLE_HIP
// #include "VelocityVerletGPU.h"
// //#include "SuspendedObjectIntegratorGPU.h"
// // #include "RigidBodyIntegratorGPU.h"
// #include "SinglePhaseFlowGPU.h"
// // #include "TwoPhaseFlowGPU.h"
// #endif

// // ParticleFilter objects
// #include "hoomd/filter/export_filters.h"

#include <pybind11/pybind11.h>
// #include <pybind11/stl_bind.h>

// #include <fstream>
// #include <iostream>
// #include <sstream>

// #ifdef ENABLE_TBB
// #include <tbb/task_arena.h>
// #endif


// Include boost.python to do the exporting
// #include <boost/python.hpp>
// using namespace boost::python;

namespace hoomd 
{
namespace sph
{
namespace detail
{



    void export_SPHIntegratorTwoStep(pybind11::module& m);
    void export_SPHIntegrationMethodTwoStep(pybind11::module& m);
    void export_VelocityVerlet(pybind11::module& m);
    // void export_SuspendedObjectIntegrator(pybind11::module& m);
    // void export_RigidBodyIntegrator(pybind11::module& m);
    void export_WendlandC2(pybind11::module& m);
    void export_WendlandC4(pybind11::module& m);
    void export_WendlandC6(pybind11::module& m);
    void export_Quintic(pybind11::module& m);
    void export_CubicSpline(pybind11::module& m);
    void export_SinglePhaseFlow(pybind11::module& m);
    void export_StateEquations(pybind11::module& m);
    // void export_TwoPhaseFlow(pybind11::module& m);
    void export_SPHBasesClass(pybind11::module& m);

#ifdef ENABLE_HIP
    void export_VelocityVerletGPU(pybind11::module& m);
    //void export_SuspendedObjectIntegratorGPU();
    // void export_RigidBodyIntegratorGPU(pybind11::module& m);
    void export_SinglePhaseFlowGPU(pybind11::module& m);
    // void export_TwoPhaseFlowGPU(pybind11::module& m);
#endif



// PYBIND11_MODULE(_hoomd, m){
//     export_SPHIntegratorTwoStep(m);
//     export_SPHIntegrationMethodTwoStep(m);
//     export_VelocityVerlet(m);
//     // export_SuspendedObjectIntegrator(m);
//     // export_RigidBodyIntegrator(m);
//     export_WendlandC2(m);
//     export_WendlandC4(m);
//     export_WendlandC6(m);
//     export_Quintic(m);
//     export_CubicSpline(m);
//     export_SinglePhaseFlow(m);
//     export_StateEquations(m);
//     // export_TwoPhaseFlow(m);
//     export_SPHBasesClass(m);

//     #ifdef ENABLE_HIP
//     export_VelocityVerletGPU(m);
//     //export_SuspendedObjectIntegratorGPU();
//     // export_RigidBodyIntegratorGPU(m);
//     export_SinglePhaseFlowGPU(m);
//     // export_TwoPhaseFlowGPU(m);
//     #endif
// }


} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 


// BOOST_PYTHON_MODULE(_sph)
//     {
//     export_SPHIntegratorTwoStep();
//     export_SPHIntegrationMethodTwoStep();
//     export_VelocityVerlet();
//     export_SuspendedObjectIntegrator();
//     export_RigidBodyIntegrator();
//     export_WendlandC2();
//     export_WendlandC4();
//     export_WendlandC6();
//     export_Quintic();
//     export_CubicSpline();
//     export_SinglePhaseFlow();
//     export_StateEquations();
//     export_TwoPhaseFlow();
//     export_SPHBasesClass();
// #ifdef ENABLE_HIP
//     export_VelocityVerletGPU();
//     //export_SuspendedObjectIntegratorGPU();
//     export_RigidBodyIntegratorGPU();
//     export_SinglePhaseFlowGPU();
//     export_TwoPhaseFlowGPU();
// #endif

//     // boost 1.60.0 compatibility
//     #if (BOOST_VERSION == 106000)
//     register_ptr_to_python< boost::shared_ptr< SPHIntegratorTwoStep > >();
//     register_ptr_to_python< boost::shared_ptr< SPHIntegrationMethodTwoStep > >();
//     register_ptr_to_python< boost::shared_ptr< VelocityVerlet > >();
//     register_ptr_to_python< boost::shared_ptr< SuspendedObjectIntegrator > >();
//     register_ptr_to_python< boost::shared_ptr< RigidBodyIntegrator > >();
//     register_ptr_to_python< boost::shared_ptr< WendlandC2 > >();
//     register_ptr_to_python< boost::shared_ptr< WendlandC4 > >();
//     register_ptr_to_python< boost::shared_ptr< WendlandC6 > >();
//     register_ptr_to_python< boost::shared_ptr< Quintic > >();
//     register_ptr_to_python< boost::shared_ptr< CubicSpline > >();
//     register_ptr_to_python< boost::shared_ptr< SinglePhaseFlow > >();
//     register_ptr_to_python< boost::shared_ptr< TwoPhaseFlow > >();
//     register_ptr_to_python< boost::shared_ptr< StateEquation > >();
//     register_ptr_to_python< boost::shared_ptr< Tait > >();
//     register_ptr_to_python< boost::shared_ptr< Linear > >();

//     #ifdef ENABLE_HIP
//     register_ptr_to_python< boost::shared_ptr< VelocityVerletGPU > >();
//     //register_ptr_to_python< boost::shared_ptr< SuspendedObjectIntegratorGPU > >();
//     register_ptr_to_python< boost::shared_ptr< RigidBodyIntegratorGPU > >();
//     register_ptr_to_python< boost::shared_ptr< SinglePhaseFlowGPU > >();
//     #endif

//     #endif

//     }
