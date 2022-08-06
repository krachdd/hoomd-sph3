// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Maintainer: David Krach

// Include the defined classes that are to be exported to python
// #include <hoomd/HOOMDMath.h>

#include "StateEquations.h"
#include "SmoothingKernel.h"
#include "SPHBaseClass.h"
// #include "SPHIntegratorTwoStep.h"
// #include "SPHIntegrationMethodTwoStep.h"
// #include "VelocityVerlet.h"
// // #include "SuspendedObjectIntegrator.h"
// // #include "RigidBodyIntegrator.h"
#include "SinglePhaseFlow.h"
// // #include "TwoPhaseFlow.h"
#include "CustomForceCompute.h"

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



    // void export_SPHIntegratorTwoStep(pybind11::module& m);
    // void export_SPHIntegrationMethodTwoStep(pybind11::module& m);
    // void export_VelocityVerlet(pybind11::module& m);
    // void export_SuspendedObjectIntegrator(pybind11::module& m);
    // void export_RigidBodyIntegrator(pybind11::module& m);
    // void export_WendlandC2(pybind11::module& m);
    // void export_WendlandC4(pybind11::module& m);
    // void export_WendlandC6(pybind11::module& m);
    // void export_Quintic(pybind11::module& m);
    // void export_CubicSpline(pybind11::module& m);
    // void export_SinglePhaseFlow(pybind11::module& m);
    // void export_StateEquations(pybind11::module& m);
    // void export_TwoPhaseFlow(pybind11::module& m);
    // void export_SPHBaseClass(pybind11::module& m);
    // void export_CustomForceCompute(pybind11::module& m);

// #ifdef ENABLE_HIP
//     void export_VelocityVerletGPU(pybind11::module& m);
//     //void export_SuspendedObjectIntegratorGPU();
//     // void export_RigidBodyIntegratorGPU(pybind11::module& m);
//     void export_SinglePhaseFlowGPU(pybind11::module& m);
//     // void export_TwoPhaseFlowGPU(pybind11::module& m);
// #endif


} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 


using namespace hoomd;
using namespace hoomd::sph;
using namespace hoomd::sph::detail;


PYBIND11_MODULE(_sph, m){
    // export_SPHIntegratorTwoStep(m);
    // export_SPHIntegrationMethodTwoStep(m);
    // export_VelocityVerlet(m);
    // export_SuspendedObjectIntegrator(m);
    // export_RigidBodyIntegrator(m);
    // export_WendlandC2(m);
    // export_WendlandC4(m);
    // export_WendlandC6(m);
    // export_Quintic(m);
    // export_CubicSpline(m);
    // export_SinglePhaseFlow(m);
    // export_StateEquation_tait(m);
    // export_StateEquation_linear(m);
    // export_TwoPhaseFlow(m);
    // export_SPHBaseClass(m);
    // export_SPHBaseClass<wendlandc2, linear>(m, "SPHBaseClass_WC2_L");
    // export_SPHBaseClass<wendlandc2, tait>(m, "SPHBaseClass_WC2_T");
    // export_SPHBaseClass<wendlandc4, linear>(m, "SPHBaseClass_WC4_L");
    // export_SPHBaseClass<wendlandc4, tait>(m, "SPHBaseClass_WC4_T");
    // export_SPHBaseClass<wendlandc6, linear>(m, "SPHBaseClass_WC6_L");
    // export_SPHBaseClass<wendlandc6, tait>(m, "SPHBaseClass_WC6_T");
    // export_SPHBaseClass<quintic, linear>(m, "SPHBaseClass_Q_L");
    // export_SPHBaseClass<quintic, tait>(m, "SPHBaseClass_Q_T");
    // export_SPHBaseClass<cubicspline, linear>(m, "SPHBaseClass_CS_L");
    // export_SPHBaseClass<cubicspline, tait>(m, "SPHBaseClass_CS_T");
    // export_SinglePhaseFlow<wendlandc2, linear>(m, "SinglePF_WC2_L");
    // export_SinglePhaseFlow<wendlandc2, tait>(m, "SinglePF_WC2_T");
    // export_SinglePhaseFlow<wendlandc4, linear>(m, "SinglePF_WC4_L");
    // export_SinglePhaseFlow<wendlandc4, tait>(m, "SinglePF_WC4_T");
    // export_SinglePhaseFlow<wendlandc6, linear>(m, "SinglePF_WC6_L");
    // export_SinglePhaseFlow<wendlandc6, tait>(m, "SinglePF_WC6_T");
    // export_SinglePhaseFlow<quintic, linear>(m, "SinglePF_Q_L");
    // export_SinglePhaseFlow<quintic, tait>(m, "SinglePF_Q_T");
    // export_SinglePhaseFlow<cubicspline, linear>(m, "SinglePF_CS_L");
    // export_SinglePhaseFlow<cubicspline, tait>(m, "SinglePF_CS_T");
    export_CustomForceCompute(m);

// #ifdef ENABLE_HIP
//     export_VelocityVerletGPU(m);
//     //export_SuspendedObjectIntegratorGPU();
//     // export_RigidBodyIntegratorGPU(m);
//     export_SinglePhaseFlowGPU(m);
//     // export_TwoPhaseFlowGPU(m);
// #endif

}

