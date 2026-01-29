/* ---------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

// Include the defined classes that are to be exported to python
// #include <hoomd/HOOMDMath.h>

// #include "StateEquations.h"
#include "SmoothingKernel.h"
#include "SPHBaseClass.h"
// #include "SPHIntegratorTwoStep.h"
// #include "SPHIntegrationMethodTwoStep.h"
#include "VelocityVerlet.h"
#include "VelocityVerletBasic.h"
#include "KickDriftKickTV.h"
// // #include "SuspendedObjectIntegrator.h"
// // #include "RigidBodyIntegrator.h"
#include "SinglePhaseFlow.h"
#include "SinglePhaseFlowTV.h"
#include "TwoPhaseFlow.h"
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

// #ifdef ENABLE_TBB
// #include <tbb/task_arena.h>
// #endif


namespace hoomd 
{
namespace sph
{
namespace detail
{



    void export_SPHIntegratorTwoStep(pybind11::module& m);
    void export_SPHIntegrationMethodTwoStep(pybind11::module& m);
    // void export_VelocityVerlet(pybind11::module& m);
    // void export_SuspendedObjectIntegrator(pybind11::module& m);
    // void export_RigidBodyIntegrator(pybind11::module& m);
    // void export_SinglePhaseFlow(pybind11::module& m);
    // void export_StateEquations(pybind11::module& m);
    // void export_TwoPhaseFlow(pybind11::module& m);
    // void export_SPHBaseClass(pybind11::module& m);
    void export_CustomForceCompute(pybind11::module& m);

    void export_WendlandC2(pybind11::module& m);
    void export_WendlandC4(pybind11::module& m);
    void export_WendlandC6(pybind11::module& m);
    void export_Quintic(pybind11::module& m);
    void export_CubicSpline(pybind11::module& m);

    // void export_StateEquation_Tait(pybind11::module& m);
    // void export_StateEquation_Linear(pybind11::module& m);

    void export_ComputeSPFMechanicalProperties(pybind11::module& m);
    void export_ComputeSolidProperties(pybind11::module& m);
    // void export_LocalNeighborListDataHost(pybind11::module& m);
    void export_HalfStepHook(pybind11::module& m);


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
    export_SPHIntegratorTwoStep(m);
    export_SPHIntegrationMethodTwoStep(m);
    export_VelocityVerlet(m);
    export_VelocityVerletBasic(m);
    export_KickDriftKickTV(m);
    // export_SuspendedObjectIntegrator(m);
    // export_RigidBodyIntegrator(m);
    export_WendlandC2(m);
    export_WendlandC4(m);
    export_WendlandC6(m);
    export_Quintic(m);
    export_CubicSpline(m);
    // export_SinglePhaseFlow(m);
    // export_StateEquation_Tait(m);
    // export_StateEquation_Linear(m);
    export_StateEquation<tait>(m, "Tait");
    export_StateEquation<linear>(m, "Linear");

    // export_TwoPhaseFlow(m);
    // export_SPHBaseClass(m);
    export_SPHBaseClass<wendlandc2, linear>(m, "SPHBaseClass_WC2_L");
    export_SPHBaseClass<wendlandc2, tait>(m, "SPHBaseClass_WC2_T");
    export_SPHBaseClass<wendlandc4, linear>(m, "SPHBaseClass_WC4_L");
    export_SPHBaseClass<wendlandc4, tait>(m, "SPHBaseClass_WC4_T");
    export_SPHBaseClass<wendlandc6, linear>(m, "SPHBaseClass_WC6_L");
    export_SPHBaseClass<wendlandc6, tait>(m, "SPHBaseClass_WC6_T");
    export_SPHBaseClass<quintic, linear>(m, "SPHBaseClass_Q_L");
    export_SPHBaseClass<quintic, tait>(m, "SPHBaseClass_Q_T");
    export_SPHBaseClass<cubicspline, linear>(m, "SPHBaseClass_CS_L");
    export_SPHBaseClass<cubicspline, tait>(m, "SPHBaseClass_CS_T");

    export_SinglePhaseFlow<wendlandc2, linear>(m, "SinglePF_WC2_L");
    export_SinglePhaseFlow<wendlandc2, tait>(m, "SinglePF_WC2_T");
    export_SinglePhaseFlow<wendlandc4, linear>(m, "SinglePF_WC4_L");
    export_SinglePhaseFlow<wendlandc4, tait>(m, "SinglePF_WC4_T");
    export_SinglePhaseFlow<wendlandc6, linear>(m, "SinglePF_WC6_L");
    export_SinglePhaseFlow<wendlandc6, tait>(m, "SinglePF_WC6_T");
    export_SinglePhaseFlow<quintic, linear>(m, "SinglePF_Q_L");
    export_SinglePhaseFlow<quintic, tait>(m, "SinglePF_Q_T");
    export_SinglePhaseFlow<cubicspline, linear>(m, "SinglePF_CS_L");
    export_SinglePhaseFlow<cubicspline, tait>(m, "SinglePF_CS_T");

    export_SinglePhaseFlowTV<wendlandc2, linear>(m, "SinglePFTV_WC2_L");
    export_SinglePhaseFlowTV<wendlandc2, tait>(m, "SinglePFTV_WC2_T");
    export_SinglePhaseFlowTV<wendlandc4, linear>(m, "SinglePFTV_WC4_L");
    export_SinglePhaseFlowTV<wendlandc4, tait>(m, "SinglePFTV_WC4_T");
    export_SinglePhaseFlowTV<wendlandc6, linear>(m, "SinglePFTV_WC6_L");
    export_SinglePhaseFlowTV<wendlandc6, tait>(m, "SinglePFTV_WC6_T");
    export_SinglePhaseFlowTV<quintic, linear>(m, "SinglePFTV_Q_L");
    export_SinglePhaseFlowTV<quintic, tait>(m, "SinglePFTV_Q_T");
    export_SinglePhaseFlowTV<cubicspline, linear>(m, "SinglePFTV_CS_L");
    export_SinglePhaseFlowTV<cubicspline, tait>(m, "SinglePFTV_CS_T");

    export_TwoPhaseFlow<wendlandc2, linear, linear>(m, "TwoPF_WC2_LL");
    export_TwoPhaseFlow<wendlandc2, linear, tait>(m, "TwoPF_WC2_LT");
    export_TwoPhaseFlow<wendlandc2, tait, linear>(m, "TwoPF_WC2_TL");
    export_TwoPhaseFlow<wendlandc2, tait, tait>(m, "TwoPF_WC2_TT");
    
    export_TwoPhaseFlow<wendlandc4, linear, linear>(m, "TwoPF_WC4_LL");
    export_TwoPhaseFlow<wendlandc4, linear, tait>(m, "TwoPF_WC4_LT");
    export_TwoPhaseFlow<wendlandc4, tait, linear>(m, "TwoPF_WC4_TL");
    export_TwoPhaseFlow<wendlandc4, tait, tait>(m, "TwoPF_WC4_TT");
    
    export_TwoPhaseFlow<wendlandc6, linear, linear>(m, "TwoPF_WC6_LL");
    export_TwoPhaseFlow<wendlandc6, linear, tait>(m, "TwoPF_WC6_LT");
    export_TwoPhaseFlow<wendlandc6, tait, linear>(m, "TwoPF_WC6_TL");
    export_TwoPhaseFlow<wendlandc6, tait, tait>(m, "TwoPF_WC6_TT");
    
    export_TwoPhaseFlow<quintic, linear, linear>(m, "TwoPF_Q_LL");
    export_TwoPhaseFlow<quintic, linear, tait>(m, "TwoPF_Q_LT");
    export_TwoPhaseFlow<quintic, tait, linear>(m, "TwoPF_Q_TL");
    export_TwoPhaseFlow<quintic, tait, tait>(m, "TwoPF_Q_TT");
    
    export_TwoPhaseFlow<cubicspline, linear, linear>(m, "TwoPF_CS_LL");
    export_TwoPhaseFlow<cubicspline, linear, tait>(m, "TwoPF_CS_LT");
    export_TwoPhaseFlow<cubicspline, tait, linear>(m, "TwoPF_CS_TL");
    export_TwoPhaseFlow<cubicspline, tait, tait>(m, "TwoPF_CS_TT");

    export_CustomForceCompute(m);

    export_ComputeSPFMechanicalProperties(m);
    export_ComputeSolidProperties(m);

    export_DensityMethod(m);
    export_ViscosityMethod(m);
    export_ColorGradientMethod(m);

    // export_LocalNeighborListDataHost(m);
    export_HalfStepHook(m);

// #ifdef ENABLE_HIP
//     export_VelocityVerletGPU(m);
//     //export_SuspendedObjectIntegratorGPU();
//     // export_RigidBodyIntegratorGPU(m);
//     export_SinglePhaseFlowGPU(m);
//     // export_TwoPhaseFlowGPU(m);
// #endif

}

