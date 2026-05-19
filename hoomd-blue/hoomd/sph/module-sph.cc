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
#include "SinglePhaseFlowGDGD.h"
#include "SinglePhaseFlowFS.h"
#include "TwoPhaseFlow.h"
#include "TwoPhaseFlowTV.h"
#include "CustomForceCompute.h"


// include GPU classes
#ifdef ENABLE_HIP
#include "ComputeSPFBasicPropertiesGPU.h"
#include "SinglePhaseFlowGPU.h"
#include "SinglePhaseFlowTVGPU.h"
#include "SinglePhaseFlowGDGDGPU.h"
#include "SinglePhaseFlowFSGPU.h"
#include "TwoPhaseFlowGPU.h"
#include "TwoPhaseFlowTVGPU.h"
#include "VelocityVerletGPU.h"
#endif

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
    // export_SinglePhaseFlowFS declared via template in SinglePhaseFlowFS.cc

    void export_WendlandC2(pybind11::module& m);
    void export_WendlandC4(pybind11::module& m);
    void export_WendlandC6(pybind11::module& m);
    void export_Quintic(pybind11::module& m);
    void export_CubicSpline(pybind11::module& m);

    // void export_StateEquation_Tait(pybind11::module& m);
    // void export_StateEquation_Linear(pybind11::module& m);

    void export_ComputeSPFMechanicalProperties(pybind11::module& m);
    void export_ComputeSolidProperties(pybind11::module& m);

#ifdef ENABLE_HIP
    void export_ComputeSPFMechanicalPropertiesGPU(pybind11::module& m);
    void export_VelocityVerletGPU(pybind11::module& m);
#endif
    // void export_LocalNeighborListDataHost(pybind11::module& m);
    void export_HalfStepHook(pybind11::module& m);


#ifdef ENABLE_HIP
    template<SmoothingKernelType KT_, StateEquationType SET_>
    void export_SinglePhaseFlowGPU(pybind11::module& m, std::string name);

    template<SmoothingKernelType KT_, StateEquationType SET_>
    void export_SinglePhaseFlowTVGPU(pybind11::module& m, std::string name);

    template<SmoothingKernelType KT_, StateEquationType SET_>
    void export_SinglePhaseFlowGDGDGPU(pybind11::module& m, std::string name);

    template<SmoothingKernelType KT_, StateEquationType SET_>
    void export_SinglePhaseFlowFSGPU(pybind11::module& m, std::string name);

    template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
    void export_TwoPhaseFlowGPU(pybind11::module& m, std::string name);

    template<SmoothingKernelType KT_, StateEquationType SET1_, StateEquationType SET2_>
    void export_TwoPhaseFlowTVGPU(pybind11::module& m, std::string name);
#endif


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

    export_SinglePhaseFlowGDGD<wendlandc2, linear>(m, "SinglePFGDGD_WC2_L");
    export_SinglePhaseFlowGDGD<wendlandc2, tait>  (m, "SinglePFGDGD_WC2_T");
    export_SinglePhaseFlowGDGD<wendlandc4, linear>(m, "SinglePFGDGD_WC4_L");
    export_SinglePhaseFlowGDGD<wendlandc4, tait>  (m, "SinglePFGDGD_WC4_T");
    export_SinglePhaseFlowGDGD<wendlandc6, linear>(m, "SinglePFGDGD_WC6_L");
    export_SinglePhaseFlowGDGD<wendlandc6, tait>  (m, "SinglePFGDGD_WC6_T");
    export_SinglePhaseFlowGDGD<quintic,    linear>(m, "SinglePFGDGD_Q_L");
    export_SinglePhaseFlowGDGD<quintic,    tait>  (m, "SinglePFGDGD_Q_T");
    export_SinglePhaseFlowGDGD<cubicspline,linear>(m, "SinglePFGDGD_CS_L");
    export_SinglePhaseFlowGDGD<cubicspline,tait>  (m, "SinglePFGDGD_CS_T");

    export_SinglePhaseFlowFS<wendlandc2, linear>(m, "SinglePFFS_WC2_L");
    export_SinglePhaseFlowFS<wendlandc2, tait>  (m, "SinglePFFS_WC2_T");
    export_SinglePhaseFlowFS<wendlandc4, linear>(m, "SinglePFFS_WC4_L");
    export_SinglePhaseFlowFS<wendlandc4, tait>  (m, "SinglePFFS_WC4_T");
    export_SinglePhaseFlowFS<wendlandc6, linear>(m, "SinglePFFS_WC6_L");
    export_SinglePhaseFlowFS<wendlandc6, tait>  (m, "SinglePFFS_WC6_T");
    export_SinglePhaseFlowFS<quintic,    linear>(m, "SinglePFFS_Q_L");
    export_SinglePhaseFlowFS<quintic,    tait>  (m, "SinglePFFS_Q_T");
    export_SinglePhaseFlowFS<cubicspline,linear>(m, "SinglePFFS_CS_L");
    export_SinglePhaseFlowFS<cubicspline,tait>  (m, "SinglePFFS_CS_T");

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

    export_TwoPhaseFlowTV<wendlandc2, linear, linear>(m, "TwoPFTV_WC2_LL");
    export_TwoPhaseFlowTV<wendlandc2, linear, tait  >(m, "TwoPFTV_WC2_LT");
    export_TwoPhaseFlowTV<wendlandc2, tait,   linear>(m, "TwoPFTV_WC2_TL");
    export_TwoPhaseFlowTV<wendlandc2, tait,   tait  >(m, "TwoPFTV_WC2_TT");
    export_TwoPhaseFlowTV<wendlandc4, linear, linear>(m, "TwoPFTV_WC4_LL");
    export_TwoPhaseFlowTV<wendlandc4, linear, tait  >(m, "TwoPFTV_WC4_LT");
    export_TwoPhaseFlowTV<wendlandc4, tait,   linear>(m, "TwoPFTV_WC4_TL");
    export_TwoPhaseFlowTV<wendlandc4, tait,   tait  >(m, "TwoPFTV_WC4_TT");
    export_TwoPhaseFlowTV<wendlandc6, linear, linear>(m, "TwoPFTV_WC6_LL");
    export_TwoPhaseFlowTV<wendlandc6, linear, tait  >(m, "TwoPFTV_WC6_LT");
    export_TwoPhaseFlowTV<wendlandc6, tait,   linear>(m, "TwoPFTV_WC6_TL");
    export_TwoPhaseFlowTV<wendlandc6, tait,   tait  >(m, "TwoPFTV_WC6_TT");
    export_TwoPhaseFlowTV<quintic,    linear, linear>(m, "TwoPFTV_Q_LL");
    export_TwoPhaseFlowTV<quintic,    linear, tait  >(m, "TwoPFTV_Q_LT");
    export_TwoPhaseFlowTV<quintic,    tait,   linear>(m, "TwoPFTV_Q_TL");
    export_TwoPhaseFlowTV<quintic,    tait,   tait  >(m, "TwoPFTV_Q_TT");
    export_TwoPhaseFlowTV<cubicspline, linear, linear>(m, "TwoPFTV_CS_LL");
    export_TwoPhaseFlowTV<cubicspline, linear, tait  >(m, "TwoPFTV_CS_LT");
    export_TwoPhaseFlowTV<cubicspline, tait,   linear>(m, "TwoPFTV_CS_TL");
    export_TwoPhaseFlowTV<cubicspline, tait,   tait  >(m, "TwoPFTV_CS_TT");

    export_CustomForceCompute(m);

    export_ComputeSPFMechanicalProperties(m);
    export_ComputeSolidProperties(m);

#ifdef ENABLE_HIP
    export_ComputeSPFMechanicalPropertiesGPU(m);
    export_VelocityVerletGPU(m);
#endif

    export_DensityMethod(m);
    export_ViscosityMethod(m);
    export_ColorGradientMethod(m);

    // export_LocalNeighborListDataHost(m);
    export_HalfStepHook(m);

#ifdef ENABLE_HIP
    export_SinglePhaseFlowGPU<wendlandc2, linear>(m, "SinglePFGPU_WC2_L");
    export_SinglePhaseFlowGPU<wendlandc2, tait>  (m, "SinglePFGPU_WC2_T");
    export_SinglePhaseFlowGPU<wendlandc4, linear>(m, "SinglePFGPU_WC4_L");
    export_SinglePhaseFlowGPU<wendlandc4, tait>  (m, "SinglePFGPU_WC4_T");
    export_SinglePhaseFlowGPU<wendlandc6, linear>(m, "SinglePFGPU_WC6_L");
    export_SinglePhaseFlowGPU<wendlandc6, tait>  (m, "SinglePFGPU_WC6_T");
    export_SinglePhaseFlowGPU<quintic,    linear>(m, "SinglePFGPU_Q_L");
    export_SinglePhaseFlowGPU<quintic,    tait>  (m, "SinglePFGPU_Q_T");
    export_SinglePhaseFlowGPU<cubicspline,linear>(m, "SinglePFGPU_CS_L");
    export_SinglePhaseFlowGPU<cubicspline,tait>  (m, "SinglePFGPU_CS_T");

    export_SinglePhaseFlowTVGPU<wendlandc2, linear>(m, "SinglePFTVGPU_WC2_L");
    export_SinglePhaseFlowTVGPU<wendlandc2, tait>  (m, "SinglePFTVGPU_WC2_T");
    export_SinglePhaseFlowTVGPU<wendlandc4, linear>(m, "SinglePFTVGPU_WC4_L");
    export_SinglePhaseFlowTVGPU<wendlandc4, tait>  (m, "SinglePFTVGPU_WC4_T");
    export_SinglePhaseFlowTVGPU<wendlandc6, linear>(m, "SinglePFTVGPU_WC6_L");
    export_SinglePhaseFlowTVGPU<wendlandc6, tait>  (m, "SinglePFTVGPU_WC6_T");
    export_SinglePhaseFlowTVGPU<quintic,    linear>(m, "SinglePFTVGPU_Q_L");
    export_SinglePhaseFlowTVGPU<quintic,    tait>  (m, "SinglePFTVGPU_Q_T");
    export_SinglePhaseFlowTVGPU<cubicspline,linear>(m, "SinglePFTVGPU_CS_L");
    export_SinglePhaseFlowTVGPU<cubicspline,tait>  (m, "SinglePFTVGPU_CS_T");

    export_SinglePhaseFlowGDGDGPU<wendlandc2, linear>(m, "SinglePFGDGDGPU_WC2_L");
    export_SinglePhaseFlowGDGDGPU<wendlandc2, tait>  (m, "SinglePFGDGDGPU_WC2_T");
    export_SinglePhaseFlowGDGDGPU<wendlandc4, linear>(m, "SinglePFGDGDGPU_WC4_L");
    export_SinglePhaseFlowGDGDGPU<wendlandc4, tait>  (m, "SinglePFGDGDGPU_WC4_T");
    export_SinglePhaseFlowGDGDGPU<wendlandc6, linear>(m, "SinglePFGDGDGPU_WC6_L");
    export_SinglePhaseFlowGDGDGPU<wendlandc6, tait>  (m, "SinglePFGDGDGPU_WC6_T");
    export_SinglePhaseFlowGDGDGPU<quintic,    linear>(m, "SinglePFGDGDGPU_Q_L");
    export_SinglePhaseFlowGDGDGPU<quintic,    tait>  (m, "SinglePFGDGDGPU_Q_T");
    export_SinglePhaseFlowGDGDGPU<cubicspline,linear>(m, "SinglePFGDGDGPU_CS_L");
    export_SinglePhaseFlowGDGDGPU<cubicspline,tait>  (m, "SinglePFGDGDGPU_CS_T");

    export_SinglePhaseFlowFSGPU<wendlandc2, linear>(m, "SinglePFFSGPU_WC2_L");
    export_SinglePhaseFlowFSGPU<wendlandc2, tait>  (m, "SinglePFFSGPU_WC2_T");
    export_SinglePhaseFlowFSGPU<wendlandc4, linear>(m, "SinglePFFSGPU_WC4_L");
    export_SinglePhaseFlowFSGPU<wendlandc4, tait>  (m, "SinglePFFSGPU_WC4_T");
    export_SinglePhaseFlowFSGPU<wendlandc6, linear>(m, "SinglePFFSGPU_WC6_L");
    export_SinglePhaseFlowFSGPU<wendlandc6, tait>  (m, "SinglePFFSGPU_WC6_T");
    export_SinglePhaseFlowFSGPU<quintic,    linear>(m, "SinglePFFSGPU_Q_L");
    export_SinglePhaseFlowFSGPU<quintic,    tait>  (m, "SinglePFFSGPU_Q_T");
    export_SinglePhaseFlowFSGPU<cubicspline,linear>(m, "SinglePFFSGPU_CS_L");
    export_SinglePhaseFlowFSGPU<cubicspline,tait>  (m, "SinglePFFSGPU_CS_T");

    export_TwoPhaseFlowGPU<wendlandc2, tait,   tait  >(m, "TwoPFGPU_WC2_TT");
    export_TwoPhaseFlowGPU<wendlandc2, tait,   linear>(m, "TwoPFGPU_WC2_TL");
    export_TwoPhaseFlowGPU<wendlandc2, linear, tait  >(m, "TwoPFGPU_WC2_LT");
    export_TwoPhaseFlowGPU<wendlandc2, linear, linear>(m, "TwoPFGPU_WC2_LL");
    export_TwoPhaseFlowGPU<wendlandc4, tait,   tait  >(m, "TwoPFGPU_WC4_TT");
    export_TwoPhaseFlowGPU<wendlandc4, tait,   linear>(m, "TwoPFGPU_WC4_TL");
    export_TwoPhaseFlowGPU<wendlandc4, linear, tait  >(m, "TwoPFGPU_WC4_LT");
    export_TwoPhaseFlowGPU<wendlandc4, linear, linear>(m, "TwoPFGPU_WC4_LL");
    export_TwoPhaseFlowGPU<wendlandc6, tait,   tait  >(m, "TwoPFGPU_WC6_TT");
    export_TwoPhaseFlowGPU<wendlandc6, tait,   linear>(m, "TwoPFGPU_WC6_TL");
    export_TwoPhaseFlowGPU<wendlandc6, linear, tait  >(m, "TwoPFGPU_WC6_LT");
    export_TwoPhaseFlowGPU<wendlandc6, linear, linear>(m, "TwoPFGPU_WC6_LL");
    export_TwoPhaseFlowGPU<quintic,    tait,   tait  >(m, "TwoPFGPU_Q_TT");
    export_TwoPhaseFlowGPU<quintic,    tait,   linear>(m, "TwoPFGPU_Q_TL");
    export_TwoPhaseFlowGPU<quintic,    linear, tait  >(m, "TwoPFGPU_Q_LT");
    export_TwoPhaseFlowGPU<quintic,    linear, linear>(m, "TwoPFGPU_Q_LL");
    export_TwoPhaseFlowGPU<cubicspline,tait,   tait  >(m, "TwoPFGPU_CS_TT");
    export_TwoPhaseFlowGPU<cubicspline,tait,   linear>(m, "TwoPFGPU_CS_TL");
    export_TwoPhaseFlowGPU<cubicspline,linear, tait  >(m, "TwoPFGPU_CS_LT");
    export_TwoPhaseFlowGPU<cubicspline,linear, linear>(m, "TwoPFGPU_CS_LL");

    export_TwoPhaseFlowTVGPU<wendlandc2, tait,   tait  >(m, "TwoPFTVGPU_WC2_TT");
    export_TwoPhaseFlowTVGPU<wendlandc2, tait,   linear>(m, "TwoPFTVGPU_WC2_TL");
    export_TwoPhaseFlowTVGPU<wendlandc2, linear, tait  >(m, "TwoPFTVGPU_WC2_LT");
    export_TwoPhaseFlowTVGPU<wendlandc2, linear, linear>(m, "TwoPFTVGPU_WC2_LL");
    export_TwoPhaseFlowTVGPU<wendlandc4, tait,   tait  >(m, "TwoPFTVGPU_WC4_TT");
    export_TwoPhaseFlowTVGPU<wendlandc4, tait,   linear>(m, "TwoPFTVGPU_WC4_TL");
    export_TwoPhaseFlowTVGPU<wendlandc4, linear, tait  >(m, "TwoPFTVGPU_WC4_LT");
    export_TwoPhaseFlowTVGPU<wendlandc4, linear, linear>(m, "TwoPFTVGPU_WC4_LL");
    export_TwoPhaseFlowTVGPU<wendlandc6, tait,   tait  >(m, "TwoPFTVGPU_WC6_TT");
    export_TwoPhaseFlowTVGPU<wendlandc6, tait,   linear>(m, "TwoPFTVGPU_WC6_TL");
    export_TwoPhaseFlowTVGPU<wendlandc6, linear, tait  >(m, "TwoPFTVGPU_WC6_LT");
    export_TwoPhaseFlowTVGPU<wendlandc6, linear, linear>(m, "TwoPFTVGPU_WC6_LL");
    export_TwoPhaseFlowTVGPU<quintic,    tait,   tait  >(m, "TwoPFTVGPU_Q_TT");
    export_TwoPhaseFlowTVGPU<quintic,    tait,   linear>(m, "TwoPFTVGPU_Q_TL");
    export_TwoPhaseFlowTVGPU<quintic,    linear, tait  >(m, "TwoPFTVGPU_Q_LT");
    export_TwoPhaseFlowTVGPU<quintic,    linear, linear>(m, "TwoPFTVGPU_Q_LL");
    export_TwoPhaseFlowTVGPU<cubicspline,tait,   tait  >(m, "TwoPFTVGPU_CS_TT");
    export_TwoPhaseFlowTVGPU<cubicspline,tait,   linear>(m, "TwoPFTVGPU_CS_TL");
    export_TwoPhaseFlowTVGPU<cubicspline,linear, tait  >(m, "TwoPFTVGPU_CS_LT");
    export_TwoPhaseFlowTVGPU<cubicspline,linear, linear>(m, "TwoPFTVGPU_CS_LL");
#endif

}

