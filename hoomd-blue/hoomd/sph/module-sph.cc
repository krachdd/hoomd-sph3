/* ---------------------------------------------------------
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
#include "SuspendedObjectIntegrator.h"
// // #include "RigidBodyIntegrator.h"
#include "SinglePhaseFlow.h"
#include "SinglePhaseFlowNN.h"
#include "SinglePhaseFlowTV.h"
#include "SinglePhaseFlowTVNN.h"
#include "SuspensionFlow.h"
#include "SuspensionFlowNN.h"
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
    export_SuspendedObjectIntegrator(m);
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

    export_SinglePhaseFlowNN<wendlandc2, linear>(m, "SinglePFNN_WC2_L");
    export_SinglePhaseFlowNN<wendlandc2, tait>(m, "SinglePFNN_WC2_T");
    export_SinglePhaseFlowNN<wendlandc4, linear>(m, "SinglePFNN_WC4_L");
    export_SinglePhaseFlowNN<wendlandc4, tait>(m, "SinglePFNN_WC4_T");
    export_SinglePhaseFlowNN<wendlandc6, linear>(m, "SinglePFNN_WC6_L");
    export_SinglePhaseFlowNN<wendlandc6, tait>(m, "SinglePFNN_WC6_T");
    export_SinglePhaseFlowNN<quintic, linear>(m, "SinglePFNN_Q_L");
    export_SinglePhaseFlowNN<quintic, tait>(m, "SinglePFNN_Q_T");
    export_SinglePhaseFlowNN<cubicspline, linear>(m, "SinglePFNN_CS_L");
    export_SinglePhaseFlowNN<cubicspline, tait>(m, "SinglePFNN_CS_T");

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

    export_SinglePhaseFlowTVNN<wendlandc2, linear>(m, "SinglePFTVNN_WC2_L");
    export_SinglePhaseFlowTVNN<wendlandc2, tait>(m, "SinglePFTVNN_WC2_T");
    export_SinglePhaseFlowTVNN<wendlandc4, linear>(m, "SinglePFTVNN_WC4_L");
    export_SinglePhaseFlowTVNN<wendlandc4, tait>(m, "SinglePFTVNN_WC4_T");
    export_SinglePhaseFlowTVNN<wendlandc6, linear>(m, "SinglePFTVNN_WC6_L");
    export_SinglePhaseFlowTVNN<wendlandc6, tait>(m, "SinglePFTVNN_WC6_T");
    export_SinglePhaseFlowTVNN<quintic, linear>(m, "SinglePFTVNN_Q_L");
    export_SinglePhaseFlowTVNN<quintic, tait>(m, "SinglePFTVNN_Q_T");
    export_SinglePhaseFlowTVNN<cubicspline, linear>(m, "SinglePFTVNN_CS_L");
    export_SinglePhaseFlowTVNN<cubicspline, tait>(m, "SinglePFTVNN_CS_T");

    export_SuspensionFlow<wendlandc2, linear>(m, "SuspensionF_WC2_L");
    export_SuspensionFlow<wendlandc2, tait>(m, "SuspensionF_WC2_T");
    export_SuspensionFlow<wendlandc4, linear>(m, "SuspensionF_WC4_L");
    export_SuspensionFlow<wendlandc4, tait>(m, "SuspensionF_WC4_T");
    export_SuspensionFlow<wendlandc6, linear>(m, "SuspensionF_WC6_L");
    export_SuspensionFlow<wendlandc6, tait>(m, "SuspensionF_WC6_T");
    export_SuspensionFlow<quintic, linear>(m, "SuspensionF_Q_L");
    export_SuspensionFlow<quintic, tait>(m, "SuspensionF_Q_T");
    export_SuspensionFlow<cubicspline, linear>(m, "SuspensionF_CS_L");
    export_SuspensionFlow<cubicspline, tait>(m, "SuspensionF_CS_T");

    export_SuspensionFlowNN<wendlandc2, linear>(m, "SuspensionFNN_WC2_L");
    export_SuspensionFlowNN<wendlandc2, tait>(m, "SuspensionFNN_WC2_T");
    export_SuspensionFlowNN<wendlandc4, linear>(m, "SuspensionFNN_WC4_L");
    export_SuspensionFlowNN<wendlandc4, tait>(m, "SuspensionFNN_WC4_T");
    export_SuspensionFlowNN<wendlandc6, linear>(m, "SuspensionFNN_WC6_L");
    export_SuspensionFlowNN<wendlandc6, tait>(m, "SuspensionFNN_WC6_T");
    export_SuspensionFlowNN<quintic, linear>(m, "SuspensionFNN_Q_L");
    export_SuspensionFlowNN<quintic, tait>(m, "SuspensionFNN_Q_T");
    export_SuspensionFlowNN<cubicspline, linear>(m, "SuspensionFNN_CS_L");
    export_SuspensionFlowNN<cubicspline, tait>(m, "SuspensionFNN_CS_T");

    export_CustomForceCompute(m);

    export_ComputeSPFMechanicalProperties(m);
    export_ComputeSolidProperties(m);

    export_DensityMethod(m);
    export_ViscosityMethod(m);
    export_MaterialModel(m);

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

