/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/Integrator.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

#include <memory>
#include <vector>
#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "SPHBaseClass.h"
#include "SinglePhaseFlow.h"
#include "SolidFluidTypeBit.h"

#include "EvaluationMethodDefinition.h"


/*! \file SinglePhaseFlowNN.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowNN_H__
#define __SinglePhaseFlowNN_H__


namespace hoomd 
{
namespace sph
{

//! Computes SinglePhaseFlowNN forces on each particle
/*!
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowNN : public SinglePhaseFlow<KT_, SET_>
    {
    public:

        //! Constructor
        SinglePhaseFlowNN(std::shared_ptr<SystemDefinition> sysdef,
                        std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,
                        std::shared_ptr<nsearch::NeighborList> nlist,
                        std::shared_ptr<ParticleGroup> fluidgroup,
                        std::shared_ptr<ParticleGroup> solidgroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE,
                        MaterialModel mmaterialmodel=REGULARIZEDBINGHAM);

        //! Destructor
        virtual ~SinglePhaseFlowNN();


        /*! Set the parameters
         * \param k, yield stress tau, regularization parameter m
         */
        virtual void setParams(Scalar k, Scalar tau, Scalar m);

        //! Getter and Setter methods for density method
        MaterialModel getMaterialModel()
            {
            return m_material_model;
            }
        void setMaterialModel(MaterialModel materialmodel)
            {
            m_material_model = materialmodel;
            }

        //! Computes forces
        virtual void computeForces(uint64_t timestep);

    #ifdef ENABLE_MPI
        /// The system's communicator.
        std::shared_ptr<Communicator> m_comm;
    #endif

    #ifdef ENABLE_MPI
        //! Get requested ghost communication flags
        virtual CommFlags getRequestedCommFlags(uint64_t timestep)
            {
            // Request communication of all field required during ForceCompute
            CommFlags flags(0);
            flags[comm_flag::net_force] = 0;
            flags[comm_flag::position] = 1; // Stores position and type
            flags[comm_flag::velocity] = 1; // Stores velocity and mass
            flags[comm_flag::density] = 1; // Stores density 
            flags[comm_flag::pressure] = 1; // Stores pressure
            flags[comm_flag::energy] = 0; // Stores density and pressure
            flags[comm_flag::auxiliary1] = 1; // Stores fictitious velocity
            flags[comm_flag::auxiliary4] = 1; // Stores viscosity/shear stress/shear rate
            flags[comm_flag::slength] = 1; // Stores smoothing length TODO is this needed
            // Add flags requested by base class
            flags |= ForceCompute::getRequestedCommFlags(timestep);
            return flags;
            }
    #endif

        //! Returns true because we compute dpe array content
        virtual bool ComputesDPE()
            {
            return true;
            }

    protected:

        // Model parameters
        MaterialModel m_material_model; //!< Material model to use


        // Physical variables
        Scalar m_k; //!< Viscosity ( Must be set by user )
        Scalar m_tau; //!< Shear stress ( Must be set by user )
        Scalar m_m; //!< m value ( Must be set by user )


        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        virtual void forcecomputation(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 4
        * \param timestep The time step
        * \post Ghost particle auxiliary array 4 is up-to-date
        */
        void update_ghost_aux4(uint64_t timestep);

        /*! Helper function to compute particle viscosites
         *  \post Viscosity of fluid particle computed
         */
        void compute_viscosity(uint64_t timestep);
        
    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowNN(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SinglePhaseFlowNN_H__
