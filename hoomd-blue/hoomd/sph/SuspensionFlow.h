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


/*! \file SuspensionFlow.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SuspensionFlow_H__
#define __SuspensionFlow_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_TWOPI
#define M_TWOPI 6.28318530717958647692
#endif

namespace hoomd 
{
namespace sph
{

//! Computes SuspensionFlow forces on each particle
/*!
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
class PYBIND11_EXPORT SuspensionFlow : public SinglePhaseFlow<KT_, SET_>
    {
    public:

        //! Constructor
        SuspensionFlow(std::shared_ptr<SystemDefinition> sysdef,
                        std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,
                        std::shared_ptr<nsearch::NeighborList> nlist,
                        std::shared_ptr<ParticleGroup> fluidgroup,
                        std::shared_ptr<ParticleGroup> solidgroup,
                        std::shared_ptr<ParticleGroup> suspendedgroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~SuspensionFlow();

        /*! Set the parameters
         * \param mu Dynamic viscosity, solid density, repulsive force factor
         */
        virtual void setParams(Scalar mu, Scalar rhoS, Scalar f0);

        /*! Set compute solid forces option to true. This is necessary if suspended object
         *  are present or if solid drag forces are to be evaluated.
         */
        void deactivateFluidAcceleration()
            {
            m_compute_fluid_acceleration = false;
            }

        /*! Set compute solid forces option to true. This is necessary if suspended object
         *  are present or if solid drag forces are to be evaluated.
         */
        void computeSolidForces()
            {
            m_compute_solid_forces = true;
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
            flags[comm_flag::auxiliary2] = 1; // Stores fictitious repulsive force
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

        // Shared pointers
        // std::shared_ptr<ParticleGroup> m_fluidgroup; //!< Group of fluid particles
        // std::shared_ptr<ParticleGroup> m_solidgroup; //!< Group of fluid particles
        std::shared_ptr<ParticleGroup> m_suspendedgroup; //!< Group of fluid particles

        // Auxiliary variables
        // std::vector<unsigned int> m_fluidtypes; //!< Fluid type numbers
        // std::vector<unsigned int> m_solidtypes; //!< Solid type numbers
        std::vector<unsigned int> m_suspendedtypes; //!< Solid type numbers
        // GPUArray<unsigned int> m_type_property_map; //!< to check if a particle type is solid or fluid
        unsigned int m_maxSuspendedID;              //!< maximum nuumber of solid bodies        
        std::vector<Scalar4> m_centerofmasses;  //!< Vector with center of mass and typeID for each solid body besides walls  
        std::vector<Scalar3> m_repulsiveforces; //!< Vector with repulsive forces per solid body     
        std::vector<Scalar3> m_velocities;      //!< Vector with velocities of solid bodies 
        Scalar m_f0; //!< Spring stiffness of contact force
        Scalar m_rhoS; //!< Solid density
        std::vector<Scalar> m_radii;            //!< Vector with equivalent radi of all solid bodies

        // Flags
        bool m_compute_fluid_acceleration;
        bool m_compute_solid_forces;

        // Log parameters
        uint64_t m_log_computed_last_timestep; //!< Last time step where log quantities were computed


        // /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        // * \pre Ghost particle number densities (i.e. density array) must be up-to-date
        // * \pre Solid normalization constant \sum_j w_ij must be computed and stored in density array
        // * \post Fictitious particle properties are computed and stored in aux1 array
        // */
        // void compute_noslipsolid(uint64_t timestep);

        /*! Helper function to compute fictitious suspended particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. density array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in density array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
        void compute_noslipsuspended(uint64_t timestep);

        /*! Helper function to calculate center of mass of suspendedtypes
        * \post Based on SuspendedObjectIntegrator
        */
        void compute_Centerofmasses(uint64_t timestep);

        /*! Helper function to calculate equivalent radii of solidtypes
        * \post eqivalent radii stored in m_radi dependent on type id
        */
        void compute_equivalentRadii(uint64_t timestep);

        /*! Helper function to calculate repulsive force between solids
        * \post Based on approach of Bian, Ellero, Vazquez
        */
        void compute_repulsiveForce(uint64_t timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        virtual void forcecomputation(uint64_t timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputationsuspended(uint64_t timestep);

        // /*! Helper function to set communication flags and update ghosts densities
        // * \param timestep The time step
        // * \post Ghost particle density array is up-to-date
        // */
        // void update_ghost_density(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 1
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1 is up-to-date
        */
        void update_ghost_aux2(uint64_t timestep);

    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SuspensionFlow(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SuspensionFlow_H__