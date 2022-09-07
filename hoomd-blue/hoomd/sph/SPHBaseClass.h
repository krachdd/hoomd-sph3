// Copyright (c) 2009-2022 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "hoomd/Compute.h"
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/nsearch/NeighborList.h"

#ifdef ENABLE_MPI
#include "hoomd/Communicator.h"
#include "hoomd/HOOMDMPI.h"
#endif

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#include "SmoothingKernel.h"
#include "StateEquations.h"

/*! \file SPHBaseClass.cc
    \brief Contains base class for any SPH Force compute. Takes care of
           storing SmoothingKernel and NeighborList class instances.
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef __SPHBaseClass_H__
#define __SPHBaseClass_H__

namespace hoomd
{
namespace sph
{
//! Enum for various density evaluation approaches
enum DensityMethod
{
    DENSITYSUMMATION,    //!< Summation approach
    DENSITYCONTINUITY,    //!< Continuity approach
};

//! Enum for various viscosity evaluation approaches
enum ViscosityMethod
{
    HARMONICAVERAGE, //!< Viscosity operator based on inter-particle averaged shear stress
};

template<SmoothingKernelType KT_, StateEquationType SET_>
class PYBIND11_EXPORT SPHBaseClass : public ForceCompute
    {
    public:
        
        //! Constructor
        SPHBaseClass(std::shared_ptr<SystemDefinition> sysdef,
                     std::shared_ptr<SmoothingKernel<KT_> > skernel,
                     std::shared_ptr<StateEquation<SET_> > eos,
                     std::shared_ptr<nsearch::NeighborList> nlist);

        //! Destructor
        virtual ~SPHBaseClass();

        /*! Helper function to compute available type ids for a given group of particles
         * \param pgroup Group of particles to construct type id vectors for
         */
        void constructTypeVectors(std::shared_ptr<ParticleGroup> const pgroup,
                                  std::vector<unsigned int> *global_typeids);

        /*! Helper function to apply external body force to a given group of particles
         * \param pgroup Group of particles to apply body force to
         */
        void applyBodyForce(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup);
        
// #ifdef ENABLE_HIP
//         void applyBodyForceGPU(uint64_t timestep, std::shared_ptr<ParticleGroup> pgroup);
// #endif

        /*! Set the volumetric acceleration
         * \param gx Volumetric acceleration in x-Direction
         * \param gy Volumetric acceleration in y-Direction
         * \param gz Volumetric acceleration in z-Direction
         * \param damp damping time in units of time steps during which body acceleration is smoothly applied
         */
        void setAcceleration(Scalar gx, Scalar gy, Scalar gz, unsigned int damptime);

        // Get the volumetric acceleration
        Scalar3 getAcceleration(uint64_t timestep);

        
// #ifdef ENABLE_MPI
//         /// The system's communicator.
//         std::shared_ptr<Communicator> m_comm;
//         //! Set the communicator to use
//         void setCommunicator(std::shared_ptr<Communicator> comm)
//             {
//             if (!m_comm)
//                 {
//                 assert(comm);
//                 }
//             Compute::setCommunicator(comm);
//             }
// #endif

    protected:
        std::shared_ptr<SmoothingKernel<KT_> > m_skernel; //!< The kernel function class this method is associated with
        std::shared_ptr<StateEquation<SET_> > m_eos; //!< The equation of state class this method is associated with
        std::shared_ptr<nsearch::NeighborList> m_nlist; //!< The neighbor list to use for the computation

        DensityMethod m_densitymethod;
        ViscosityMethod m_viscositymethod;

        Scalar3 m_bodyforce; //!< Volumetric force
        unsigned int m_damptime; //!< Damping time
        bool m_body_acceleration; //!< True if body acceleration has been set and not null
        


    // private:
    //     //! Connection to the signal notifying when number of particle types changes
    //     boost::signals2::connection m_num_type_change_connection;

    //     //! Connection to the signal notifying when number of particle types changes
    //     boost::signals2::connection m_particle_num_change_connection;
    };


namespace detail 
{

template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SPHBaseClass(pybind11::module& m, std::string name);

void export_DensityMethod(pybind11::module& m);

void export_ViscosityMethod(pybind11::module& m);

} // end namespace detail

} // end namespace sph
} // end namespace hoomd

#endif // __SPHBaseClass_H__
