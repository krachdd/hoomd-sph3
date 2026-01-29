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

#include "EvaluationMethodDefinition.h"

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
// //! Enum for various density evaluation approaches
// enum DensityMethod
// {
//     DENSITYSUMMATION,    //!< Summation approach
//     DENSITYCONTINUITY,    //!< Continuity approach
// };

// //! Enum for various viscosity evaluation approaches
// enum ViscosityMethod
// {
//     HARMONICAVERAGE, //!< Viscosity operator based on inter-particle averaged shear stress
// };

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

        Index2D m_typpair_idx;        //!< Helper class for indexing per type pair arrays

        DensityMethod m_densitymethod;
        ViscosityMethod m_viscositymethod;
        ColorGradientMethod m_colorgradient_method;

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

void export_ColorGradientMethod(pybind11::module& m);

} // end namespace detail

} // end namespace sph
} // end namespace hoomd

#endif // __SPHBaseClass_H__
