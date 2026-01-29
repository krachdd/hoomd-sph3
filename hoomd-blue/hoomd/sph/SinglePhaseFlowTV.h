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


/*! \file SinglePhaseFlowTV.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow using Adami type transport velocity
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#ifndef __SinglePhaseFlowTV_H__
#define __SinglePhaseFlowTV_H__


namespace hoomd 
{
namespace sph
{

//! Computes SinglePhaseFlowTV forces on each particle
/*!
*/
template<SmoothingKernelType KT_,StateEquationType SET_>
class PYBIND11_EXPORT SinglePhaseFlowTV : public SinglePhaseFlow<KT_, SET_>
    {
    public:

        //! Constructor
        SinglePhaseFlowTV(std::shared_ptr<SystemDefinition> sysdef,
                        std::shared_ptr<SmoothingKernel<KT_> > skernel,
                        std::shared_ptr<StateEquation<SET_> > equationofstate,
                        std::shared_ptr<nsearch::NeighborList> nlist,
                        std::shared_ptr<ParticleGroup> fluidgroup,
                        std::shared_ptr<ParticleGroup> solidgroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~SinglePhaseFlowTV();

        //! Computes forces
        void computeForces(uint64_t timestep);

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
            flags[comm_flag::auxiliary2] = 1; // Stores background pressure contribution to tv
            flags[comm_flag::auxiliary3] = 1; // Stores transport velocity
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

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputation(uint64_t timestep);

        /*! Helper function to set communication flags and update ghosts auxiliary array 1
        * \param timestep The time step
        * \post Ghost particle auxiliary array 1, 2 and 3 is up-to-date
        */
        void update_ghost_aux123(uint64_t timestep);

    private:

    };


namespace detail 
{
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowTV(pybind11::module& m, std::string name);

} // end namespace detail
} // end namespace sph
} // end namespace hoomd

#endif // __SinglePhaseFlowTV_H__