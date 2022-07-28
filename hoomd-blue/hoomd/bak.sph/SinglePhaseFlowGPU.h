// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: schirwon

#include "SinglePhaseFlow.h"
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
#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/signals2.hpp>

#include "StateEquations.h"
#include "SPHBaseClass.h"
#include "SinglePhaseFlowGPU.cuh"
#include "SinglePhaseFlow.h"

/*! \file SinglePhaseFlow.h
    \brief Contains code for the Quasi-incompressible Navier-Stokes solver
          for Single-phase flow
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

//#ifndef __SinglePhaseFlow_H__
//#define __SinglePhaseFlow_H__

//! Computes SinglePhaseFlow forces on each particle
/*!
*/
template<SmoothingKernelType KT_, StateEquationType SET_>
class SinglePhaseFlowGPU : public SinglePhaseFlow<KT_, SET_>
    {
    public:

        //! Constructor
        SinglePhaseFlowGPU(boost::shared_ptr<SystemDefinition> sysdef,
                        boost::shared_ptr<SmoothingKernel<KT_> > skernel,
                        boost::shared_ptr<StateEquation<SET_> > equationofstate,
                        boost::shared_ptr<nsearch::NeighborList> nlist,
                        boost::shared_ptr<ParticleGroup> fluidgroup,
                        boost::shared_ptr<ParticleGroup> solidgroup,
                        DensityMethod   mdensitymethod=DENSITYSUMMATION,
                        ViscosityMethod mviscositymethod=HARMONICAVERAGE);

        //! Destructor
        virtual ~SinglePhaseFlowGPU();

        //! Computes forces
        void computeForces(unsigned int timestep);
        
    protected:

        //! Computes log properties
        void computeProperties();

#ifdef ENABLE_MPI_CUDA

        //! Reduce properties over MPI
        void reduceProperties();
#endif

        /*! Helper function to compute particle number density
         * \post For fluid particles, compute number density. For solid particles,
                 compute fluid normalization constant.
         */
        void compute_ndensity(unsigned int timestep);

        /*! Helper function to compute particle pressures
         *  \post Pressure of fluid particle computed
         */
        void compute_pressure(unsigned int timestep);

        /*! Helper function to compute fictitious solid particle properties (pressures and velocities)
        * \pre Ghost particle number densities (i.e. dpe array) must be up-to-date
        * \pre Solid normalization constant \sum_j w_ij must be computed and stored in dpe array
        * \post Fictitious particle properties are computed and stored in aux1 array
        */
         void compute_noslip(unsigned int timestep);

        /*! Helper function to apply Shepard density filter
        * \post Fluid particle densities are recomputed based on the Shepard renormalization
        */
         void renormalize_density(unsigned int timestep);

        /*! Helper function where the actual force computation takes place
         * \pre Number densities and fictitious solid particle properties must be up-to-date
         * \post h_force stores forces acting on fluid particles and .w component stores rate of change of density
         */
        void forcecomputation(unsigned int timestep);


    private:

    };

//! Exports the SinglePhaseFlow class to python
template<SmoothingKernelType KT_, StateEquationType SET_>
void export_SinglePhaseFlowGPU_templ();

void export_SinglePhaseFlowGPU();

//#endif // __SinglePhaseFlow_H__
