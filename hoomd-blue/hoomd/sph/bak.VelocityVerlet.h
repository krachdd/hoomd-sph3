// Copyright (c) 2009-2022 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "SPHIntegrationMethodTwoStep.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifndef __VELOCITY_VERLET_H__
#define __VELOCITY_VERLET_H__

/*! \file VelocityVerlet.h
    \brief Declares an Velocity-Verlet integrator for performing two-step integration on SPH particles
*/

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

namespace hoomd 
{
namespace sph
{
//! Integrates the system forward one step
/*! Implements velocity-verlet integration through the IntegrationMethodTwoStep interface
    \ingroup updaters
*/
class VelocityVerlet : public SPHIntegrationMethodTwoStep
    {
    public:
        //! Constructor
        VelocityVerlet(std::shared_ptr<SystemDefinition> sysdef,
                       std::shared_ptr<ParticleGroup> group,
                       bool skip_restart=false);

        //! Destructor
        virtual ~VelocityVerlet();

        //! Prepare for the run
        virtual void integrateStepOne(unsigned int timestep);

        //! Prepare for the run
        virtual void integrateStepTwo(unsigned int timestep);

    };

namespace detail {
//! Exports the VelocityVerlet class to python
void export_VelocityVerlet(pybind11::module& m);

}   // end namespace detail
}   // end namespace sph
}   // end namespace hoomd

#endif // #ifndef __VELOCITY_VERLET_H__
