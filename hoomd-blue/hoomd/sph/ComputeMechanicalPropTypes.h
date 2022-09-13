/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/
#ifndef _COMPUTE_MECHANICALPROP_TYPES_H_
#define _COMPUTE_MECHANICALPROP_TYPES_H_

#include "hoomd/HOOMDMath.h"
/*! \file ComputeMechanicalPropTypes.h
    \brief Data structures common to both CPU and GPU implementations of ComputeMechanicalProp
    */

namespace hoomd
    {
namespace sph
    {
//! Enum for indexing the GPUArray of computed values
struct singlephaseflow_logger_index
    {
    //! The enum
    enum Enum
        {
        sum_fluid_velocity_x=0, //!< Index for the sum of fluid x-velocity in the GPUArray
        sum_fluid_velocity_y,   //!< Index for the sum of fluid y-velocity in the GPUArray
        sum_fluid_velocity_z,   //!< Index for the sum of fluid z-velocity in the GPUArray
        kinetic_energy,         //!< Index for the overall kinetic energy of the system
        // total_fluid_particles,  //!< Total number of fluid particles
        // dt_adapt,               //!< Adaptive timestep size
        num_quantities // final element to count number of quantities
        };
    };

    } // end namespace sph
    } // end namespace hoomd

#endif