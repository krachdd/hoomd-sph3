/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/
#ifndef __COMPUTE_MECHANICALPROP_TYPES_H__
#define __COMPUTE_MECHANICALPROP_TYPES_H__

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
        sum_fluid_density,      //!< Index for the mean fluid particle density
        abs_velocity,           //!< Index for the overall absolute velocity per particle
        e_kin_fluid,            //!< Index for the overall kinetic energy of the system
        num_quantities // final element to count number of quantities
        };
    };

struct solidphase_logger_index
    {
    //! The enum
    enum Enum
        {
        total_drag_x,           //!< Index for the overall total drag on the solid phase
        total_drag_y,           //!< Index for the overall total drag on the solid phase
        total_drag_z,           //!< Index for the overall total drag on the solid phase
        num_quantities // final element to count number of quantities
        };
    };

struct suspendedphase_logger_index
    {
    //! The enum
    enum Enum
        {
        total_drag_x,           //!< Index for the overall total drag on the solid phase
        total_drag_y,           //!< Index for the overall total drag on the solid phase
        total_drag_z,           //!< Index for the overall total drag on the solid phase
        num_quantities // final element to count number of quantities
        };
    };


    } // end namespace sph
    } // end namespace hoomd

#endif