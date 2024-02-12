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
        abs_velocity,           //!< Index for the overall kinetic energy of the system
        num_quantities // final element to count number of quantities
        };
    };

//! Enum for indexing the GPUArray of computed values
struct singlephaseflownn_logger_index
    {
    //! The enum
    enum Enum
        {
        sum_fluid_velocity_x=0, //!< Index for the sum of fluid x-velocity in the GPUArray
        sum_fluid_velocity_y,   //!< Index for the sum of fluid y-velocity in the GPUArray
        sum_fluid_velocity_z,   //!< Index for the sum of fluid z-velocity in the GPUArray
        sum_fluid_density,      //!< Index for the mean fluid particle density
        abs_velocity,           //!< Index for the overall kinetic energy of the system
        mean_viscosity,         //!< Index for the mean viscosity in the system
        max_viscosity,          //!< Index for the max viscosity in the system
        max_shearrate,            //!< Index for the max shearrate in the system
        num_quantities // final element to count number of quantities
        };
    };

//! Enum for indexing the GPUArray of computed values
struct suspensionflow_logger_index
    {
    //! The enum
    enum Enum
        {
        sum_fluid_velocity_x=0, //!< Index for the sum of fluid x-velocity in the GPUArray
        sum_fluid_velocity_y,   //!< Index for the sum of fluid y-velocity in the GPUArray
        sum_fluid_velocity_z,   //!< Index for the sum of fluid z-velocity in the GPUArray
        sum_fluid_density,      //!< Index for the mean fluid particle density
        abs_velocity,           //!< Index for the overall kinetic energy of the system
        angularvel_x,
        angularvel_y,
        angularvel_z,
        translationvel_x,
        translationvel_y,
        translationvel_z,
        com_x,
        com_y,
        com_z,
        sum_force_x,
        sum_force_y,
        sum_force_z,
        num_quantities // final element to count number of quantities
        };
    };

    } // end namespace sph
    } // end namespace hoomd

#endif