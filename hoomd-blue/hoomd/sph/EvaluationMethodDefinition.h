/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/
#ifndef __EVALUATION_METHOD_DEFINITON_H__
#define __EVALUATION_METHOD_DEFINITON_H__

#include "hoomd/HOOMDMath.h"
/*! \file EvaluationMethodDefinition.h
    \brief Data structures to define Evaluation methods for density, viscosity etc
    */

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

    //! Enum for various viscosity evaluation approaches
    enum MaterialModel
    {
        REGULARIZEDBINGHAM, //!< Regularized Bingham material model
        BIVISCOUS, //!< Biviscous material model
    };


    } // end namespace sph
    } // end namespace hoomd

#endif