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


#ifndef __EVALUATION_METHOD_DEFINITON_H__
#define __EVALUATION_METHOD_DEFINITON_H__

#include "hoomd/HOOMDMath.h"
#include <cmath>
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

    //! Enum for various Colorgradient evaluation approaches
    enum ColorGradientMethod
    {
        DENSITYRATIO, //!< Method to compute TPF Color gradient
        NUMBERDENSITY, //!< Method to compute TPF Color gradient
    };

    //! Enum for non-Newtonian rheology models
    enum NonNewtonianModel
    {
        NEWTONIAN,          //!< Constant viscosity (default)
        POWERLAW,           //!< Power-law (Ostwald-de Waele): mu_eff = max(mu_min, K*|gdot|^(n-1))
        CARREAU,            //!< Carreau: mu_eff = muinf + (mu0-muinf)*(1+(lambda*|gdot|)^2)^((n-1)/2)
        BINGHAM,            //!< Bingham (Papanastasiou): mu_eff = max(mu_min, mu_p + tauy*(1-exp(-m|gdot|))/|gdot|)
        HERSCHELBULKLEY,    //!< Herschel-Bulkley (Papanastasiou): mu_eff = max(mu_min, K*|gdot|^(n-1) + tauy*(1-exp(-m|gdot|))/|gdot|)
    };

    /*! Compute effective dynamic viscosity for non-Newtonian rheology models.
     *
     *  \param mu_base   Base Newtonian dynamic viscosity [Pa·s] (used for NEWTONIAN only)
     *  \param gamma_dot Shear-rate estimate |v_ij|/r_ij  [1/s], >= 0
     *  \param model     Active NonNewtonianModel
     *  \param K         Power-law / H-B consistency index [Pa·s^n]; or plastic viscosity mu_p for BINGHAM
     *  \param n         Power-law / Carreau / H-B exponent
     *  \param mu0       Carreau zero-shear viscosity      [Pa·s]
     *  \param muinf     Carreau infinite-shear viscosity  [Pa·s]
     *  \param lambda_NN Carreau relaxation time           [s]
     *  \param tauy      Yield stress (BINGHAM / H-B)      [Pa]
     *  \param m_reg     Papanastasiou regularization param [s]
     *  \param mu_min    Lower viscosity clamp             [Pa·s]
     *  \returns         Effective dynamic viscosity [Pa·s]
     */
    inline Scalar computeNNViscosity(
        Scalar mu_base,
        Scalar gamma_dot,
        NonNewtonianModel model,
        Scalar K,
        Scalar n,
        Scalar mu0,
        Scalar muinf,
        Scalar lambda_NN,
        Scalar tauy,
        Scalar m_reg,
        Scalar mu_min)
        {
        switch (model)
            {
            case NEWTONIAN:
                return mu_base;

            case POWERLAW:
                {
                // Clamp to avoid 0^(n-1) = infinity for shear-thinning fluids (n < 1)
                Scalar gdot = gamma_dot > Scalar(1e-12) ? gamma_dot : Scalar(1e-12);
                Scalar mu_eff = K * std::pow(gdot, n - Scalar(1.0));
                return mu_eff > mu_min ? mu_eff : mu_min;
                }

            case CARREAU:
                {
                Scalar lg = lambda_NN * gamma_dot;
                Scalar mu_eff = muinf + (mu0 - muinf)
                                * std::pow(Scalar(1.0) + lg * lg,
                                           Scalar(0.5) * (n - Scalar(1.0)));
                return mu_eff;
                }

            case BINGHAM:
                {
                // K stores the plastic viscosity mu_p
                // Papanastasiou regularization; L'Hopital limit at gdot -> 0:
                //   tauy*(1-exp(-m*gdot))/gdot -> tauy*m
                Scalar mu_eff;
                if (gamma_dot < Scalar(1e-12))
                    mu_eff = K + tauy * m_reg;
                else
                    mu_eff = K + tauy * (Scalar(1.0) - std::exp(-m_reg * gamma_dot)) / gamma_dot;
                return mu_eff > mu_min ? mu_eff : mu_min;
                }

            case HERSCHELBULKLEY:
                {
                Scalar gdot = gamma_dot > Scalar(1e-12) ? gamma_dot : Scalar(1e-12);
                Scalar mu_eff = K * std::pow(gdot, n - Scalar(1.0))
                                + tauy * (Scalar(1.0) - std::exp(-m_reg * gdot)) / gdot;
                return mu_eff > mu_min ? mu_eff : mu_min;
                }

            default:
                return mu_base;
            }
        }

    } // end namespace sph
    } // end namespace hoomd

#endif
