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

/*! \file SPHDeviceFunctions.cuh
    \brief Device-side (GPU) math primitives for the SPH solver.

    This header provides __device__ __forceinline__ implementations of:
      - Smoothing-kernel wij / dwijdr / w0 / normfactor (templated by SmoothingKernelType)
      - Equation-of-state Pressure / dPressuredDensity / Density / VRD variants
        (templated by StateEquationType)
      - Non-Newtonian viscosity sph_nn_viscosity()
      - Type-check helpers sph_checksolid / sph_checkfluid
      - POD parameter structs for passing parameters to GPU kernels

    Must NOT be compiled by the host C++ compiler as a standalone header — it
    includes hip/hip_runtime.h and uses __device__ qualifiers.
    Include guards prevent double definition when the file is transitively
    included by both a .cu and a .cc translation unit.
*/

#ifndef __SPH_DEVICE_FUNCTIONS_CUH__
#define __SPH_DEVICE_FUNCTIONS_CUH__

#include "hip/hip_runtime.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"

// In device compilation (__HIPCC__) we cannot include pybind11-bearing headers.
// Forward-declare only the enum types we need from those headers.
#ifdef __HIPCC__

namespace hoomd { namespace sph {

#ifndef KERNELTYPES
#define KERNELTYPES (wendlandc2)(wendlandc4)(wendlandc6)(quintic)(cubicspline)
#endif
enum SmoothingKernelType { wendlandc2=0, wendlandc4, wendlandc6, quintic, cubicspline };

#ifndef SEQTYPES
#define SEQTYPES (linear)(tait)
#endif
enum StateEquationType { linear=0, tait };

enum NonNewtonianModel
    {
    NEWTONIAN=0,
    POWERLAW,
    CARREAU,
    BINGHAM,
    HERSCHELBULKLEY
    };

}} // namespace hoomd::sph

#else
// Host compilation: use the authoritative declarations from the project headers
#include "SmoothingKernel.h"
#include "StateEquations.h"
#include "EvaluationMethodDefinition.h"
#endif // __HIPCC__

namespace hoomd
{
namespace sph
{

// =========================================================================
// POD parameter structs — usable on both host and device
// =========================================================================

struct SPHKernelDevParams
    {
    Scalar alpha;       //!< Kernel normalisation constant
    Scalar self_density;//!< Kernel self-contribution factor (w0 = self_density * alpha/h^3)
    Scalar ch;          //!< Constant smoothing length (unused when const_slength==0)
    Scalar rcutsq;      //!< Cutoff radius squared (unused when const_slength==0)
    int    const_slength;//!< 1 if constant smoothing length, 0 if variable
    Scalar kappa;       //!< kappa factor (cutoff = kappa * h); used for rcutsq when variable
    };

struct SPHEOSDevParams
    {
    Scalar rho0; //!< Reference (rest) density
    Scalar c;    //!< Speed of sound
    Scalar bp;   //!< Background pressure
    };

struct SPHNNViscParams
    {
    Scalar mu;       //!< Newtonian dynamic viscosity
    Scalar K;        //!< Power-law / H-B consistency index (or Bingham plastic viscosity)
    Scalar n;        //!< Power-law / Carreau / H-B exponent
    Scalar mu0;      //!< Carreau zero-shear viscosity
    Scalar muinf;    //!< Carreau infinite-shear viscosity
    Scalar lambda_NN;//!< Carreau relaxation time
    Scalar tauy;     //!< Yield stress (Bingham / H-B)
    Scalar m_reg;    //!< Papanastasiou regularisation parameter
    Scalar mu_min;   //!< Lower viscosity clamp
    int    model;    //!< NonNewtonianModel cast to int
    };

namespace kernel
{

// =========================================================================
// Smoothing kernel device functions
// =========================================================================

template<SmoothingKernelType KT_>
__device__ __forceinline__ Scalar sph_normfactor(Scalar alpha, Scalar h)
    {
    return alpha / (h * h * h);
    }

template<SmoothingKernelType KT_>
__device__ __forceinline__ Scalar sph_w0(Scalar alpha, Scalar self_density, Scalar h)
    {
    return self_density * sph_normfactor<KT_>(alpha, h);
    }

// --- wij specialisations ---

template<SmoothingKernelType KT_>
__device__ __forceinline__ Scalar sph_wij(Scalar alpha, Scalar h, Scalar r);

template<>
__device__ __forceinline__ Scalar sph_wij<wendlandc2>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return norm * t * t * t * t * (Scalar(1.0) + Scalar(2.0) * q);
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_wij<wendlandc4>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return norm * t * t * t * t * t * t * (Scalar(8.75) * q * q + Scalar(9.0) * q + Scalar(3.0));
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_wij<wendlandc6>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return norm * t * t * t * t * t * t * t * t
               * (Scalar(4.0) * q * q * q + Scalar(6.25) * q * q + Scalar(4.0) * q + Scalar(1.0));
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_wij<quintic>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    Scalar norm = alpha / (h * h * h);
    if (q >= Scalar(0) && q <= Scalar(1))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0 * t0;
        Scalar t1 = (Scalar(2.0) - q); t1 = t1 * t1 * t1 * t1 * t1;
        Scalar t2 = (Scalar(1.0) - q); t2 = t2 * t2 * t2 * t2 * t2;
        return norm * (t0 - Scalar(6.0) * t1 + Scalar(15.0) * t2);
        }
    else if (q > Scalar(1) && q <= Scalar(2))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0 * t0;
        Scalar t1 = (Scalar(2.0) - q); t1 = t1 * t1 * t1 * t1 * t1;
        return norm * (t0 - Scalar(6.0) * t1);
        }
    else if (q > Scalar(2) && q <= Scalar(3))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0 * t0;
        return norm * t0;
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_wij<cubicspline>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    Scalar norm = alpha / (h * h * h);
    if (q >= Scalar(0) && q <= Scalar(1))
        return norm * (Scalar(1.0) - Scalar(1.5) * q * q + Scalar(0.75) * q * q * q);
    else if (q > Scalar(1) && q <= Scalar(2))
        {
        Scalar t = Scalar(2.0) - q;
        return norm * Scalar(0.25) * t * t * t;
        }
    return Scalar(0.0);
    }

// --- dwijdr specialisations ---

template<SmoothingKernelType KT_>
__device__ __forceinline__ Scalar sph_dwijdr(Scalar alpha, Scalar h, Scalar r);

template<>
__device__ __forceinline__ Scalar sph_dwijdr<wendlandc2>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return -(norm / h) * t * t * t * Scalar(5.0) * q;
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_dwijdr<wendlandc4>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return -(norm / h) * t * t * t * t * t * (Scalar(35.0) * q * q + Scalar(14.0) * q);
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_dwijdr<wendlandc6>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    if (q >= Scalar(0) && q <= Scalar(2))
        {
        Scalar norm = alpha / (h * h * h);
        Scalar t = Scalar(1.0) - Scalar(0.5) * q;
        return -(norm / h) * t * t * t * t * t * t * t
               * Scalar(11.0) * q * (Scalar(2.0) * q * q + Scalar(1.75) * q + Scalar(0.5));
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_dwijdr<quintic>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    Scalar norm = alpha / (h * h * h);
    if (q >= Scalar(0) && q <= Scalar(1))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0;
        Scalar t1 = (Scalar(2.0) - q); t1 = t1 * t1 * t1 * t1;
        Scalar t2 = (Scalar(1.0) - q); t2 = t2 * t2 * t2 * t2;
        return -(norm / h) * (Scalar(5.0) * t0 - Scalar(30.0) * t1 + Scalar(75.0) * t2);
        }
    else if (q > Scalar(1) && q <= Scalar(2))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0;
        Scalar t1 = (Scalar(2.0) - q); t1 = t1 * t1 * t1 * t1;
        return -(norm / h) * (Scalar(5.0) * t0 - Scalar(30.0) * t1);
        }
    else if (q > Scalar(2) && q <= Scalar(3))
        {
        Scalar t0 = (Scalar(3.0) - q); t0 = t0 * t0 * t0 * t0;
        return -(norm / h) * Scalar(5.0) * t0;
        }
    return Scalar(0.0);
    }

template<>
__device__ __forceinline__ Scalar sph_dwijdr<cubicspline>(Scalar alpha, Scalar h, Scalar r)
    {
    Scalar q = r / h;
    Scalar norm = alpha / (h * h * h);
    if (q >= Scalar(0) && q <= Scalar(1))
        return -(norm / h) * (Scalar(3.0) * q - Scalar(2.25) * q * q);
    else if (q > Scalar(1) && q <= Scalar(2))
        {
        Scalar t = Scalar(2.0) - q;
        return -(norm / h) * Scalar(0.75) * t * t;
        }
    return Scalar(0.0);
    }

// =========================================================================
// EOS device functions
// =========================================================================

template<StateEquationType SET_>
__device__ __forceinline__ Scalar sph_pressure(Scalar rho0, Scalar c, Scalar bp, Scalar rho);

template<>
__device__ __forceinline__ Scalar sph_pressure<tait>(Scalar rho0, Scalar c, Scalar bp, Scalar rho)
    {
    return ((rho0 * c * c) / Scalar(7.0)) * (pow(rho / rho0, Scalar(7.0)) - Scalar(1.0)) + bp;
    }

template<>
__device__ __forceinline__ Scalar sph_pressure<linear>(Scalar rho0, Scalar c, Scalar bp, Scalar rho)
    {
    return c * c * (rho - rho0) + bp;
    }

template<StateEquationType SET_>
__device__ __forceinline__ Scalar sph_dpressuredrho(Scalar rho0, Scalar c, Scalar rho);

template<>
__device__ __forceinline__ Scalar sph_dpressuredrho<tait>(Scalar rho0, Scalar c, Scalar rho)
    {
    return c * c * pow(rho / rho0, Scalar(6.0));
    }

template<>
__device__ __forceinline__ Scalar sph_dpressuredrho<linear>(Scalar rho0, Scalar c, Scalar rho)
    {
    return c * c;
    }

template<StateEquationType SET_>
__device__ __forceinline__ Scalar sph_density_from_p(Scalar rho0, Scalar c, Scalar bp, Scalar p);

template<>
__device__ __forceinline__ Scalar sph_density_from_p<tait>(Scalar rho0, Scalar c, Scalar bp, Scalar p)
    {
    return rho0 * pow((p - bp) * (Scalar(7.0) / (rho0 * c * c)) + Scalar(1.0),
                      Scalar(0.14285714285714285));
    }

template<>
__device__ __forceinline__ Scalar sph_density_from_p<linear>(Scalar rho0, Scalar c, Scalar bp, Scalar p)
    {
    return (p - bp) / (c * c) + rho0;
    }

// VRD variants (variable reference density, per-particle rho0_local)
template<StateEquationType SET_>
__device__ __forceinline__ Scalar sph_pressure_vrd(Scalar c, Scalar bp,
                                                    Scalar rho0_local, Scalar rho);

template<>
__device__ __forceinline__ Scalar sph_pressure_vrd<tait>(Scalar c, Scalar bp,
                                                          Scalar rho0_local, Scalar rho)
    {
    return ((rho0_local * c * c) / Scalar(7.0))
           * (pow(rho / rho0_local, Scalar(7.0)) - Scalar(1.0)) + bp;
    }

template<>
__device__ __forceinline__ Scalar sph_pressure_vrd<linear>(Scalar c, Scalar bp,
                                                            Scalar rho0_local, Scalar rho)
    {
    return c * c * (rho - rho0_local) + bp;
    }

template<StateEquationType SET_>
__device__ __forceinline__ Scalar sph_dpressure_vrd_drho(Scalar c, Scalar rho0_local, Scalar rho);

template<>
__device__ __forceinline__ Scalar sph_dpressure_vrd_drho<tait>(Scalar c, Scalar rho0_local,
                                                                Scalar rho)
    {
    return c * c * pow(rho / rho0_local, Scalar(6.0));
    }

template<>
__device__ __forceinline__ Scalar sph_dpressure_vrd_drho<linear>(Scalar c, Scalar rho0_local,
                                                                   Scalar rho)
    {
    return c * c;
    }

// =========================================================================
// Non-Newtonian viscosity device function
// =========================================================================

__device__ __forceinline__ Scalar sph_nn_viscosity(
    Scalar mu_base,
    Scalar gamma_dot,
    int    model,
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
        case 0: // NEWTONIAN
            return mu_base;

        case 1: // POWERLAW
            {
            Scalar gdot = gamma_dot > Scalar(1e-12) ? gamma_dot : Scalar(1e-12);
            Scalar mu_eff = K * pow(gdot, n - Scalar(1.0));
            return mu_eff > mu_min ? mu_eff : mu_min;
            }

        case 2: // CARREAU
            {
            Scalar lg = lambda_NN * gamma_dot;
            return muinf + (mu0 - muinf)
                   * pow(Scalar(1.0) + lg * lg, Scalar(0.5) * (n - Scalar(1.0)));
            }

        case 3: // BINGHAM
            {
            Scalar mu_eff;
            if (gamma_dot < Scalar(1e-12))
                mu_eff = K + tauy * m_reg;
            else
                mu_eff = K + tauy * (Scalar(1.0) - exp(-m_reg * gamma_dot)) / gamma_dot;
            return mu_eff > mu_min ? mu_eff : mu_min;
            }

        case 4: // HERSCHELBULKLEY
            {
            Scalar gdot = gamma_dot > Scalar(1e-12) ? gamma_dot : Scalar(1e-12);
            Scalar mu_eff = K * pow(gdot, n - Scalar(1.0))
                            + tauy * (Scalar(1.0) - exp(-m_reg * gdot)) / gdot;
            return mu_eff > mu_min ? mu_eff : mu_min;
            }

        default:
            return mu_base;
        }
    }

// =========================================================================
// Type-check helpers
// =========================================================================

__device__ __forceinline__ bool sph_checksolid(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & 2; // SolidFluidTypeBit::SOLID = 1<<1 = 2
    }

__device__ __forceinline__ bool sph_checkfluid(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & 4; // SolidFluidTypeBit::FLUID = 1<<2 = 4
    }

__device__ __forceinline__ bool sph_checkfluid1(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & 8; // SolidFluidTypeBit::FLUID1 = FLUID<<1 = 8
    }

__device__ __forceinline__ bool sph_checkfluid2(const unsigned int* type_props, Scalar mytype)
    {
    return type_props[__scalar_as_int(mytype)] & 16; // SolidFluidTypeBit::FLUID2 = FLUID<<2 = 16
    }

} // namespace kernel
} // namespace sph
} // namespace hoomd

#endif // __SPH_DEVICE_FUNCTIONS_CUH__
