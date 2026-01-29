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


#include <hoomd/HOOMDMath.h>
#include <hoomd/VectorMath.h>

#include <memory>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>


#ifndef __SPH_STATE_EQUATIONS_H__
#define __SPH_STATE_EQUATIONS_H__

/*! \file StateEquations.h
    \brief Declares state equations.
*/
#ifdef __HIPCC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif

namespace hoomd
{
namespace sph
{

#ifndef SEQTYPES
#define SEQTYPES (linear)(tait)
#endif


enum StateEquationType
{
    linear,
    tait
};

template<StateEquationType SET_>
struct StateEquation
    {
    public:
        //! Construct the state equation
        StateEquation();
        virtual ~StateEquation() {};

        /*! Set the parameters
         * \param rho0 Initial density, e.g. rest density
         * \param c Speed of sound
         * \param bpfactor Back pressure factor
         */
        void setParams(Scalar rho0, Scalar c, Scalar bpfactor, Scalar tvpfactor);

        /*! Set the parameters
         * \param bp Back pressure
         */
        void setBackPressure(Scalar bp);

        /*! Set the parameters
         * \param tv Transport Velocity Pressure
         */
        void setTransportVelocityPressure(Scalar tvp);


        // Getter and setter methods
        HOSTDEVICE Scalar getRestDensity()
            {
            return m_rho0;
            }
        HOSTDEVICE Scalar getBackgroundPressure()
            {
            return m_bp;
            }
        HOSTDEVICE Scalar getTransportVelocityPressure()
            {
            return m_tvp;
            }
        HOSTDEVICE Scalar getSpeedOfSound()
            {
            return m_c;
            }

        /*! Equation of state
         * \param rho Density
         */
        HOSTDEVICE Scalar Pressure(const Scalar rho);


        /*! Inverse equation of state
         * \param p Pressure
         */
        HOSTDEVICE Scalar Density(const Scalar p);

    protected:
        Scalar m_bpfactor; //!< Back pressure scaling factor
        Scalar m_tvpfactor; //!< Back pressure scaling factor
        Scalar m_bp; //!< Back pressure
        Scalar m_tvp; //!< Back pressure
        Scalar m_rho0; //!< Reference density
        Scalar m_c; //!< Numerical speed of sound
        bool m_params_set; //!< True if parameters are set
    };

// template<StateEquationType SET_> std::string get_SE_name();

namespace detail
{
template<StateEquationType SET_>
void export_StateEquation(pybind11::module& m, std::string name);
} // end namespace detail 
} // end namespace sph 
} // end namespace hoomd 

// undefine HOSTDEVICE so we don't interfere with other headers
#undef HOSTDEVICE


#endif // #ifndef __SPH_STATE_EQUATIONS_H__

