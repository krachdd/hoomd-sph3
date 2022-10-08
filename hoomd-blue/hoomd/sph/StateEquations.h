/* ---------------------------------------------------------
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
        void setParams(Scalar rho0, Scalar c, Scalar bpfactor);

        /*! Set the parameters
         * \param bp Back pressure
         */
        void setBackPressure(Scalar bp);

        // Getter and setter methods
        HOSTDEVICE Scalar getRestDensity()
            {
            return m_rho0;
            }
        HOSTDEVICE Scalar getBackgroundPressure()
            {
            return m_bp;
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
        Scalar m_bp; //!< Back pressure
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

