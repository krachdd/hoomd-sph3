// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Maintainer: David Krach 

#include <hoomd/HOOMDMath.h>
#include <hoomd/VectorMath.h>

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

// #include <boost/shared_ptr.hpp>
// #include <boost/preprocessor/seq.hpp>


#ifndef __SPH_SMOOTHING_KERNEL_H__
#define __SPH_SMOOTHING_KERNEL_H__

/*! \file SmoothingKernel.h
    \brief Declares a base class for smoothing kernels.
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

#ifndef KERNELTYPES
#define KERNELTYPES (wendlandc2)(wendlandc4)(wendlandc6)(quintic)(cubicspline)
#endif



enum SmoothingKernelType
{
    wendlandc2,
    wendlandc4,
    wendlandc6,
    quintic,
    cubicspline
};

template<SmoothingKernelType KT_>
struct PYBIND11_EXPORT SmoothingKernel
    {
    public:
        //! Construct the smoothing kernel and associate it with the neighbor list method
        SmoothingKernel();
        virtual ~SmoothingKernel() {};

        //! Return kernel evaluation
        /*! \param h Smoothing length
            \param rij Particle distance
        */
        HOSTDEVICE Scalar wij(const Scalar h, const Scalar rij);

        //! Return kernel derivative
        /*! \param h Smoothing length
            \param rij Particle distance
        */
        HOSTDEVICE Scalar dwijdr(const Scalar h, const Scalar rij);

        //! Set kernel kappa
        void setKernelKappa(const Scalar kappa);

        //! Set kernel normalization factor
        void setAlpha(const Scalar alpha);

        //! Set kernel self-density
        void setSelfDensity(const Scalar self_density);

        //! Get kernel kappa
        Scalar getKernelKappa();

        //!Set neighbor list instance
        //void setNeighborList(const boost::shared_ptr<nsearch::NeighborList> nlist);

        //! Return kernel self density
        /*! \param h Smoothing length
        */
        // HOSTDEVICE Scalar w0(const Scalar h);
        Scalar w0(const Scalar h);

        //! Return kernel normalization factor
        /*! \param h Smoothing length
        */
        // HOSTDEVICE Scalar normalizationfactor(const Scalar h);
        Scalar normalizationfactor(const Scalar h);
/*
    protected:
        const boost::shared_ptr<SystemDefinition> m_sysdef; //!< The system definition this method is associated with
        boost::shared_ptr<nsearch::NeighborList> m_nlist; //!< The neighbor list to use for the computation
        const boost::shared_ptr<ParticleData> m_pdata;  //!< The particle data this method is associated with
        boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< Stored shared ptr to the execution configuration
*/
    private:
        Scalar m_kappa; //!< Kernel size scaling factor
        Scalar m_self_density; //!< Kernel self-density, i.e. w(0)
        Scalar m_alpha; //!< Kernel renormalization factor
    };


// template<SmoothingKernelType KT_>std::string get_kernel_name();
// namespace detail
// {
// //! Exports the SPHIntegrationMethodTwoStep class to python

// void export_WendlandC2(pybind11::module& m);
// void export_WendlandC4(pybind11::module& m);
// void export_WendlandC6(pybind11::module& m);
// void export_Quintic(pybind11::module& m);
// void export_CubicSpline(pybind11::module& m);


// } // end namespace detail
} // end namespace sph
} // end namespace hoomd

// undefine HOSTDEVICE so we don't interfere with other headers
#undef HOSTDEVICE

#endif // #ifndef __SPH_SMOOTHING_KERNEL_H__
