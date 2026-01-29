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
#include <pybind11/stl_bind.h>


#ifndef __SPH_SMOOTHING_KERNEL_H__
#define __SPH_SMOOTHING_KERNEL_H__

/*! \file SmoothingKernel.h
    \brief Declares a base class for smoothing kernels.

    DK: Included the self density directly in alpha. 
    e.g. WendlandC4:
    actually: m_alpha = 165./(256. * PI)
    changed : m_alpha = 495./(256. * PI)
    which is m_self_density * 165./(256. * PI)
    Therefore the w0 has to equal the normalisation factor
    normalisationfactor = m_alpha/(h*h*h) 
    instead of :
    normalisationfactor = m_self_density * m_alpha/(h*h*h) 
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

    private:
        Scalar m_kappa; //!< Kernel size scaling factor
        Scalar m_self_density; //!< Kernel self-density, i.e. w(0)
        Scalar m_alpha; //!< Kernel renormalization factor
    };
    
} // end namespace sph
} // end namespace hoomd

// undefine HOSTDEVICE so we don't interfere with other headers
#undef HOSTDEVICE

#endif // #ifndef __SPH_SMOOTHING_KERNEL_H__
