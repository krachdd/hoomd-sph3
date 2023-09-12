/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/

#include "SmoothingKernel.h"

#include <pybind11/pybind11.h>

using namespace std;

namespace hoomd
{
namespace sph
{
/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param nlist Neighbor list instance. Must not be NULL.
    \post SmoothingKernel base class is constructed.
*/

template<SmoothingKernelType KT_>
SmoothingKernel<KT_>::SmoothingKernel()
    : m_kappa(Scalar(0.0)), m_self_density(Scalar(0.0)), m_alpha(Scalar(0.0))
    {
    // sanity check
    assert(this->m_sysdef);
    assert(this->m_pdata);
    }

template<>
SmoothingKernel<wendlandc2>::SmoothingKernel()
    : m_kappa(Scalar(2.0)), m_self_density(Scalar(1.0)), m_alpha(Scalar(0.41778174))
    {
    }

/*
Wendland C4 Kernel
m_alpha = 495./(256. * PI),
m_kappa = 2.0 
m_self_density = 3.0
see https://pysph.readthedocs.io/en/latest/reference/kernels.html
*/
// template<>
// SmoothingKernel<wendlandc4>::SmoothingKernel()
//     : m_kappa(Scalar(2.0)), m_self_density(Scalar(3.0)), m_alpha(Scalar(0.6154820064881891))
//     {
//     }
template<>
SmoothingKernel<wendlandc4>::SmoothingKernel()
    : m_kappa(Scalar(2.0)), m_self_density(Scalar(3.0)), m_alpha(Scalar(0.20516066882939635))
    {
    }
template<>
SmoothingKernel<wendlandc6>::SmoothingKernel()
    : m_kappa(Scalar(2.0)), m_self_density(Scalar(1.0)), m_alpha(Scalar(0.8486191301579575))
    {
    }
template<>
SmoothingKernel<quintic>::SmoothingKernel()
    : m_kappa(Scalar(3.0)), m_self_density(Scalar(66.0)), m_alpha(Scalar(0.0026525825))
    {
    }

/*
Cubic Spline
m_alpha = 1./PI, 
m_kappa = 2.0,
m_self_density check properties if needed for n_density or density renormalization
see: J. Monaghan, Smoothed Particle Hydrodynamics, “Annual Review of Astronomy and Astrophysics”, 30 (1992), pp. 543-574.
and https://pysph.readthedocs.io/en/latest/reference/kernels.html
*/
template<>
SmoothingKernel<cubicspline>::SmoothingKernel()
    : m_kappa(Scalar(2.0)), m_self_density(Scalar(1.0)), m_alpha(Scalar(0.31830987))
    {
    }



template<SmoothingKernelType KT_>
Scalar SmoothingKernel<KT_>::wij(const Scalar h, const Scalar rij)
{
    return Scalar(0.0);
}

template<SmoothingKernelType KT_>
Scalar SmoothingKernel<KT_>::dwijdr(const Scalar h, const Scalar rij)
{
    return Scalar(0.0);
}

//! Set kernel normalization factor
template<SmoothingKernelType KT_>
void SmoothingKernel<KT_>::setAlpha(const Scalar alpha)
{
    m_alpha = alpha;
}

//! Set kernel self-density
template<SmoothingKernelType KT_>
void SmoothingKernel<KT_>::setSelfDensity(const Scalar self_density)
{
    m_self_density = self_density;
}

//! Set kernel kappa
template<SmoothingKernelType KT_>
void SmoothingKernel<KT_>::setKernelKappa(const Scalar kappa)
{
    m_kappa = kappa;
}

//!Set neighbor list instance
// template<SmoothingKernelType KT_>
// void SmoothingKernel<KT_>::setNeighborList(const boost::shared_ptr<nsearch::NeighborList> nlist)
// {
//     m_nlist = nlist;
// }

//! Get kernel kappa
template<SmoothingKernelType KT_>
Scalar SmoothingKernel<KT_>::getKernelKappa()
{
    return m_kappa;
}

//! Return kernel self density
/*! \param h Smoothing length parameter
*/
template<SmoothingKernelType KT_>
Scalar SmoothingKernel<KT_>::w0(const Scalar h)
{
    return m_self_density * normalizationfactor(h);
}

//! Return kernel normalization factor
/*! \param h Smoothing length parameter
*/
template<SmoothingKernelType KT_>
Scalar SmoothingKernel<KT_>::normalizationfactor(const Scalar h)
{
    return m_alpha / (h*h*h);
}

template<>
Scalar SmoothingKernel<wendlandc2>::wij(const Scalar h, const Scalar rij)
{
    Scalar q = rij / h;
    Scalar w = Scalar(0.0);
    if ( q >= 0 && q < Scalar(2) )
    {
        // (alpha/h^D) * (1-0.5q)^4 * ( 1.0 + 2*q )
        w = normalizationfactor(h)*( pow(Scalar(1)-Scalar(0.5)*q,Scalar(4))*(Scalar(1) + Scalar(2)*q) );
    }
    return w;
}
/* checked DK 10/2022 */
// template<>
// Scalar SmoothingKernel<wendlandc4>::wij(const Scalar h, const Scalar rij)
// {
//         Scalar q = rij / h;
//         Scalar w = Scalar(0.0);
//         if ( q >= 0 && q <= 2.0)
//             {
//             // (alpha/h^D) * (1-0.5q)^6 * (35/12*q^2 +3*q + 1)
//             w = normalizationfactor(h)*pow(1.0-(0.5*q),6.0)*(2.9166666*q*q + 3.0*q + 1.0);
//             }
//         return w;
// }
template<>
Scalar SmoothingKernel<wendlandc4>::wij(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar w = Scalar(0.0);
        if ( q >= 0 && q < Scalar(2) )
            {
            // (alpha/h^D) * (1-0.5q)^6 * ( 8.75*q^2 + 9*q + 3 )
            w = normalizationfactor(h)*pow(Scalar(1)-Scalar(0.5)*q,Scalar(6))*(Scalar(8.75)*q*q+Scalar(9)*q+Scalar(3));
            }
        return w;
}

template<>
Scalar SmoothingKernel<wendlandc6>::wij(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar w = Scalar(0.0);
        if ( q >= 0 && q < Scalar(2) )
            {
            // (alpha/h^D) * (1-0.5q)^8 * ( 4*q^3 + 6.25*q^2 + 4*q + 1 )
            w = normalizationfactor(h)*( pow(Scalar(1)-Scalar(0.5)*q,Scalar(8)) * (Scalar(4)*q*q*q+Scalar(6.25)*q*q+Scalar(4)*q+Scalar(1)) );
            }
        return w;
}
template<>
Scalar SmoothingKernel<quintic>::wij(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar w = Scalar(0.0);
        if( q >= Scalar(0.0) && q < Scalar(1) )
        {
            // (alpha/h^D)*( (3.0-q)^5 - 6.0*(2.0-q)^5 + 15.0*(1.0-q)^5 );
            w = normalizationfactor(h)* ( pow(Scalar(3)-q,Scalar(5)) - Scalar(6)*pow(Scalar(2)-q,Scalar(5)) + Scalar(15)*pow(Scalar(1)-q,Scalar(5)) );
        }else if ( q >= Scalar(1) && q < Scalar(2) )
        {
            // (alpha/h^D)*( (3.0-q)^5 - 6.0*(2.0-q)^5 ;
            w = normalizationfactor(h)* ( pow(Scalar(3)-q,Scalar(5)) - Scalar(6)*pow(Scalar(2)-q,Scalar(5)) );
        }else if ( q >= Scalar(2) && q < Scalar(3) )
        {
            // (alpha/h^D)*( (3.0-q)^5 ) ;
            w = normalizationfactor(h)* ( pow(Scalar(3)-q,Scalar(5)) );
        }

        return w;
}

/* checked DK 10/2022 */
template<>
Scalar SmoothingKernel<cubicspline>::wij(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar w = Scalar(0.0);
        if( q >= 0.0 && q <= 1.0)
        {
            // (alpha/h^3)*( 1.0 - 1.5*q^2 + 0.75*q^3 );
            w = normalizationfactor(h)* ( 1.0 - 1.5*q*q + 0.75*q*q*q);
        }
        else if ( q > 1.0 && q <= 2.0 )
        {
            // (alpha/h^D)*( 0.25*(2.0-q)^3 );
            w = normalizationfactor(h)* ( 0.25* (2.0-q) * (2.0-q) * (2.0-q));
        }

        return w;
}

template<>
Scalar SmoothingKernel<wendlandc2>::dwijdr(const Scalar h, const Scalar rij)
{
    Scalar q = rij / h;
    Scalar dwdr = Scalar(0.0);
    if ( q >= 0 && q < Scalar(2) )
    {
        // (alpha/h^D+1) * 5*(1-0.5)^3
        dwdr = -(normalizationfactor(h)/h)*( Scalar(5)*pow(Scalar(1)-Scalar(0.5)*q,Scalar(3)) );
    }
    return dwdr;
}

/* checked DK 10/2022 
adapted limits, removed additional datatyping
*/
// template<>
// Scalar SmoothingKernel<wendlandc4>::dwijdr(const Scalar h, const Scalar rij)
// {
//         Scalar q = rij / h;
//         Scalar dwdr = 0.0;
//         if ( q >= 0 && q <=  2.0 )
//             {
//             // - (alpha/h^D+1) * 7./96. * (2-q)^5 *q*(5q +2)
//             dwdr = -(normalizationfactor(h)/h)* 7./96. * (2.0-q) * (2.0-q) * (2.0-q) * (2.0-q) * (2.0-q) * q * (5.0*q +2);
//             }
//         return dwdr;
// }
template<>
Scalar SmoothingKernel<wendlandc4>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);
        if ( q >= 0 && q < Scalar(2) && rij > 1e-12)
            {
            // (alpha/h^D+1) * (1-0.5*q)^5 * q * -(35*q + 14)
            dwdr = -(normalizationfactor(h)/h)*pow(Scalar(1)-Scalar(0.5)*q,Scalar(5))*q*(Scalar(35)*q+Scalar(14));
            }
        return dwdr;
}

template<>
Scalar SmoothingKernel<wendlandc6>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);
        if ( q >= 0 && q < Scalar(2) )
            {
            // (alpha/h^D+1) * 11*(1-0.5)^7 * (2*q^2 + 1.75*q + 0.5)
            dwdr = -(normalizationfactor(h)/h)* Scalar(11)*q*( pow(Scalar(1)-Scalar(0.5)*q,Scalar(7))*(Scalar(2)*q*q+Scalar(1.75)*q+Scalar(0.5)) );
            }
        return dwdr;
}
template<>
Scalar SmoothingKernel<quintic>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);

        if( q >= Scalar(0.0) && q < Scalar(1) )
        {
            // (alpha/(h^D+1))*-5*( (3.0-q)^4 - 6.0*(2.0-q)^4 + 15.0*(1.0-q)^4 );
            dwdr = -(normalizationfactor(h)/h)*Scalar(5)*( pow(Scalar(3)-q,Scalar(4)) - Scalar(6)*( pow(Scalar(2)-q,Scalar(4)) + Scalar(15)*pow(Scalar(1)-q,Scalar(3)) ) );
        }else if ( q >= Scalar(1) && q < Scalar(2) )
        {
            // (alpha/(h^D+1))*-5*( (3.0-q)^4 - 6.0*(2.0-q)^4 );
            dwdr = -(normalizationfactor(h)/h)*Scalar(5)*( pow(Scalar(3)-q,Scalar(4)) - Scalar(6)*( pow(Scalar(2)-q,Scalar(4)) ) );
        }else if ( q >= Scalar(2) && q < Scalar(3) )
        {
            // (alpha/(h^D+1))*-5*( (3.0-q)^4 );
            dwdr = -(normalizationfactor(h)/h)*Scalar(5)*( pow(Scalar(3)-q,Scalar(4)) );
        }

        return dwdr;
}


/* checked DK 10/2022 
adapted limits, removed additional datatyping
*/
template<>
Scalar SmoothingKernel<cubicspline>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);

        if( q >= 0.0 && q <= 1.0 )
        {
            // (alpha/(h^D+1))*-1*( (3.0-q) + 2.25*q^2 );
            dwdr = -(normalizationfactor(h)/h)*( 3*q + 2.25*q*q );
        }
        else if ( q > 1.0 && q <= 2.0 )
        {
            // (alpha/(h^D+1))*-1*( 0.75*(2.0-q)^2 );
            dwdr = -(normalizationfactor(h)/h)*( 0.75*(2.0-q)*(2.0-q) );
        }

        return dwdr;
}

namespace detail
{
void export_WendlandC2(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc2>, std::shared_ptr<SmoothingKernel<wendlandc2>>>(m, "WendlandC2")
         .def(pybind11::init<>())
         .def("getKernelKappa", &SmoothingKernel<wendlandc2>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc2>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc2>::dwijdr)
         .def("setAlpha", &SmoothingKernel<wendlandc2>::setAlpha)
         .def("setSelfDensity", &SmoothingKernel<wendlandc2>::setSelfDensity)
         .def("setKernelKappa", &SmoothingKernel<wendlandc2>::setKernelKappa)
         .def("w0", &SmoothingKernel<wendlandc2>::w0)
         .def("normalizationfactor", &SmoothingKernel<wendlandc2>::normalizationfactor)

         ;
     }

void export_WendlandC4(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc4>, std::shared_ptr<SmoothingKernel<wendlandc4>>>(m, "WendlandC4")
         .def(pybind11::init<>())
         .def("getKernelKappa", &SmoothingKernel<wendlandc4>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc4>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc4>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc4>::dwijdr)
         .def("setAlpha", &SmoothingKernel<wendlandc4>::setAlpha)
         .def("setSelfDensity", &SmoothingKernel<wendlandc4>::setSelfDensity)
         .def("setKernelKappa", &SmoothingKernel<wendlandc4>::setKernelKappa)
         .def("w0", &SmoothingKernel<wendlandc4>::w0)
         .def("normalizationfactor", &SmoothingKernel<wendlandc4>::normalizationfactor)

         ;
     }

void export_WendlandC6(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc6>, std::shared_ptr<SmoothingKernel<wendlandc6>>>(m, "WendlandC6")
         .def(pybind11::init<>())
         .def("getKernelKappa", &SmoothingKernel<wendlandc6>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc6>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc6>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc6>::dwijdr)
         .def("setAlpha", &SmoothingKernel<wendlandc6>::setAlpha)
         .def("setSelfDensity", &SmoothingKernel<wendlandc6>::setSelfDensity)
         .def("setKernelKappa", &SmoothingKernel<wendlandc6>::setKernelKappa)
         .def("w0", &SmoothingKernel<wendlandc6>::w0)
         .def("normalizationfactor", &SmoothingKernel<wendlandc6>::normalizationfactor)

         ;
     }

void export_Quintic(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<quintic>, std::shared_ptr<SmoothingKernel<quintic>>>(m, "Quintic")
         .def(pybind11::init<>())
         .def("getKernelKappa", &SmoothingKernel<quintic>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<quintic>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<quintic>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<quintic>::dwijdr)
         .def("setAlpha", &SmoothingKernel<quintic>::setAlpha)
         .def("setSelfDensity", &SmoothingKernel<quintic>::setSelfDensity)
         .def("setKernelKappa", &SmoothingKernel<quintic>::setKernelKappa)
         .def("w0", &SmoothingKernel<quintic>::w0)
         .def("normalizationfactor", &SmoothingKernel<quintic>::normalizationfactor)

         ;
     }

void export_CubicSpline(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<cubicspline>, std::shared_ptr<SmoothingKernel<cubicspline>>>(m, "CubicSpline")
         .def(pybind11::init<>())
         .def("getKernelKappa", &SmoothingKernel<cubicspline>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<cubicspline>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<cubicspline>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<cubicspline>::dwijdr)
         .def("setAlpha", &SmoothingKernel<cubicspline>::setAlpha)
         .def("setSelfDensity", &SmoothingKernel<cubicspline>::setSelfDensity)
         .def("setKernelKappa", &SmoothingKernel<cubicspline>::setKernelKappa)
         .def("w0", &SmoothingKernel<cubicspline>::w0)
         .def("normalizationfactor", &SmoothingKernel<cubicspline>::normalizationfactor)
         ;
     }

} // end namespace detail
} // end namespace sph
} // end namespace hoomd


