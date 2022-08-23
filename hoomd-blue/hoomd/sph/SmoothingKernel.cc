// Copyright (c) 2009-2022 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "SmoothingKernel.h"

// #include <boost/python.hpp>
// using namespace boost::python;

// #include <boost/bind.hpp>
// using namespace boost;

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
template<>
SmoothingKernel<wendlandc4>::SmoothingKernel()
    : m_kappa(Scalar(2.0)), m_self_density(Scalar(3.0)), m_alpha(Scalar(0.198192948))
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
    return m_self_density*normalizationfactor(h);
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
template<>
Scalar SmoothingKernel<cubicspline>::wij(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar w = Scalar(0.0);
        if( q >= Scalar(0.0) && q < Scalar(1) )
        {
            // (alpha/h^D)*( 1.0 - 1.5*q^2 + 0.75*q^3 );
            w = normalizationfactor(h)* ( Scalar(1) - Scalar(1.5)*q*q + Scalar(0.75)*q*q*q );
        }
        else if ( q >= Scalar(1) && q < Scalar(2) )
        {
            // (alpha/h^D)*( 0.25*(2.0-q)^3 );
            w = normalizationfactor(h)* ( Scalar(0.25)* (pow((Scalar(2.0)-q),Scalar(3)))  );
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
template<>
Scalar SmoothingKernel<wendlandc4>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);
        if ( q >= 0 && q < Scalar(2) )
            {
            // (alpha/h^D+1) * (1-0.5)^5 * q * -(35*q + 14)
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

template<>
Scalar SmoothingKernel<cubicspline>::dwijdr(const Scalar h, const Scalar rij)
{
        Scalar q = rij / h;
        Scalar dwdr = Scalar(0.0);

        if( q >= Scalar(0.0) && q < Scalar(1) )
        {
            // (alpha/(h^D+1))*-5*( (3.0-q) + 2.25*q^2 );
            dwdr = -(normalizationfactor(h)/h)*( Scalar(3)*q + Scalar(2.25)*q*q );
        }else if ( q >= Scalar(1) && q < Scalar(2) )
        {
            // (alpha/(h^D+1))*-5*( 0.75*(2.0-q)^2 );
            dwdr = -(normalizationfactor(h)/h)*( Scalar(0.75)*(Scalar(2.0)-q)*(Scalar(2.0)-q) );

        }

        return dwdr;
}

//! Explicit template instantiations
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc2>::setAlpha;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc2>::setSelfDensity;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc2>::setKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc2>::getKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc2>::w0;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc2>::normalizationfactor;

// template void PYBIND11_EXPORT SmoothingKernel<wendlandc4>::setAlpha;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc4>::setSelfDensity;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc4>::setKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc4>::getKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc4>::w0;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc4>::normalizationfactor;

// template void PYBIND11_EXPORT SmoothingKernel<wendlandc6>::setAlpha;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc6>::setSelfDensity;
// template void PYBIND11_EXPORT SmoothingKernel<wendlandc6>::setKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc6>::getKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc6>::w0;
// template Scalar PYBIND11_EXPORT SmoothingKernel<wendlandc6>::normalizationfactor;

// template void PYBIND11_EXPORT SmoothingKernel<quintic>::setAlpha;
// template void PYBIND11_EXPORT SmoothingKernel<quintic>::setSelfDensity;
// template void PYBIND11_EXPORT SmoothingKernel<quintic>::setKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<quintic>::getKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<quintic>::w0;
// template Scalar PYBIND11_EXPORT SmoothingKernel<quintic>::normalizationfactor;

// template void PYBIND11_EXPORT SmoothingKernel<cubicspline>::setAlpha;
// template void PYBIND11_EXPORT SmoothingKernel<cubicspline>::setSelfDensity;
// template void PYBIND11_EXPORT SmoothingKernel<cubicspline>::setKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<cubicspline>::getKernelKappa;
// template Scalar PYBIND11_EXPORT SmoothingKernel<cubicspline>::w0;
// template Scalar PYBIND11_EXPORT SmoothingKernel<cubicspline>::normalizationfactor;

// template<> std::string get_kernel_name<wendlandc2>()
// {return "WC2";}
// template<> std::string get_kernel_name<wendlandc4>()
// {return "WC4";}
// template<> std::string get_kernel_name<wendlandc6>()
// {return "WC6";}
// template<> std::string get_kernel_name<quintic>()
// {return "Q";}
// template<> std::string get_kernel_name<cubicspline>()
// {return "CS";}

//! Explicit template instantiations

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
     pybind11::class_<SmoothingKernel<cubicspline>, std::shared_ptr<SmoothingKernel<cubicspline>>>(m, "Cubicspline")
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





 // void export_WendlandC4(pybind11::module& m)
 //     {
 //     pybind11::class_<SmoothingKernel<wendlandc4>, std::shared_ptr<SmoothingKernel<wendlandc4>>>(m, "WendlandC4")
 //         .def("getKernelKappa", &SmoothingKernel<wendlandc4>::getKernelKappa)
 //         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
 //         .def("EvalKernel", &SmoothingKernel<wendlandc4>::wij)
 //         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc4>::dwijdr)
 //         ;
 //     }
 // void export_WendlandC6(pybind11::module& m)
 //     {
 //     pybind11::class_<SmoothingKernel<wendlandc6>, std::shared_ptr<SmoothingKernel<wendlandc6>>>(m, "WendlandC6")
 //         .def("getKernelKappa", &SmoothingKernel<wendlandc6>::getKernelKappa)
 //         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
 //         .def("EvalKernel", &SmoothingKernel<wendlandc6>::wij)
 //         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc6>::dwijdr)
 //         ;
 //     }
 // void export_Quintic(pybind11::module& m)
 //     {
 //     pybind11::class_<SmoothingKernel<quintic>, std::shared_ptr<SmoothingKernel<quintic>>>(m, "Quintic")
 //         .def("getKernelKappa", &SmoothingKernel<quintic>::getKernelKappa)
 //         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
 //         .def("EvalKernel", &SmoothingKernel<quintic>::wij)
 //         .def("EvalKernelDerivative", &SmoothingKernel<quintic>::dwijdr)
 //         ;
 //     }
 // void export_CubicSpline(pybind11::module& m)
 //     {
 //     pybind11::class_<SmoothingKernel<cubicspline>, std::shared_ptr<SmoothingKernel<cubicspline>>>(m, "CubicSpline")
 //         .def("getKernelKappa", &SmoothingKernel<cubicspline>::getKernelKappa)
 //         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
 //         .def("EvalKernel", &SmoothingKernel<cubicspline>::wij)
 //         .def("EvalKernelDerivative", &SmoothingKernel<cubicspline>::dwijdr)
 //         ;
 //     }
 } // end namespace detail

// #define INST_getKernelKappa(r,data,k) template Scalar SmoothingKernel<k>::getKernelKappa();
// BOOST_PP_SEQ_FOR_EACH(INST_getKernelKappa, ~, KERNELTYPES);
// #undef INST_getKernelKappa

// #define INST_w0(r,data,k) template Scalar SmoothingKernel<k>::w0(const Scalar h);
// BOOST_PP_SEQ_FOR_EACH(INST_w0, ~, KERNELTYPES);
// #undef INST_w0

} // end namespace sph
} // end namespace hoomd


