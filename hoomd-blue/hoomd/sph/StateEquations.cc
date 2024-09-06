/* ---------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
----------------------------------------------------------*/ 

#include "StateEquations.h"

// #include <boost/python.hpp>
// using namespace boost::python;

// #include <boost/bind.hpp>
// using namespace boost;
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;


namespace hoomd
{
namespace sph
{

template<StateEquationType SET_>
StateEquation<SET_>::StateEquation()
    : m_bp(Scalar(0.0)), m_tvp(Scalar(0.0)), m_rho0(Scalar(0.0)), m_c(Scalar(0.0))
    {
        // Default values
        m_params_set = false;
    }

/*
This is called two times in a simulation setup. First initially to 
construct EOS with c = 0.1, secondly with c computed in sphmodel.py
*/
template<StateEquationType SET_>
void StateEquation<SET_>::setParams(Scalar rho0, Scalar c, Scalar bpfactor, Scalar tvpfactor)
    {
        m_rho0 = rho0;
        m_c = c;
        m_bpfactor = bpfactor;
        m_bp = m_bpfactor*m_rho0*m_c*m_c;
        m_tvpfactor = tvpfactor;
        m_tvp = m_tvp*m_rho0*m_c*m_c;

        m_params_set = true;
    }
    
template<StateEquationType SET_>
void StateEquation<SET_>::setBackPressure(Scalar bp)
    {
        m_bp = bp;
        m_params_set = true;
    }

template<StateEquationType SET_>
void StateEquation<SET_>::setTransportVelocityPressure(Scalar tvp)
    {
        m_tvp = tvp;
        m_params_set = true;
    }


template<>
Scalar StateEquation<tait>::Pressure(const Scalar rho)
    {
        // p = (r0*c^2/7)*( (rho/r0)^7 - 1 )  + backp*rho*c^2
        return (((m_rho0*m_c*m_c)/Scalar(7))*(pow((rho/m_rho0), Scalar(7))-Scalar(1)))+m_bp;
    }

template<>
Scalar StateEquation<tait>::Density(const Scalar p)
    {
        // rho = rho0 * [ (p-backp) * (7/rho0*c^2) +1 ]^(1/7)
        return m_rho0 * pow((p-m_bp)*(Scalar(7)/(m_rho0*m_c*m_c))+1, Scalar(0.14285714285714285));
    }


template<>
Scalar StateEquation<linear>::Pressure(const Scalar rho)
    {
        // p = c^2*(rho - r0) + backp*rho*c^2
        return m_c*m_c*(rho-m_rho0)+m_bp;
    }

template<>
Scalar StateEquation<linear>::Density(const Scalar p)
    {
        // rho = (p-backp)/c^2 + rho0
        return (p-m_bp)/(m_c*m_c) + m_rho0;
    }

namespace detail 
{


// void export_StateEquation_Tait(pybind11::module& m)
//     {
//     pybind11::class_<StateEquation<tait>, std::shared_ptr<StateEquation<tait>>>(m, "Tait")
//         .def(pybind11::init<>())
//         .def("Pressure", &StateEquation<tait>::Pressure)
//         .def("Density", &StateEquation<tait>::Density)
//         .def("setParams", &StateEquation<tait>::setParams);
//     }

// void export_StateEquation_Linear(pybind11::module& m)
//     {
//     pybind11::class_<StateEquation<linear>, std::shared_ptr<StateEquation<linear>>>(m, "Linear")
//         .def(pybind11::init<>())
//         .def("Pressure", &StateEquation<linear>::Pressure)
//         .def("Density", &StateEquation<linear>::Density)
//         .def("setParams", &StateEquation<linear>::setParams);
//     }

template<StateEquationType SET_>
void export_StateEquation(pybind11::module& m, std::string name)
{
    pybind11::class_<StateEquation<SET_>, std::shared_ptr<StateEquation<SET_>>>(m, name.c_str())
        .def(pybind11::init<>())
        .def("Pressure", &StateEquation<SET_>::Pressure)
        .def("Density", &StateEquation<SET_>::Density)
        .def("setParams", &StateEquation<SET_>::setParams);
}
} // end namespace detail

// template class PYBIND11_EXPORT StateEquation<tait>;
// template class PYBIND11_EXPORT StateEquation<linear>;

namespace detail
{
    template void export_StateEquation<tait>(pybind11::module& m, std::string name = "Tait");
    template void export_StateEquation<linear>(pybind11::module& m, std::string name = "Linear");
} // end namespace detail
} // end namespace sph
} // end namespace hoomd