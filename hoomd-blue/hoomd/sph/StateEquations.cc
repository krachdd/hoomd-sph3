// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// Maintainer: David krach 

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
    : m_bp(Scalar(0.0)), m_rho0(Scalar(0.0)), m_c(Scalar(0.0))
    {
        // Default values
        m_params_set = false;
    }

template<StateEquationType SET_>
void StateEquation<SET_>::setParams(Scalar rho0, Scalar c, Scalar bpfactor)
    {
        m_rho0 = rho0;
        m_c = c;
        m_bpfactor = bpfactor;
        m_bp = m_bpfactor*m_rho0*m_c*m_c;
        m_params_set = true;

    }
template<StateEquationType SET_>
void StateEquation<SET_>::setBackPressure(Scalar bp)
    {
        m_bp = bp;
        m_params_set = true;
    }


template<>
Scalar StateEquation<tait>::Pressure(const Scalar rho)
    {
        // p = (r0*c^2/7)*( (rho/r0)^7 - 1 )  + backp*rho*c^2
        return (((m_rho0*m_c*m_c)/Scalar(7))*(pow(rho/m_rho0,Scalar(7))-Scalar(1)))+m_bp;
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

template<> std::string get_SE_name<linear>()
{return "L";}
template<> std::string get_SE_name<tait>()
{return "T";}

// template StateEquation<tait>::StateEquation();
// template StateEquation<linear>::StateEquation();

// template void StateEquation<tait>::setParams(Scalar rho0, Scalar c, Scalar bpfactor);
// template void StateEquation<linear>::setParams(Scalar rho0, Scalar c, Scalar bpfactor);

// template void StateEquation<tait>::setBackPressure(Scalar bp);
// template void StateEquation<linear>::setBackPressure(Scalar bp);


namespace detail 
{

void export_StateEquation_Tait(pybind11::module& m)
    {
    pybind11::class_<StateEquation<tait>, std::shared_ptr<StateEquation<tait>>>(m, "Tait")
        .def("Pressure", &StateEquation<tait>::Pressure)
        .def("Density", &StateEquation<tait>::Density)
        .def("setParams", &StateEquation<tait>::setParams);
    }

void export_StateEquation_Linear(pybind11::module& m)
    {
    pybind11::class_<StateEquation<linear>, std::shared_ptr<StateEquation<linear>>>(m, "Linear")
        .def("Pressure", &StateEquation<linear>::Pressure)
        .def("Density", &StateEquation<linear>::Density)
        .def("setParams", &StateEquation<linear>::setParams);
    }


} // end namespace detail
} // end namespace sph
} // end namespace hoomd