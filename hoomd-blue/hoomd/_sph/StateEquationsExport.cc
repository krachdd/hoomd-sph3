// Copyright (c) 2009-2016 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "StateEquations.h"
#include <sstream>
#include <string>

// #include <boost/python.hpp>
// using namespace boost::python;

// #include <boost/bind.hpp>
// using namespace boost;

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

// namespace hoomd 
// {
// namespace sph
// {
// template<> std::string get_SE_name<linear>()
// {return "L";}
// template<> std::string get_SE_name<tait>()
// {return "T";}

// namespace detail 
// {

// void export_StateEquations(pybind11::module& m)
//     {
//     pybind11::class_<StateEquation<tait>, std::shared_ptr<StateEquation<tait>>(m, "Tait")
//         .def("Pressure", &StateEquation<tait>::Pressure)
//         .def("Density", &StateEquation<tait>::Density)
//         .def("setParams", &StateEquation<tait>::setParams);
//     pybind11::class_<StateEquation<linear>, std::shared_ptr<StateEquation<linear>>(m, "Linear")
//         .def("Pressure", &StateEquation<linear>::Pressure)
//         .def("Density", &StateEquation<linear>::Density)
//         .def("setParams", &StateEquation<linear>::setParams);
//     }

// }    // end namespace detail 
// }    // end namespace sph 
// }    // end namespace hoomd 
