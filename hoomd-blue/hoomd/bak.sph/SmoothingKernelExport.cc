// Copyright (c) 2009-2022 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: David Krach

#include "SmoothingKernel.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

using namespace std;

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param nlist Neighbor list instance. Must not be NULL.
    \post SmoothingKernel base class is constructed.
*/

namespace hoomd
{
namespace sph
{
template<> std::string get_kernel_name<wendlandc2>()
{return "WC2";}
template<> std::string get_kernel_name<wendlandc4>()
{return "WC4";}
template<> std::string get_kernel_name<wendlandc6>()
{return "WC6";}
template<> std::string get_kernel_name<quintic>()
{return "Q";}
template<> std::string get_kernel_name<cubicspline>()
{return "CS";}

namespace detail
{
 void export_WendlandC2(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc2>, std::shared_ptr<SmoothingKernel<wendlandc2>>>(m, "WendlandC2")
         .def("getKernelKappa", &SmoothingKernel<wendlandc2>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc2>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc2>::dwijdr)
         ;
     }
 void export_WendlandC4(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc4>, std::shared_ptr<SmoothingKernel<wendlandc4>>>(m, "WendlandC4")
         .def("getKernelKappa", &SmoothingKernel<wendlandc4>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc4>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc4>::dwijdr)
         ;
     }
 void export_WendlandC6(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<wendlandc6>, std::shared_ptr<SmoothingKernel<wendlandc6>>>(m, "WendlandC6")
         .def("getKernelKappa", &SmoothingKernel<wendlandc6>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<wendlandc6>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<wendlandc6>::dwijdr)
         ;
     }
 void export_Quintic(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<quintic>, std::shared_ptr<SmoothingKernel<quintic>>>(m, "Quintic")
         .def("getKernelKappa", &SmoothingKernel<quintic>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<quintic>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<quintic>::dwijdr)
         ;
     }
 void export_CubicSpline(pybind11::module& m)
     {
     pybind11::class_<SmoothingKernel<cubicspline>, std::shared_ptr<SmoothingKernel<cubicspline>>>(m, "CubicSpline")
         .def("getKernelKappa", &SmoothingKernel<cubicspline>::getKernelKappa)
         //.def("setNeighborList", &SmoothingKernel<wendlandc2>::setNeighborList)
         .def("EvalKernel", &SmoothingKernel<cubicspline>::wij)
         .def("EvalKernelDerivative", &SmoothingKernel<cubicspline>::dwijdr)
         ;
     }
 } // end namespace detail 
 } // end namespace sph
 } // end namespace hoomd 