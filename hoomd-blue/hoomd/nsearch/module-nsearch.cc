// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace nsearch
    {
namespace detail
    {


void export_NeighborList(pybind11::module& m);
void export_NeighborListBinned(pybind11::module& m);
void export_NeighborListStencil(pybind11::module& m);
void export_NeighborListTree(pybind11::module& m);


#ifdef ENABLE_HIP


void export_NeighborListGPU(pybind11::module& m);
void export_NeighborListGPUBinned(pybind11::module& m);
void export_NeighborListGPUStencil(pybind11::module& m);
void export_NeighborListGPUTree(pybind11::module& m);

#endif
    } // namespace detail
    } // namespace nsearch
    } // namespace hoomd

using namespace hoomd;
using namespace hoomd::nsearch;
using namespace hoomd::nsearch::detail;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
    create the md python module and define the exports here.
*/
PYBIND11_MODULE(_nsearch, m)
    {

    // export_CustomForceCompute(m);
    export_NeighborList(m);
    export_NeighborListBinned(m);
    export_NeighborListStencil(m);
    export_NeighborListTree(m);

#ifdef ENABLE_HIP
    export_NeighborListGPU(m);
    export_NeighborListGPUBinned(m);
    export_NeighborListGPUStencil(m);
    export_NeighborListGPUTree(m);

#endif
    }
