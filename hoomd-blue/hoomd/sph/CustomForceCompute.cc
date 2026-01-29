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


#include "CustomForceCompute.h"
#include "hoomd/PythonLocalDataAccess.h"

namespace py = pybind11;

using namespace std;

/*! \file CustomForceCompute.cc
    \brief Contains code for the CustomForceCompute class
*/

namespace hoomd
    {
namespace sph
    {
/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
 */
CustomForceCompute::CustomForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                                       pybind11::object py_setForces
                                       // , bool aniso
                                       )
    : ForceCompute(sysdef) //, m_aniso(aniso)
    {
    m_exec_conf->msg->notice(5) << "Constructing ConstForceCompute" << endl;
    m_setForces = py_setForces;
    m_buffers_writeable = true;
    }

CustomForceCompute::~CustomForceCompute()
    {
    m_exec_conf->msg->notice(5) << "Destroying ConstForceCompute" << endl;
    }

/*! This function calls the python set_forces method.
    \param timestep Current timestep
*/
void CustomForceCompute::computeForces(uint64_t timestep)
    {
        // zero necessary arrays
        {
        ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::overwrite);
        m_force.zeroFill();
        }
    // execute python callback to update the forces, if present
    m_setForces(timestep);
    }

namespace detail
    {
void export_CustomForceCompute(py::module& m)
    {
    py::class_<CustomForceCompute, ForceCompute, std::shared_ptr<CustomForceCompute>>(
        m,
        "CustomForceCompute")
        .def(py::init<std::shared_ptr<SystemDefinition>, pybind11::object>());
    }

    } // end namespace detail
    } // end namespace sph
    } // end namespace hoomd
