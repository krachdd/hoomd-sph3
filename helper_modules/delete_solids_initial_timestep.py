#!/usr/bin/env python3

"""
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

"""

# --- HEADER -------------------------------------------------
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd
# ------------------------------------------------------------


def delete_solids(sim, device, kernel, dt, mu, DX, rho0):
    """
    

    Parameters
    ----------
    sim : hoomd simulation type 
        DESCRIPTION.
    device : hoomd device type
        DESCRIPTION.
    kernel : str
        DESCRIPTION.
    dt : float
        DESCRIPTION.
    mu : float
        DESCRIPTION.
    DX : float
        DESCRIPTION.
    rho0 : float
        DESCRIPTION.

    Returns
    -------
    hoomd simulation type.

    """

    # Some additional parameters
    densitymethod = 'CONTINUITY'


    # Kernel
    KERNEL  = str(kernel)
    H       = hoomd.sph.kernel.OptimalH[KERNEL]*DX       # m
    RCUT    = hoomd.sph.kernel.Kappa[KERNEL]*H           # m
    Kernel = hoomd.sph.kernel.Kernels[KERNEL]()
    Kappa  = Kernel.Kappa() 

    # Neighbor list
    NList = hoomd.nsearch.nlist.Cell(buffer = RCUT*0.05, rebuild_check_delay = 1, kappa = Kappa)
    
    # Setup all necessary simulation inputs
    EOS = hoomd.sph.eos.Linear()
    EOS.set_params(rho0 ,0.01)

    # Define groups/filters
    filterFLUID  = hoomd.filter.Type(['F']) # is zero
    filterSOLID  = hoomd.filter.Type(['S']) # is one
    filterAll    = hoomd.filter.All()

    # Set up SPH solver
    model = hoomd.sph.sphmodel.SinglePhaseFlow(kernel = Kernel,
                                               eos    = EOS,
                                               nlist  = NList,
                                               fluidgroup_filter = filterFLUID,
                                               solidgroup_filter = filterSOLID)

    model.mu = mu
    model.densitymethod = densitymethod
    model.gx = 0.0000001
    model.damp = 1000
    model.artificialviscosity = False 
    model.densitydiffusion = False
    model.shepardrenormanlization = False 

    # denfine integrator
    integrator = hoomd.sph.Integrator(dt=dt)
    VelocityVerlet = hoomd.sph.methods.VelocityVerletBasic(filter=filterFLUID, densitymethod = densitymethod)

    sim.operations.integrator = integrator
    integrator.methods.append(VelocityVerlet)
    integrator.forces.append(model)

    sim.run(1, write_at_start=False)

    # Identify solid particles with zero charge and delete them ( redundant solid particles )
    tags    = []
    deleted = 0

    # if device.communicator.rank == 0:
    #     data = sim.state.get_snapshot()
    #     for i in range(len(data.particles.position)):
    #         if data.particles.typeid[i] == 1 and data.particles.energy[i] == 1:
    #             tags.append(data.particles.tag[i])
    #             print(f'Rank: {device.communicator.rank} -> Delete Particle {data.particles.tag[i]}')
    #             deleted += 1

    with sim.state.cpu_local_snapshot as data:
        for i in range(len(data.particles.position)):
            # print(data.particles.mass[i])
            if data.particles.typeid[i] == 1 and data.particles.mass[i] == -999:
                tags.append(data.particles.tag[i])
                # print(f'Rank: {device.communicator.rank} -> Delete Particle {data.particles.tag[i]}')
                deleted += 1

    # if device.communicator.rank == 0:
    for t in tags:
        # print(f'Rank: {device.communicator.rank} --> Remove particle {t} of {deleted}')
        sim.state.removeParticle(t)

    device.communicator.barrier_all()

    print(f'Rank {device.communicator.rank}: {deleted} unnecessary solid particles deleted.')


    return sim, deleted

