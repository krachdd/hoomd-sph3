#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

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
    EOS.set_params(rho0 ,0.05)

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
    with sim.state.cpu_local_snapshot as data:
        for i in range(len(data.particles.position)):
            if data.particles.typeid[i] == 1 and data.particles.energy[i] == 1:
                tags.append(data.particles.tag[i])
                # print(f'Rank: {device.communicator.rank} -> Delete Particle {data.particles.tag[i]}')
                deleted += 1

    for t in tags:
        # print(f'Rank: {device.communicator.rank} --> Remove particle {t} of {deleted}')
        sim.state.removeParticle(t)

    # if device.communicator.rank == 0:
    print(f'Rank {device.communicator.rank}: {deleted} unnecessary solid particles deleted.')

    return sim, deleted

