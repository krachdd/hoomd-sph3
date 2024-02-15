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


def print_decomp_info(sim, device):
    """
    

    Parameters
    ----------
    sim : hoomd simulation type 
        DESCRIPTION.
    device : hoomd device type
        DESCRIPTION.
    
    Returns
    -------
    Nothing.

    """

    # Print the domain decomposition.
    domain_decomposition = sim.state.domain_decomposition
    if device.communicator.rank == 0:
        print(f'Domain Decomposition: {domain_decomposition}')

    # Print the location of the split planes.
    split_fractions = sim.state.domain_decomposition_split_fractions
    if device.communicator.rank == 0:
        print(f'Locations of SplitPlanes: {split_fractions}')

    # Print the number of particles on each rank.
    with sim.state.cpu_local_snapshot as snap:
        N = len(snap.particles.position)
        print(f'{N} particles on rank {device.communicator.rank}')