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


def set_max_sl(sim, device, snapshot, model):
    """
    

    Parameters
    ----------
    sim : hoomd simulation type 
        DESCRIPTION.
    device : hoomd device type
        DESCRIPTION.
    snapshot : hoomd snapshot
        DESCRIPTION.
    model : sph model
        DESCRIPTION.
    Returns
    -------
    maximum_smoothing_length

    """

    maximum_smoothing_length = 0.0
    # Call get_snapshot on all ranks.
    snapshot = sim.state.get_snapshot()
    # Access particle data on rank 0 only.
    if snapshot.communicator.rank == 0:
        maximum_smoothing_length = np.max(snapshot.particles.slength)

    maximum_smoothing_length = device.communicator.bcast_double(maximum_smoothing_length)
    model.max_sl = maximum_smoothing_length

    return maximum_smoothing_length

