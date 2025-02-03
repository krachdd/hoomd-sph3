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


def set_max_sl(sim, device, model):
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


def get_c0_bf(lref, bforce, cfactor):
    return cfactor * np.sqrt( bforce * lref )

def get_c0_umax(uref, cfactor):
    return cfactor * uref

def update_min_c0(device, model, c, mode = 'uref', lref = 0.0, uref = 0.0, bforce = 0.0, cfactor = 10.0):
    """
    
    Parameters
    ----------
    model : sph model
        DESCRIPTION.
    model : string
        DESCRIPTION.
    lref : 
        DESCRIPTION.
    uref : 
        DESCRIPTION.
    bforce : 
        DESCRIPTION.
    cfactor : 
        DESCRIPTION.
    Returns
    -------
    Nothing

    """

    if mode == 'uref':
        if uref <= 0.0:
            raise ValueError('Give correct uref!')
        c0 = get_c0_umax(uref, cfactor)
        if c0 <= 0.0:
            raise ValueError('c0 must not be smaller or equal to 0.')
    elif mode == 'bforce':
        if bforce <= 0.0 or lref <= 0.0 or uref <= 0.0:
            raise ValueError('Give correct bforce and lref!')
        c0 = get_c0_bf(lref, bforce, cfactor)
        if c0 <= 0.0:
            raise ValueError('c0 must not be smaller or equal to 0.') 
    elif mode == 'both':
        if bforce <= 0.0 or lref <= 0.0 or uref <= 0.0:
            raise ValueError('Give correct bforce, lref and uref!')
        c0 = np.max(np.asarray([get_c0_bf(lref, bforce, cfactor), get_c0_umax(uref, cfactor)]))
        if c0 <= 0.0:
            raise ValueError('c0 must not be smaller or equal to 0.') 
    else:
        raise ValueError('Give correct mode')

    Ma = uref/c0
    if Ma > 0.01:
        c0 *= 0.01/Ma
        Ma = uref/c0
    if c > c0:
        if device.communicator.rank == 0: 
            print(f'c0 not updated, Ma = {uref/c}')
    else:
        model.set_speedofsound(c0)
        if device.communicator.rank == 0:
            print(f'Increase Speed of Sound: {model.get_speedofsound()}, Ma = {Ma}')







