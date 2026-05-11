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

SPH simulation helper utilities.

This module provides helper functions used across SPH run scripts:
  - Setting the maximum smoothing length on the SPH model (required before
    the neighbor-list cutoff can be determined).
  - Computing and adaptively updating the speed-of-sound parameter c0 for
    both single-phase (SPF) and two-phase (TPF) flow models, enforcing the
    weak-compressibility condition :math:`\mathrm{Ma} = u_\mathrm{ref}/c_0 \leq 0.01`.
"""

# --- HEADER -------------------------------------------------
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd
# ------------------------------------------------------------


def set_max_sl(sim, device, model):
    """Read the maximum smoothing length from the simulation state and set it
    on the SPH model.

    The maximum smoothing length is needed to determine the neighbor-list
    cutoff radius (r_cut = kappa * max_sl).  This must be called after the
    simulation has been created but before the first run step.

    Parameters
    ----------
    sim : hoomd.Simulation
        Active HOOMD simulation object.
    device : hoomd.device.Device
        HOOMD device (CPU or GPU); used for MPI broadcast.
    model : hoomd.sph.sphmodel.SPHModel
        SPH model instance whose ``max_sl`` attribute will be updated.

    Returns
    -------
    float
        The maximum smoothing length found across all particles (broadcast to
        all MPI ranks).
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
    """Estimate the reference speed of sound from a body-force scale.

    Uses the relation :math:`c_0 = \text{cfactor} \cdot \sqrt{g \cdot l_\mathrm{ref}}`,
    which ensures that pressure fluctuations driven by gravity/body-force remain
    small compared to the background pressure.

    Parameters
    ----------
    lref : float
        Reference length scale (e.g. channel half-width) [m].
    bforce : float
        Body-force magnitude (e.g. gravitational acceleration g) [:math:`\mathrm{m/s^2}`].
    cfactor : float
        Multiplier, typically 10, so that :math:`\mathrm{Ma} \leq 0.1` (then clamped to 0.01).

    Returns
    -------
    float
        Estimated speed of sound [m/s].
    """
    return cfactor * np.sqrt(bforce * lref)


def get_c0_umax(uref, cfactor):
    """Estimate the reference speed of sound from a velocity scale.

    Uses the relation :math:`c_0 = \text{cfactor} \cdot u_\mathrm{ref}`, which
    enforces weak compressibility:
    :math:`\mathrm{Ma} = u_\mathrm{ref} / c_0 = 1 / \text{cfactor}`.

    Parameters
    ----------
    uref : float
        Reference (maximum) flow velocity [m/s].
    cfactor : float
        Multiplier, typically 10, yielding :math:`\mathrm{Ma} = 0.1` (then clamped to 0.01).

    Returns
    -------
    float
        Estimated speed of sound [m/s].
    """
    return cfactor * uref


def update_min_c0(device, model, c, mode='uref', lref=0.0, uref=0.0, bforce=0.0, cfactor=10.0):
    """Adaptively increase the speed of sound for a single-phase SPH model.

    Computes a target :math:`c_0` from the flow scales, clamps it to satisfy
    :math:`\mathrm{Ma} = u_\mathrm{ref}/c_0 \leq 0.01`, and calls
    ``model.set_speedofsound(c0)`` only if
    the new value exceeds the current value ``c``.  Does nothing if the model
    is already fast enough.

    Parameters
    ----------
    device : hoomd.device.Device
        HOOMD device; used to restrict console output to rank 0.
    model : hoomd.sph.sphmodel.SPHModel
        Single-phase SPH model instance.
    c : float
        Current speed of sound stored in the model [m/s].
    mode : {'uref', 'bforce', 'both'}
        Strategy for estimating c0:
        - ``'uref'``   — c0 from reference velocity only (requires ``uref``).
        - ``'bforce'`` — c0 from body-force scale (requires ``bforce``,
          ``lref``, and ``uref`` for the Mach-number check).
        - ``'both'``   — maximum of the two estimates (requires all three).
    lref : float, optional
        Reference length scale [m].  Required for ``'bforce'`` and ``'both'``.
    uref : float, optional
        Reference velocity [m/s].  Required for all modes.
    bforce : float, optional
        Body-force magnitude [:math:`\mathrm{m/s^2}`].  Required for ``'bforce'`` and ``'both'``.
    cfactor : float, optional
        Speed-of-sound multiplier, default 10.

    Returns
    -------
    None
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

    # Clamp to Ma <= 0.01 for weak compressibility.
    Ma = uref / c0
    if Ma > 0.01:
        c0 *= 0.01 / Ma
        Ma = uref / c0
    if c > c0:
        if device.communicator.rank == 0:
            print(f'c0 not updated, Ma = {uref/c}')
    else:
        model.set_speedofsound(c0)
        if device.communicator.rank == 0:
            print(f'Increase Speed of Sound: {model.get_speedofsound()}, Ma = {Ma}')


def get_dt_cfl(h, c0, umax=0.0, mu=0.0, rho0=1.0, bforce=0.0, CFL=0.25):
    """Compute the maximum stable time step from CFL, viscous, and body-force conditions.

    Evaluates three independent time-step restrictions and returns the minimum:

    1. **CFL (acoustic)**:
       :math:`\\Delta t_\\mathrm{CFL} = \\mathrm{CFL} \\cdot h / (c_0 + u_\\mathrm{max})`

    2. **Viscous** (only when ``mu > 0``):
       :math:`\\Delta t_\\mathrm{visc} = 0.125 \\cdot h^2 / \\nu`,
       where :math:`\\nu = \\mu / \\rho_0`.

    3. **Body force** (only when ``bforce > 0``):
       :math:`\\Delta t_\\mathrm{body} = 0.5 \\cdot \\sqrt{h / |g|}`.

    These criteria match the expressions used in DualSPHysics and are standard
    in weakly-compressible SPH (Monaghan 1992, Morris et al. 1997).

    Parameters
    ----------
    h : float
        Smoothing length [m].
    c0 : float
        Speed of sound [m/s].
    umax : float, optional
        Maximum flow velocity [m/s].  Default 0.
    mu : float, optional
        Dynamic viscosity [Pa·s].  Set to 0 to skip the viscous criterion.
    rho0 : float, optional
        Reference density [kg/m³].  Used to compute kinematic viscosity.
    bforce : float, optional
        Body-force magnitude [m/s²].  Set to 0 to skip the body-force criterion.
    CFL : float, optional
        CFL number, typically 0.1–0.25.  Default 0.25.

    Returns
    -------
    dict
        ``{'dt': float, 'dt_cfl': float, 'dt_visc': float or None,
           'dt_body': float or None, 'limiting': str}``
        where ``'limiting'`` identifies which criterion is active.

    Examples
    --------
    >>> info = get_dt_cfl(h=0.01, c0=10.0, umax=0.5, mu=1e-3, rho0=1000.0, bforce=9.81)
    >>> print(f"Use dt = {info['dt']:.4e}  (limited by {info['limiting']})")
    """
    if c0 <= 0.0:
        raise ValueError("c0 must be positive.")
    if h <= 0.0:
        raise ValueError("h must be positive.")

    dt_cfl  = CFL * h / (c0 + umax)
    dt_visc = None
    dt_body = None

    candidates = {'CFL': dt_cfl}

    if mu > 0.0:
        nu = mu / rho0
        dt_visc = 0.125 * h * h / nu
        candidates['viscous'] = dt_visc

    if bforce > 0.0:
        dt_body = 0.5 * np.sqrt(h / bforce)
        candidates['body_force'] = dt_body

    limiting = min(candidates, key=candidates.get)
    dt = candidates[limiting]

    return {
        'dt':       dt,
        'dt_cfl':   dt_cfl,
        'dt_visc':  dt_visc,
        'dt_body':  dt_body,
        'limiting': limiting,
    }


def update_min_c0_tpf(device, model, c1, c2, mode='plain', lref=0.0, uref=0.0, bforce=0.0, cfactor=10.0):
    """Adaptively increase the speed of sound for a two-phase SPH model.

    Computes target speeds of sound c01 and c02 for the two fluid phases and
    calls ``model.set_speedofsound(c01, c02)`` only if at least one current
    value is below its target.  Currently only ``mode='plain'`` is
    implemented; velocity- and body-force-based modes raise
    ``NotImplementedError``.

    Parameters
    ----------
    device : hoomd.device.Device
        HOOMD device; used to restrict console output to rank 0.
    model : hoomd.sph.sphmodel.TwoPhaseFlow
        Two-phase SPH model instance.
    c1 : float
        Current speed of sound for phase 1 [m/s].
    c2 : float
        Current speed of sound for phase 2 [m/s].
    mode : {'plain', 'uref', 'bforce', 'both'}
        Strategy for estimating c0.  Only ``'plain'`` is currently
        implemented: c0i = ci * cfactor.
    lref : float, optional
        Reference length scale [m].  Reserved for future modes.
    uref : float, optional
        Reference velocity [m/s].  Used only for the printed Mach number.
    bforce : float, optional
        Body-force magnitude [:math:`\mathrm{m/s^2}`].  Reserved for future modes.
    cfactor : float, optional
        Speed-of-sound multiplier, default 10.

    Returns
    -------
    None
    """

    if mode == 'uref':
        raise NotImplementedError
    elif mode == 'bforce':
        raise NotImplementedError
    elif mode == 'both':
        raise NotImplementedError
    elif mode == 'plain':
        c01 = c1 * cfactor
        if c01 <= 0.0:
            raise ValueError('c0 must not be smaller or equal to 0.')
        c02 = c2 * cfactor
        if c02 <= 0.0:
            raise ValueError('c0 must not be smaller or equal to 0.')
    else:
        raise ValueError('Give correct mode')

    Ma1 = uref / c01
    Ma2 = uref / c02
    if c1 > c01 or c2 > c02:
        if device.communicator.rank == 0:
            print(f'c0 not updated, Ma1 = {Ma1}, Ma2 = {Ma2}')
    else:
        model.set_speedofsound(c01, c02)
        if device.communicator.rank == 0:
            print(f'Increase Speed of Sound: {model.get_speedofsound()}, Ma1 = {Ma1}, Ma2 = {Ma2}')







