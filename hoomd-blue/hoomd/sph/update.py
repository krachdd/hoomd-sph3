"""----------------------------------------------------------
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved. BSD-3-Clause license, see repository LICENSE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------

Adaptive time step control for SPH simulations.

Example::

    integrator = hoomd.sph.Integrator(dt=dt0)
    ...
    adapt = hoomd.sph.update.AdaptiveTimestep(model=model, courant=0.25)
    sim.operations.updaters.append(adapt.as_updater(
        trigger=hoomd.trigger.Periodic(50)))
"""

import math

import hoomd


class AdaptiveTimestep(hoomd.custom.Action):
    r"""Recompute the integrator time step from the current flow state.

    Every time the wrapping updater triggers, the new time step is

    .. math::

        \Delta t = C \cdot \min\!\left(
            \frac{h}{c + v_\mathrm{max}},\;
            \frac{\rho_0 h^2}{8 \mu},\;
            \frac{1}{4}\sqrt{h / g}\right)

    i.e. the acoustic CFL condition (with the maximum fluid speed of the last
    force computation, MPI-reduced), the viscous diffusion condition, and the
    body-force condition (skipped when no body force is set). For GDGD models
    with a diffusing scalar (``kappa_s > 0``) the scalar-diffusion condition
    :math:`h^2/(8\kappa_s)` is included as well. The result is
    clamped to ``[dt_min, dt_max]`` and the growth per update is limited by
    ``max_growth`` to avoid oscillations; shrinking is applied immediately.

    Args:
        model: an attached ``hoomd.sph.sphmodel`` model (single-phase family).
            Supplies ``getMaxVelocity()`` plus the material parameters unless
            they are overridden explicitly.
        courant (float): safety factor :math:`C` applied to the minimum.
        dt_min (float): lower clamp; defaults to no lower bound.
        dt_max (float): upper clamp; defaults to no upper bound.
        max_growth (float): maximum factor by which dt may grow per update.
        c (float): speed of sound override; default ``model.eos.SpeedOfSound``.
        h (float): smoothing length override; default ``model.max_sl``.
        mu (float): dynamic viscosity override; default ``model.mu``.
        rho0 (float): rest density override; default ``model.eos.RestDensity``.
    """

    def __init__(self, model, courant=0.25, dt_min=0.0, dt_max=float("inf"),
                 max_growth=1.1, c=None, h=None, mu=None, rho0=None):
        super().__init__()
        if courant <= 0.0:
            raise ValueError("courant must be > 0")
        if max_growth <= 1.0:
            raise ValueError("max_growth must be > 1")
        self._model = model
        self._courant = float(courant)
        self._dt_min = float(dt_min)
        self._dt_max = float(dt_max)
        self._max_growth = float(max_growth)
        self._c = c
        self._h = h
        self._mu = mu
        self._rho0 = rho0

    def as_updater(self, trigger):
        """Wrap this action in a ``hoomd.update.CustomUpdater``."""
        return hoomd.update.CustomUpdater(trigger=trigger, action=self)

    def _param(self, override, getter):
        return float(override) if override is not None else float(getter())

    def act(self, timestep):
        model = self._model
        if getattr(model, "_cpp_obj", None) is None:
            return  # model not attached yet

        c    = self._param(self._c,    lambda: model.eos.SpeedOfSound)
        h    = self._param(self._h,    lambda: model.max_sl)
        mu   = self._param(self._mu,   lambda: model.mu)
        rho0 = self._param(self._rho0, lambda: model.eos.RestDensity)
        if c <= 0.0 or h <= 0.0:
            raise RuntimeError(
                "AdaptiveTimestep: speed of sound and smoothing length must be "
                "positive (did you call compute_speedofsound()?)")

        # Maximum fluid speed of the last force computation, MPI-reduced.
        vmax = float(model._cpp_obj.getMaxVelocity())

        candidates = [h / (c + vmax)]
        if mu > 0.0:
            candidates.append(rho0 * h * h / (8.0 * mu))
        gmag = model.get_GMAG() if hasattr(model, "get_GMAG") else 0.0
        if gmag > 0.0:
            candidates.append(0.25 * math.sqrt(h / gmag))

        # Scalar diffusion condition (GDGD models transporting a scalar
        # with diffusivity kappa_s): dt <= h^2 / (8 kappa_s).
        try:
            kappa_s = float(model._param_dict["kappa_s"])
        except (AttributeError, KeyError):
            kappa_s = 0.0
        if kappa_s > 0.0:
            candidates.append(h * h / (8.0 * kappa_s))

        dt_new = self._courant * min(candidates)

        integrator = self._state._simulation.operations.integrator
        dt_old = integrator.dt
        # limit growth, allow immediate shrink
        dt_new = min(dt_new, dt_old * self._max_growth)
        dt_new = min(max(dt_new, self._dt_min), self._dt_max)

        if dt_new != dt_old:
            integrator.dt = dt_new


__all__ = ["AdaptiveTimestep"]
