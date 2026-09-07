#!/usr/bin/env python3
"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved. BSD-3-Clause, see repository LICENSE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

Bingham slug mobilization benchmark (06) — postprocessing.

Reads mobilization_results.txt (tau_y g v_mid v_end u_creep u_newt), and
classifies each point by the late-time growth ratio v_end/v_mid:
an ARRESTED slug has plateaued at the Papanastasiou creep floor (ratio ~1),
a MOBILIZED slug is still accelerating on the mu_p viscous time scale
(ratio >> 1). The measured critical body force is bracketed between the last
arrested and first mobilized point and compared against the analytic
quasi-static criterion g_crit = 2 tau_y l / (rho r L) (the periodic pressure
field concentrates the column drive onto the slug in the arrested limit).
"""

import numpy as np

dx = 5.0e-4
r_in, rho0 = 8 * dx, 1000.0
l_s, L = 16 * dx, 48 * dx
GROWTH_THRESHOLD = 1.25

data = np.loadtxt("mobilization_results.txt", ndmin=2)
print(f"{'tau_y':>6} {'g':>8} {'g/gcrit':>8} {'v_mid':>11} {'v_end':>11} "
      f"{'growth':>7} {'creep':>10}  state")
for tau_y in np.unique(data[:, 0]):
    rows = data[data[:, 0] == tau_y]
    rows = rows[rows[:, 1].argsort()]
    g_crit = 2.0 * tau_y * l_s / (rho0 * r_in * L)
    g_arr, g_mob = None, None
    for tau, g, vm, ve, ucreep, unewt in rows:
        growth = ve / vm if vm > 0 else float("inf")
        # mobilized if still accelerating late OR clearly above the creep
        # floor (calibrated: arrested plateaus at ~1.2x creep, mobilized >2.6x)
        mobilized = (growth > GROWTH_THRESHOLD) or (ve > 2.0 * ucreep)
        state = "MOBILIZED" if mobilized else "arrested"
        print(f"{tau:6.2f} {g:8.4f} {g/g_crit:8.2f} {vm:11.4e} {ve:11.4e} "
              f"{growth:7.2f} {ucreep:10.3e}  {state}")
        if not mobilized:
            g_arr = g
        elif g_mob is None:
            g_mob = g
    if g_arr is not None and g_mob is not None and g_mob > g_arr:
        g_meas = 0.5 * (g_arr + g_mob)
        err = (g_meas - g_crit) / g_crit * 100.0
        print(f"  tau_y={tau_y:g}: measured g_crit in [{g_arr:.3f}, {g_mob:.3f}] "
              f"-> {g_meas:.3f} m/s²  |  analytic 2·tau_y·l/(rho·r·L) = {g_crit:.3f} m/s²  "
              f"({err:+.0f}%)")
    else:
        print(f"  tau_y={tau_y:g}: transition not bracketed by the sweep")
    print()
