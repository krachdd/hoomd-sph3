#!/usr/bin/env python3
"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

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
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

Creates the initial GSD file for the Bingham slug mobilization benchmark (06).

BENCHMARK DESCRIPTION
---------------------
A cylindrical tube (axis x, periodic) contains a Bingham slug (phase B) of
length l embedded in a Newtonian carrier fluid (phase A). A body force g_x
drives the fluid column. Quasi-static force balance on the column gives the
analytic mobilization criterion

    rho * g_x * L * pi r^2  >  2 pi r * l * tau_y
    =>  g_crit = 2 tau_y l / (rho r L)

In the quasi-static (arrested) limit the Newtonian sections carry no gradient,
so the periodic pressure field concentrates the whole column's body-force
drive onto the slug (G_slug = rho g L / l) — the slug's wall stress is
G_slug r/2 and mobilization occurs when it exceeds tau_y, giving the criterion
above. Below g_crit the slug arrests the column up to the Papanastasiou creep
floor; above it the column accelerates toward yielded flow. This is the
single-throat anchor for the yield-arrested invasion criterion (WP1 of the
viscoplastic-displacement study).

No interfacial tension is applied (sigma12 = 0), so the measured transition
isolates the yield criterion exactly.

Usage:
    python3 create_slug_tube.py [out.gsd]
"""

import sys, os
import numpy as np
import gsd.hoomd

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph

out = sys.argv[1] if len(sys.argv) > 1 else "slugtube_init.gsd"

# ─── Geometry parameters ─────────────────────────────────────────────────────
dx      = 5.0e-4                 # particle spacing              [m]
r_vox   = 8                      # inner tube radius             [dx]
n_wall  = 3                      # wall thickness                [dx]
L_vox   = 48                     # tube length (periodic x)      [dx]
l_vox   = 16                     # Bingham slug length           [dx]
rho0    = 1000.0                 # both phases density-matched   [kg/m³]

r_in  = r_vox * dx
L     = L_vox * dx
l_s   = l_vox * dx
half  = (r_vox + n_wall) * dx + dx  # transverse half-width of the box

kernel  = "WendlandC4"
slength = sph.kernel.OptimalH[kernel] * dx

ny = 2 * (r_vox + n_wall) + 2
xs = np.linspace(-L/2 + dx/2, L/2 - dx/2, L_vox)
ys = np.linspace(-half + dx/2, half - dx/2, ny)

xg, yg, zg = np.meshgrid(xs, ys, ys, indexing="ij")
pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])

# keep only particles inside the outer tube radius (walls included)
rad = np.sqrt(pos[:, 1]**2 + pos[:, 2]**2)
keep = rad < (r_vox + n_wall) * dx
pos = pos[keep]
rad = rad[keep]
N = pos.shape[0]

typeid = np.zeros(N, dtype=np.int32)                 # A (carrier)
typeid[rad >= r_in] = 2                               # S (wall)
slug = (rad < r_in) & (np.abs(pos[:, 0]) < l_s / 2)
typeid[slug] = 1                                      # B (slug)

nA = int((typeid == 0).sum()); nB = int((typeid == 1).sum()); nS = int((typeid == 2).sum())

frame = gsd.hoomd.Frame()
frame.configuration.box = [L, 2*half, 2*half, 0, 0, 0]
frame.particles.N = N
frame.particles.types = ["A", "B", "S"]
frame.particles.typeid = typeid
frame.particles.position = pos.astype(np.float32)
frame.particles.velocity = np.zeros((N, 3), dtype=np.float32)
frame.particles.mass = np.full(N, rho0 * dx**3, dtype=np.float32)
frame.particles.slength = np.full(N, slength, dtype=np.float32)
frame.particles.density = np.full(N, rho0, dtype=np.float32)

with gsd.hoomd.open(name=out, mode="w") as f:
    f.append(frame)

print(f"Written {out}: N={N}  A(carrier)={nA}  B(slug)={nB}  S(wall)={nS}")
print(f"  r={r_in*1e3:.2f} mm  L={L*1e3:.2f} mm  slug l={l_s*1e3:.2f} mm  dx={dx*1e3:.2f} mm")
print(f"  analytic: g_crit = 2 tau_y l / (rho r L) = tau_y * "
      f"{2*l_s/(rho0*r_in*L):.4f}  [m/s^2 per Pa]")
