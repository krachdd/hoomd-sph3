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

Creates the initial GSD file for the viscoplastic invasion benchmark (05).

BENCHMARK DESCRIPTION
---------------------
Two-phase immiscible displacement in a periodic simple-cubic (SC) sphere
packing. A Newtonian invader (phase A) initially occupies the pore space in a
slab at the $-x$ end of the box; the defender (phase B — Newtonian baseline or
Bingham at run time) saturates the remaining pore space:

    lattice constant  a = a_vox * dx
    sphere radius     R = R_frac * a           (default R_frac = 0.45)
    pore-channel radius  r_c = (sqrt(2)/2 - R_frac) * a
    porosity          phi = 1 - (4 pi / 3) R_frac^3   (~0.62 at default)

Both phases are density-matched (rho_0 = 1000 kg/m^3) so that any difference
between a Newtonian-defender and a Bingham-defender run is pure rheology.
Types: A = 0 (invader), B = 1 (defender), S = 2 (solid grains).

Usage:
    python3 create_packing.py [cells] [a_vox] [out.gsd]
      cells  : SC unit cells per direction        (default 2)
      a_vox  : lattice constant in particle units (default 16)
      out    : output file name                   (default packing_init.gsd)
"""

import sys
import numpy as np
import gsd.hoomd

import os
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph

cells  = int(sys.argv[1]) if len(sys.argv) > 1 else 2
a_vox  = int(sys.argv[2]) if len(sys.argv) > 2 else 16
out    = sys.argv[3] if len(sys.argv) > 3 else "packing_init.gsd"

# ─── Geometry parameters ─────────────────────────────────────────────────────
dx      = 2.0e-4                 # particle spacing            [m]
R_frac  = 0.45                   # sphere radius / a
rho0    = 1000.0                 # both phases density-matched [kg/m³]
slab_x  = 0.375                  # invader slab fraction of Lx

n       = cells * a_vox          # particles per direction
L       = n * dx
a       = a_vox * dx
R       = R_frac * a

kernel  = "WendlandC4"
slength = sph.kernel.OptimalH[kernel] * dx

# particle lattice, cell-centered
xs = np.linspace(-L/2 + dx/2, L/2 - dx/2, n)
xg, yg, zg = np.meshgrid(xs, xs, xs, indexing="ij")
pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
N = pos.shape[0]

# SC sphere centers (periodic): at cell corners
centers = np.array([[i*a - L/2, j*a - L/2, k*a - L/2]
                    for i in range(cells+1)
                    for j in range(cells+1)
                    for k in range(cells+1)])

# minimum periodic distance to any center
solid = np.zeros(N, dtype=bool)
for c in centers:
    d = pos - c
    d -= L * np.round(d / L)
    solid |= (d**2).sum(axis=1) < R*R

typeid = np.ones(N, dtype=np.int32)          # default: B (defender)
typeid[solid] = 2                             # S
invader = (~solid) & (pos[:, 0] < -L/2 + slab_x*L)
typeid[invader] = 0                           # A

nA = int((typeid == 0).sum()); nB = int((typeid == 1).sum()); nS = int((typeid == 2).sum())
phi = 1.0 - nS / N
r_pore = (np.sqrt(2.0)/2.0 - R_frac) * a

frame = gsd.hoomd.Frame()
frame.configuration.box = [L, L, L, 0, 0, 0]
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

print(f"Written {out}: N={N}  A(invader)={nA}  B(defender)={nB}  S={nS}")
print(f"  L={L*1e3:.2f} mm  a={a*1e3:.2f} mm  R={R*1e3:.2f} mm  porosity={phi:.3f}")
print(f"  pore-channel radius r_c ~ {r_pore*1e3:.3f} mm  ({r_pore/dx:.1f} dx)")
