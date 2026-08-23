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

Creates the pseudo-3D piston-driven invasion domain (WP2 variant).

DOMAIN LAYOUT (x = flow direction, all boundaries periodic)
-----------------------------------------------------------
    | piston P -> | reservoir A (invader) | cylinder matrix (B-saturated) |
    | reservoir B (defender) ... wraps around to the piston's rear |

  * Pseudo-3D: thin periodic slab in z (default 8 dx > 2 r_cut), grains are
    CYLINDERS with axis z -> a quasi-2D micromodel (Lenormand-type geometry).
  * The piston is a solid slab of its own type 'P' with prescribed velocity
    U_p x^: the Adami no-slip BC drives the fluid; the run script advances the
    piston positions kinematically. Velocity-controlled displacement: Ca is
    imposed exactly, and the piston drag force (ComputeSolidProperties on 'P')
    measures the driving pressure — the observable for threshold-gradient
    studies.
  * Valid until the piston has traveled ~ one reservoir length.

Types: A = 0 (invader), B = 1 (defender), S = 2 (matrix grains), P = 3 (piston).

Usage:
    python3 create_piston_domain.py [--U_p 0.02] [--res 8] [--out FILE]
                                    [--resA_mm 12.0] [--resB_mm 12.6] [--no-piston]
      --no-piston : omit the piston slab (body-force-driven variant); the domain
                    is then reservoir A | matrix | reservoir B, periodic, with
                    the second A|B interface across the wrap.
      U_piston : piston velocity [m/s] written into the P particles (default 0.02)
      res      : throat resolution in particles (default 8; physical size fixed)
"""

import sys, os, argparse
import numpy as np
import gsd.hoomd

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import hoomd
from hoomd import sph

ap = argparse.ArgumentParser()
ap.add_argument("--U_p", type=float, default=0.02, help="piston velocity [m/s]")
ap.add_argument("--res", type=int, default=8, help="throat width in particles")
ap.add_argument("--out", type=str, default=None)
ap.add_argument("--resA_mm", type=float, default=12.0, help="invader reservoir length [mm]")
ap.add_argument("--resB_mm", type=float, default=12.6, help="defender reservoir length [mm]")
ap.add_argument("--no-piston", dest="piston", action="store_false")
args = ap.parse_args()

U_p, res = args.U_p, args.res
out = args.out or (f"piston_domain_res{res}_init.gsd" if args.piston
                   else f"bodyforce_domain_res{res}_init.gsd")

# ─── Geometry parameters (units of dx) ──────────────────────────────────────
# The PHYSICAL domain is fixed; `res` sets the throat resolution in particles
# (dx scales as 8/res). nz stays at 8 dx: the pseudo-3D thickness only needs
# to exceed the minimum-image bound, which is resolution-independent in dx.
scale     = res / 8.0
dx        = 2.0e-4 * 8.0 / res   # particle spacing                [m]
nz        = 8                    # pseudo-3D slab thickness        [dx]
ny        = round(96 * scale)    # cross-section height            [dx]
n_piston  = max(6, round(6 * scale)) if args.piston else 0  # piston  [dx]
n_resA    = round(args.resA_mm * 1e-3 / dx)  # invader reservoir length [dx]
n_matrix  = round(96 * scale)    # porous matrix length            [dx]
n_resB    = round(args.resB_mm * 1e-3 / dx)  # defender reservoir length [dx]
R_vox     = res                  # cylinder radius                  [dx]
a_vox     = 3 * res              # cylinder lattice constant        [dx]
rho0      = 1000.0               # both phases density-matched      [kg/m³]

nx = n_piston + n_resA + n_matrix + n_resB
Lx, Ly, Lz = nx * dx, ny * dx, nz * dx

kernel  = "WendlandC4"
slength = sph.kernel.OptimalH[kernel] * dx

xs = np.linspace(-Lx/2 + dx/2, Lx/2 - dx/2, nx)
ys = np.linspace(-Ly/2 + dx/2, Ly/2 - dx/2, ny)
zs = np.linspace(-Lz/2 + dx/2, Lz/2 - dx/2, nz)
xg, yg, zg = np.meshgrid(xs, ys, zs, indexing="ij")
pos = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
N = pos.shape[0]

# section boundaries in x
x_piston0 = -Lx/2
x_piston1 = x_piston0 + n_piston * dx
x_resA1   = x_piston1 + n_resA * dx
x_mat1    = x_resA1 + n_matrix * dx

x, y = pos[:, 0], pos[:, 1]
typeid = np.full(N, 1, dtype=np.int32)               # default: B (defender)
typeid[(x >= x_piston1) & (x < x_resA1)] = 0          # A (invader reservoir)
typeid[x < x_piston1] = 3                             # P (piston)

# staggered cylinder lattice inside the matrix section
n_cols = n_matrix // a_vox
n_rows = ny // a_vox
R = R_vox * dx
a = a_vox * dx
centers = []
for c in range(n_cols):
    cx = x_resA1 + (c + 0.5) * a
    y_off = 0.5 * a if (c % 2 == 1) else 0.0
    for r in range(n_rows):
        cy = -Ly/2 + (r + 0.5) * a + y_off
        centers.append((cx, cy))

in_matrix = (x >= x_resA1) & (x < x_mat1)
for cx, cy in centers:
    dxx = x - cx
    dyy = y - cy
    dyy -= Ly * np.round(dyy / Ly)                    # periodic y
    grain = in_matrix & (dxx*dxx + dyy*dyy < R*R)
    typeid[grain] = 2                                 # S

nA = int((typeid == 0).sum()); nB = int((typeid == 1).sum())
nS = int((typeid == 2).sum()); nP = int((typeid == 3).sum())
n_mat_total = int(in_matrix.sum())
phi = 1.0 - nS / max(n_mat_total, 1)
throat = a - 2*R

vel = np.zeros((N, 3), dtype=np.float32)
vel[typeid == 3, 0] = U_p                             # prescribed piston velocity

frame = gsd.hoomd.Frame()
frame.configuration.box = [Lx, Ly, Lz, 0, 0, 0]
frame.particles.N = N
frame.particles.types = ["A", "B", "S", "P"]
frame.particles.typeid = typeid
frame.particles.position = pos.astype(np.float32)
frame.particles.velocity = vel
frame.particles.mass = np.full(N, rho0 * dx**3, dtype=np.float32)
frame.particles.slength = np.full(N, slength, dtype=np.float32)
frame.particles.density = np.full(N, rho0, dtype=np.float32)

with gsd.hoomd.open(name=out, mode="w") as f:
    f.append(frame)

print(f"Written {out}: N={N}")
print(f"  A(invader)={nA}  B(defender)={nB}  S(grains)={nS}  P(piston)={nP}")
print(f"  box {Lx*1e3:.1f} x {Ly*1e3:.1f} x {Lz*1e3:.1f} mm  (nx={nx}, ny={ny}, nz={nz} dx)")
print(f"  matrix: {len(centers)} cylinders R={R*1e3:.1f} mm, a={a*1e3:.1f} mm, "
      f"throat={throat*1e3:.1f} mm ({throat/dx:.0f} dx), porosity={phi:.3f}")
print(f"  piston: U_p={U_p} m/s, max travel ~ reservoir A = {n_resA*dx*1e3:.1f} mm "
      f"-> t_max ~ {n_resA*dx/U_p:.3f} s")
