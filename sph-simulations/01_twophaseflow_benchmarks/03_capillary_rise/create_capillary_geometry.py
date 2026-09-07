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

Capillary rise — 3-D geometry generator.

Creates a square-box domain with a round (circular cross-section) capillary
tube aligned with the y-axis.  One geometry file serves all contact-angle
cases in capillary_rise_params.txt.

GEOMETRY
--------
  Outer domain : $lx \times ly \times lz$   (x, z periodic;  y bounded by solid walls)
  Tube         : inner radius R_cap,  wall thickness n_wall * dx
                 centred at (0, 0) in x-z,  spanning the full fluid height
  Initial fluid: liquid 'W' fills from y_fluid_bot to y_fluid_bot + h_res
                 gas   'N' fills above h_res  (inside and outside tube)

                          ← lx = 4 R_cap (periodic x) →
                        ┌────────────────────────────────┐  ← solid top wall
                        │ N  N  N  N  N  N  N  N  N  N  │  }
                        │ N  N  N  │   │  N  N  N  N  N  │  } gas
                        │ N  N  N  │ N │  N  N  N  N  N  │  }
                        │ N  N  N  ├───┤  N  N  N  N  N  │ ← y_liq_surf
                        │ W  W  W  │ W │  W  W  W  W  W  │  }
                        │ W  W  W  ╞═══╡  W  W  W  W  W  │  } liquid
                        │ W  W  W  │ W │  W  W  W  W  W  │  }
                        └────────────────────────────────┘  ← solid bot wall
                                   ↑ tube wall (solid)

PHYSICAL DEFAULTS  (edit constants below)
  R_cap   = 1 mm  (inner tube radius = reference length)
  n_wall  = 2     (solid layers: top/bottom walls and tube ring)
  lx = lz = 4 R_cap  (square reservoir, periodic)
  ly_fluid = 12 R_cap  (fluid column height)
  h_res   = 3 R_cap  (initial liquid height above bottom wall)
  rho01   = 1000 kg/m$^3$ (liquid)
  rho02   = 100  kg/m$^3$ (gas)

CAPILLARY RISE (Jurin's law)
  $h\_rise = 2 \sigma \cos(\theta) / (\rho_1 g R\_cap)$
  with $\sigma$=0.01 N/m, g=9.81 m/s$^2$, $\rho_1$=1000 kg/m$^3$, R_cap=1 mm:
    $\theta$=30° $\rightarrow$ h_rise $\approx$ 1.77 mm   (rise)
    $\theta$=60° $\rightarrow$ h_rise $\approx$ 1.02 mm   (rise)
    $\theta$=90° $\rightarrow$ h_rise = 0          (no rise)
    $\theta$=120° $\rightarrow$ h_rise $\approx$ -1.02 mm  (depression)

FINITE-RESERVOIR CORRECTION
  The measured height difference (meniscus above current reservoir surface)
  equals h_Jurin directly, because the Jurin formula uses the net head.
  If measuring from the INITIAL flat interface, apply the correction:
    $h\_rise\_from\_init = h\_Jurin \times A\_res / (A\_tube + A\_res)$

Usage:
    python3 create_capillary_geometry.py [num_length]
      num_length : particles across R_cap  (default: 10)

Output:
    caprise_<NX>_<NY>_<NZ>_vs_<dx:.6f>_init.gsd
"""

import sys, os
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, os.path.join(_ROOT, 'helper_modules'))

import numpy as np
import gsd.hoomd
import hoomd
from hoomd import sph

# ─── Parameters ──────────────────────────────────────────────────────────────
num_length = int(sys.argv[1]) if len(sys.argv) > 1 else 10

R_cap  = 0.001          # inner capillary radius = reference length  [m]
rho01  = 1000.0         # liquid W rest density                      [kg/m³]
rho02  = 100.0          # gas    N rest density                      [kg/m³]
dx     = R_cap / num_length  # particle spacing                      [m]
n_wall = 2              # solid layers: outer walls and tube ring
kernel = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel] * dx

# ─── Domain dimensions ───────────────────────────────────────────────────────
# x, z: periodic; y: closed by solid wall layers
NX       = 4 * num_length          # reservoir width in x:  4 R_cap
NZ       = 4 * num_length          # reservoir width in z:  4 R_cap
NY_FLUID = 12 * num_length         # fluid column height:  12 R_cap
NY_TOTAL = NY_FLUID + 2 * n_wall   # total (+ top & bottom solid layers)

lx = NX       * dx
ly = NY_TOTAL * dx
lz = NZ       * dx

# Fluid region y-bounds (first/last fluid particle centres)
y_fluid_bot = -ly / 2 + (n_wall + 0.5) * dx    # bottom of lowest fluid layer
y_fluid_top = +ly / 2 - (n_wall + 0.5) * dx    # top    of highest fluid layer

# Boundary between solid and fluid (used for type assignment)
y_wall_top_boundary = +ly / 2 - n_wall * dx     # above this: solid top wall
y_wall_bot_boundary = -ly / 2 + n_wall * dx     # below this: solid bottom wall

# Initial liquid–gas interface (flat, same inside and outside the tube)
h_res       = 3.0 * R_cap
y_liq_surf  = y_wall_bot_boundary + h_res    # absolute y-coordinate

# Tube radii
R_inner = R_cap
R_outer = R_cap + n_wall * dx               # outer edge of tube wall

# ─── Build particle lattice ──────────────────────────────────────────────────
x_arr = (np.arange(NX) + 0.5) * dx - lx / 2
y_arr = (np.arange(NY_TOTAL) + 0.5) * dx - ly / 2
z_arr = (np.arange(NZ) + 0.5) * dx - lz / 2

xg, yg, zg  = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
positions   = np.column_stack((xg.ravel(), yg.ravel(), zg.ravel())).astype(np.float32)

x_p  = positions[:, 0]
y_p  = positions[:, 1]
z_p  = positions[:, 2]
r_xz = np.sqrt(x_p**2 + z_p**2)     # radial distance from tube axis (y-axis)

N_total   = len(positions)
typeids   = np.zeros(N_total, dtype=np.int32)          # default: 'W' (0)
masses    = np.full(N_total, rho01 * dx**3, dtype=np.float32)
densities = np.full(N_total, rho01, dtype=np.float32)
slengths  = np.full(N_total, slength, dtype=np.float32)

# ─── Assign particle types ───────────────────────────────────────────────────
# 1. Solid outer walls (top and bottom in y)
is_outer_wall = (y_p < y_wall_bot_boundary) | (y_p >= y_wall_top_boundary)

# 2. Cylindrical tube wall: R_inner <= r_xz < R_outer, OPEN at both ends
# (fixed 2026-08-23): the original tube spanned the full fluid column and
# was welded to BOTH the floor and the ceiling -- a sealed cylinder whose
# trapped gas (and disconnected liquid) makes capillary rise impossible
# regardless of the wetting model. The tube now dips 1.5 R_cap into the
# reservoir with a free gap to the floor, and ends 1.5 R_cap below the top
# wall so the displaced gas can escape around the rim.
y_tube_bot = y_liq_surf - 1.5 * R_cap
y_tube_top = y_wall_top_boundary - 1.5 * R_cap
is_tube_wall  = ( (~is_outer_wall) & (r_xz >= R_inner) & (r_xz < R_outer)
                  & (y_p >= y_tube_bot) & (y_p < y_tube_top) )

# 3. Gas 'N': above the initial liquid surface, not solid
is_gas = (~is_outer_wall) & (~is_tube_wall) & (y_p >= y_liq_surf)

# Set typeids:  0 = 'W',  1 = 'N',  2 = 'S'
typeids[is_gas]                              = 1    # gas N
typeids[is_outer_wall | is_tube_wall]        = 2    # solid S

# Gas particles get gas density / mass
masses[typeids == 1]    = rho02 * dx**3
densities[typeids == 1] = rho02

# ─── Write GSD ───────────────────────────────────────────────────────────────
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box  = [lx, ly, lz, 0, 0, 0]
snapshot.particles.N         = N_total
snapshot.particles.types     = ['W', 'N', 'S']
snapshot.particles.typeid    = typeids
snapshot.particles.position  = positions
snapshot.particles.velocity  = np.zeros((N_total, 3), dtype=np.float32)
snapshot.particles.mass      = masses
snapshot.particles.slength   = slengths
snapshot.particles.density   = densities

outname = f'caprise_{NX}_{NY_TOTAL}_{NZ}_vs_{dx:.6f}_init.gsd'
with gsd.hoomd.open(name=outname, mode='w') as f:
    f.append(snapshot)

# ─── Summary ─────────────────────────────────────────────────────────────────
n_W = int(np.sum(typeids == 0))
n_N = int(np.sum(typeids == 1))
n_S = int(np.sum(typeids == 2))
A_tube   = np.pi * R_inner**2
A_domain = lx * lz
A_res    = A_domain - np.pi * R_outer**2   # reservoir cross-section (outside tube)
corr     = A_res / (A_tube + A_res)        # finite-reservoir Jurin correction

print(f'Created: {outname}')
print(f'  Total particles   : {N_total}')
print(f'  Liquid W          : {n_W}')
print(f'  Gas    N          : {n_N}')
print(f'  Solid  S          : {n_S}')
print(f'  Domain            : {lx*1e3:.1f} × {ly*1e3:.1f} × {lz*1e3:.1f} mm')
print(f'  Tube inner radius : {R_inner*1e3:.2f} mm  ({num_length} particles)')
print(f'  Tube outer radius : {R_outer*1e3:.2f} mm  ({num_length + n_wall} particles)')
print(f'  Initial liq h     : {h_res*1e3:.1f} mm  ({h_res/R_cap:.0f} R_cap)')
print(f'  A_tube / A_domain : {A_tube/A_domain:.3f}')
print(f'  Jurin correction  : ×{corr:.3f}  (rise from init level = h_J × {corr:.3f})')
sigma_ref, g_ref = 0.01, 9.81
for theta_deg in [30, 45, 60, 90, 120, 135, 150]:
    h_J = 2 * sigma_ref * np.cos(np.radians(theta_deg)) / (rho01 * g_ref * R_inner)
    print(f'  θ={theta_deg:3d}°  h_Jurin = {h_J*1e3:+6.2f} mm', end='')
    if theta_deg < 90:
        print('  (rise)')
    elif theta_deg == 90:
        print('  (neutral)')
    else:
        print('  (depression)')
