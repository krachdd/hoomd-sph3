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

Sessile Droplet Shear + Snap-back — initial geometry creation.

BENCHMARK DESCRIPTION
---------------------
A cube of fluid 'W' (droplet phase) is placed in contact with the bottom
solid wall, immersed in the ambient fluid 'N'.  Surface tension (σ) and
the prescribed contact angle (θ = 90°) morph the cube into a sessile
spherical-cap droplet during the relaxation phase of the simulation.

After equilibration the upper solid wall is moved at velocity U_wall
(applied at runtime in run_sessile_droplet_snapback.py) to drive Couette-like
shear flow that deforms the sessile droplet.  Contact-angle hysteresis
(θ_rec=75°, θ_adv=105°) pins the contact line so that the droplet recovers
elastically when the shear is removed in the third (snap-back) phase.

Domain:
  x : 4·lref    periodic  (room for elongation under shear)
  y : 2·lref     solid walls  (bottom stationary, top moving in shear phase)
  z : 2·lref     periodic

Droplet cube:
  side length  a = (2π/3)^(1/3)·lref  ≈ 1.28·lref  (volume-conserving for hemisphere)
  centred at   x = 0, z = 0
  bottom face  at y = −lref  (in contact with the bottom solid wall)

Solid walls:
  n_solid = 3 layers at bottom and top  (Adami 2012 no-slip BC)
  All solid velocities are zero in the init file.
  The top wall velocity is set at runtime before the shear phase.

Expected sessile droplet shape (θ = 90°, hemisphere):
  R_sessile = H_flu/2 = lref  (by construction)

Key dimensionless numbers (default parameters):
  Bond number  Bo = ρ g R² / σ   ≈ 1.0    (gravity and surface tension comparable)
  Capillary    Ca = μ U_wall / σ  = 0.001  (surface-tension dominated)
  Reynolds     Re = ρ U_wall lref / μ       = 10

Usage:
    python3 create_input_geometry.py <num_length>

    num_length : integer resolution — particles across lref (e.g. 20, 40)

Output:
    sessile_snapback_<nx>_<ny>_<nz>_vs_<dx>_init.gsd
"""

import sys
import numpy as np
import gsd.hoomd
import hoomd
from hoomd import sph

# ─── Parameters ──────────────────────────────────────────────────────────────
num_length = int(sys.argv[1]) if len(sys.argv) > 1 else 20

lref    = 0.001                    # reference length                      [m]
rho01   = 1000.0                   # rest density phase W (droplet)        [kg/m³]
rho02   = 1000.0                   # rest density phase N (ambient)        [kg/m³]
dx      = lref / num_length        # particle spacing                      [m]
mass1   = rho01 * dx**3            # particle mass phase W                 [kg]
mass2   = rho02 * dx**3            # particle mass phase N                 [kg]
n_solid = 3                        # solid wall layers per side

kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel] * dx   # smoothing length     [m]

# ─── Box geometry ────────────────────────────────────────────────────────────
# Fluid region : 4·lref (x) × 2·lref (y) × 2·lref (z)
# Solid layers : n_solid rows at bottom and top in y.
# x, z         : periodic.   y : bounded by solid walls.
nx     = 4 * num_length            # wider in x for shear-driven elongation
ny_flu = 2 * num_length            # fluid rows in y  (channel height = 2·lref)
ny     = ny_flu + 2 * n_solid      # total rows including solid walls
nz     = 2 * num_length

lx = nx * dx                       # 4·lref
ly = ny * dx
lz = nz * dx                       # 2·lref

H_flu = ny_flu * dx                # fluid channel height = 2·lref         [m]

# Droplet cube: volume-conserving for a hemisphere of radius R = H_flu/2
R_drop_target = H_flu / 2                              # = lref             [m]
a_cube = (2 * np.pi * R_drop_target**3 / 3) ** (1/3)  # ≈ 1.28·lref        [m]
y_bot = -(H_flu / 2 + n_solid * dx)
y_top =  (H_flu / 2 + n_solid * dx)

# ─── Particle lattice ────────────────────────────────────────────────────────
x_arr = np.linspace(-lx / 2 + dx / 2, lx / 2 - dx / 2, nx, endpoint=True)
y_arr = np.linspace(y_bot + dx / 2,   y_top - dx / 2,   ny, endpoint=True)
z_arr = np.linspace(-lz / 2 + dx / 2, lz / 2 - dx / 2, nz, endpoint=True)

xg, yg, zg = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
positions = np.stack((xg.ravel(), yg.ravel(), zg.ravel()), axis=1).astype(np.float32)

n_particles = len(positions)

# Default: all 'N' (ambient fluid, index 1)
typeids    = np.ones(n_particles, dtype=np.int32)
velocities = np.zeros((n_particles, 3), dtype=np.float32)
masses     = np.full(n_particles, mass2, dtype=np.float32)
slengths   = np.full(n_particles, slength, dtype=np.float32)
densities  = np.full(n_particles, rho02, dtype=np.float32)

xp = positions[:, 0]
yp = positions[:, 1]
zp = positions[:, 2]

# ─── Droplet cube ('W', index 0) ─────────────────────────────────────────────
# Centred at x = 0, z = 0; bottom face at y = −H_flu/2.
half_a      = a_cube / 2
inside_cube = (
    (np.abs(xp) <= half_a) &
    (yp >= -H_flu / 2) &
    (yp <  -H_flu / 2 + a_cube) &
    (np.abs(zp) <= half_a)
)
typeids[inside_cube]   = 0    # 'W'
masses[inside_cube]    = mass1
densities[inside_cube] = rho01

# ─── Solid walls ('S', index 2) ──────────────────────────────────────────────
# Both walls initialised with zero velocity.
# The top wall velocity is set at runtime before the shear phase.
solid_mask = (yp < -H_flu / 2) | (yp >= H_flu / 2)
typeids[solid_mask]   = 2
masses[solid_mask]    = mass1
densities[solid_mask] = rho01
# velocities already zero

# ─── GSD snapshot ────────────────────────────────────────────────────────────
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box  = [lx, ly, lz, 0, 0, 0]
snapshot.particles.N         = n_particles
snapshot.particles.types     = ['W', 'N', 'S']
snapshot.particles.typeid    = typeids
snapshot.particles.position  = positions
snapshot.particles.velocity  = velocities
snapshot.particles.mass      = masses
snapshot.particles.slength   = slengths
snapshot.particles.density   = densities

init_filename = f'sessile_snapback_{nx}_{ny}_{nz}_vs_{dx:.6f}_init.gsd'
with gsd.hoomd.open(name=init_filename, mode='w') as f:
    f.append(snapshot)

# ─── Summary ─────────────────────────────────────────────────────────────────
n_drop   = int(np.sum(typeids == 0))
n_amb    = int(np.sum(typeids == 1))
n_sol    = int(np.sum(typeids == 2))
V_cube   = n_drop * dx**3
R_hemi   = (3 * V_cube / (2 * np.pi)) ** (1 / 3)
# Reference dimensionless numbers
sigma_r  = 0.01    # surface tension used in run script  [N/m]
gy_r     = 9.81    # gravity magnitude                   [m/s²]
U_wall_r = 0.01    # shear velocity used in run script   [m/s]
mu_r     = 0.001   # viscosity used in run script        [Pa·s]
Bo = rho01 * gy_r * R_hemi**2 / sigma_r
Ca = mu_r * U_wall_r / sigma_r
Re = rho01 * U_wall_r * lref / mu_r

print(f'Created {init_filename}')
print(f'  Total particles   : {n_particles}')
print(f'  Droplet W (cube)  : {n_drop}  (side a = {a_cube*1e3:.3f} mm, '
      f'≈{int(round(a_cube/dx))} layers per edge)')
print(f'  Ambient N         : {n_amb}')
print(f'  Solid S           : {n_sol}  ({n_solid} layers each wall)')
print(f'  Domain (fluid)    : {lx*1e3:.2f} × {H_flu*1e3:.2f} × {lz*1e3:.2f} mm  '
      f'(4·lref × 2·lref × 2·lref)')
print(f'  dx                : {dx*1e6:.1f} µm   (lref/dx = {num_length})')
print(f'  R_target          : {R_drop_target*1e3:.3f} mm  (= H_flu/2 = lref)')
print(f'  R_sessile ≈ {R_hemi*1e3:.3f} mm  (hemisphere, θ = 90°)')
print(f'  Bond number  Bo   = {Bo:.4f}  (σ={sigma_r} N/m, g={gy_r} m/s²)')
print(f'  Capillary    Ca   = {Ca:.4f}  (U_wall={U_wall_r} m/s)')
print(f'  Reynolds     Re   = {Re:.1f}')
