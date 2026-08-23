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

Create initial GSD file(s) for the pore-scale density-driven fingering benchmark
(brine over water in a square array of cylindrical grains).

GEOMETRY
--------
Quasi-2D vertical cell filled with a square lattice of z-aligned solid
cylinders (grains):

    grain spacing    a = 0.8 mm
    grain radius     R = 0.2 mm   (diameter d = 0.4 mm)
    solid fraction   c = pi R^2 / a^2 ~= 0.196,  porosity phi ~= 0.80

Fingering domain (default):
    width   Lx = 5 a = 4 mm   -> one perturbation wavelength, periodic in x
    fluid height  10 a = 8 mm (y in [-4, 4] mm) + solid walls top/bottom
    thin periodic z slice (pseudo-3D)
    Grain centres at x_c = -Lx/2 + (i+1/2) a, y_c = (j+1/2) a (j in Z), so the
    initial brine/water interface at y = 0 runs through a pore-body row.

Permeability unit cell (--perm):
    fully periodic 2a x 2a x thin-z box containing 2x2 grains, same a, R, dx.
    Used by run_permeability.py to measure the Darcy permeability K of the
    identical lattice (Drummond & Tahir 1984 analytic estimate as reference).

The scalar field (brine concentration C in aux4.x) is assigned in the run
script via cpu_local_snapshot.

Usage:
    python3 create_input_geometry.py <num_diameter> [--perm] [--cells-x N]

    num_diameter : fluid particles across one grain diameter d (e.g. 10)
    --perm       : write the fully periodic permeability unit cell instead
    --cells-x N  : domain width in grain columns (default 5) -> perturbation
                   wavelength Lx = N*a; use 3 and 7 for the sigma(k) study
                   (non-default widths are tagged _wN in the filename)
"""

import sys, math
import numpy as np
import gsd.hoomd

import hoomd
from hoomd import sph

# ─── Parameters ──────────────────────────────────────────────────────────────
num_diameter = int(sys.argv[1])
flags        = sys.argv[2:]
perm_mode    = '--perm' in flags

# Domain width in grain columns (fingering mode): Lx = cells_x * a sets the
# perturbation wavelength — vary it for the sigma(k) dispersion study
# (e.g. --cells-x 3 and --cells-x 7 around the default 5).
cells_x = 5
if '--cells-x' in flags:
    cells_x = int(flags[flags.index('--cells-x') + 1])

a_grain = 0.8e-3           # grain lattice spacing                    [m]
R_grain = 0.2e-3           # grain radius                             [m]
d_grain = 2.0 * R_grain    # grain diameter                           [m]
rho0    = 1000.0           # rest density                             [kg/m³]

dx   = d_grain / num_diameter
mass = rho0 * dx**3

kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut    = hoomd.sph.kernel.Kappa[kernel] * slength

n_solid    = math.ceil(rcut / dx)   # wall thickness [particle layers]
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut / dx)

n_cells_x = 2 if perm_mode else cells_x   # grain columns
n_cells_y = 2 if perm_mode else 10        # grain rows in the fluid region

lx = n_cells_x * a_grain
ly_fluid = n_cells_y * a_grain

nx = round(lx / dx)
ny_fluid = round(ly_fluid / dx)
ny = ny_fluid if perm_mode else ny_fluid + 2 * n_solid
nz = part_depth

ly = ny * dx
lz = nz * dx

# ─── Particle positions (half-integer offsets — no particles on box faces) ───
xs = np.linspace(-lx / 2 + dx / 2, lx / 2 - dx / 2, nx)
ys = np.linspace(-ly / 2 + dx / 2, ly / 2 - dx / 2, ny)
zs = np.linspace(-lz / 2 + dx / 2, lz / 2 - dx / 2, nz)

xg, yg, zg = np.meshgrid(xs, ys, zs, indexing='ij')
positions   = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
n_particles = positions.shape[0]

x_pos = positions[:, 0]
y_pos = positions[:, 1]

# ─── Solid mask: walls (fingering mode) + cylindrical grains ─────────────────
is_solid = np.zeros(n_particles, dtype=bool)
if not perm_mode:
    is_solid |= (y_pos < -ly_fluid / 2) | (y_pos > ly_fluid / 2)

# Grain centres: x_c = -lx/2 + (i+1/2) a  (periodic in x),
#                y_c = (j+1/2) a          (symmetric about y=0, no grain on
#                                          the interface row y=0)
x_centres = -lx / 2 + (np.arange(n_cells_x) + 0.5) * a_grain
j_lo = math.floor((-ly / 2) / a_grain - 0.5) - 1
j_hi = math.ceil((ly / 2) / a_grain - 0.5) + 1
y_centres = (np.arange(j_lo, j_hi + 1) + 0.5) * a_grain
# keep grains whose centres lie inside the fluid region (with margin R)
y_centres = y_centres[np.abs(y_centres) < ly_fluid / 2 + R_grain]

for yc in y_centres:
    dy2 = (y_pos - yc)**2
    for xc in x_centres:
        dxr = x_pos - xc
        # periodic images in x
        dxr = dxr - lx * np.round(dxr / lx)
        is_solid |= (dxr**2 + dy2) < R_grain**2

typeid = np.where(is_solid, 1, 0).astype(np.int32)

n_fluid   = int(np.sum(typeid == 0))
n_solid_p = int(np.sum(typeid == 1))

# ─── Porosity & Drummond–Tahir reference permeability ────────────────────────
# Porosity of the grain lattice (excluding wall layers)
in_lattice = np.abs(y_pos) < ly_fluid / 2 if not perm_mode else np.ones(n_particles, bool)
phi = float(np.sum((typeid == 0) & in_lattice)) / float(np.sum(in_lattice))

c_solid = np.pi * R_grain**2 / a_grain**2
# Drummond & Tahir (1984), square array of cylinders, transverse flow:
#   K = R^2/(8c) * ( -ln c - 1.476 + 2c - 1.774 c^2 + 4.076 c^3 )
K_DT = R_grain**2 / (8.0 * c_solid) * (
    -np.log(c_solid) - 1.476 + 2.0 * c_solid
    - 1.774 * c_solid**2 + 4.076 * c_solid**3)

# ─── Write GSD ───────────────────────────────────────────────────────────────
snapshot = gsd.hoomd.Frame()
snapshot.configuration.box  = [lx, ly, lz, 0, 0, 0]
snapshot.particles.N        = n_particles
snapshot.particles.types    = ['F', 'S']
snapshot.particles.typeid   = typeid
snapshot.particles.position = positions.astype(np.float32)
snapshot.particles.velocity = np.zeros((n_particles, 3), dtype=np.float32)
snapshot.particles.mass     = np.full(n_particles, mass,    dtype=np.float32)
snapshot.particles.slength  = np.full(n_particles, slength, dtype=np.float32)
snapshot.particles.density  = np.full(n_particles, rho0,    dtype=np.float32)

# Non-default widths get a _w{N} tag so sigma(k) variants coexist with the
# default geometry (job scripts glob on 'porous_finger_<nd>_dx_*').
tag = 'perm' if perm_mode else ('finger' if cells_x == 5 else f'finger_w{cells_x}')
init_filename = f'porous_{tag}_{num_diameter}_dx_{dx:.2e}_init.gsd'
with gsd.hoomd.open(name=init_filename, mode='w') as f:
    f.append(snapshot)

print(f'Written {init_filename}: {n_fluid} fluid, {n_solid_p} solid ({n_particles} total)')
print(f'  Mode: {"permeability unit cell (fully periodic)" if perm_mode else "fingering domain (walls in y)"}')
print(f'  Domain:  {lx*1e3:.2f} × {ly*1e3:.2f} × {lz*1e3:.3f} mm  (nx={nx}, ny={ny}, nz={nz})')
print(f'  dx = {dx*1e6:.1f} µm,  slength = {slength*1e6:.1f} µm,  rcut = {rcut*1e6:.1f} µm')
print(f'  Grains: a = {a_grain*1e3:.2f} mm, d = {d_grain*1e3:.2f} mm '
      f'({num_diameter} dx per diameter), solid fraction c = {c_solid:.3f}')
print(f'  Porosity (particle count) φ = {phi:.3f}   (analytic 1-c = {1.0-c_solid:.3f})')
print(f'  Drummond–Tahir K = {K_DT:.3e} m²')
if not perm_mode:
    print(f'  Fluid height = {ly_fluid*1e3:.1f} mm,  wall layers = {n_solid}')
    print(f'  Concentration assigned in run script via cpu_local_snapshot.')
