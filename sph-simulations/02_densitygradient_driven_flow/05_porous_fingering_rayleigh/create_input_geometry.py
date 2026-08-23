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

Create initial GSD file(s) for the multi-mode density-driven fingering
benchmark: a DISORDERED pack of z-aligned cylindrical grains.

GEOMETRY
--------
Random sequential addition (fixed RNG seed -> reproducible) of cylinders with
radius jitter R in [0.16, 0.24] mm (mean d = 0.4 mm) until the target solid
fraction is reached. Placement rules:
  - periodic overlap checks in x (and y for --perm),
  - minimum surface-to-surface gap of 4 dx (every pore throat resolved by
    at least 4 particles).
NOTE: because the 4 dx gap scales with resolution, the achievable solid
fraction is resolution-dependent (RSA saturates near c ~ 0.25 at 10 dx per
diameter, lower at coarser dx). The achieved fraction is printed and stored
in the sidecar .npz — always use the printed porosity, not the target.

Fingering domain (default):
    Lx = 8 mm (periodic x), fluid height 12 mm + solid walls top/bottom,
    thin periodic z slice.  Target solid fraction 0.30 (phi ~= 0.70).
    A grain-free margin of one mean diameter is kept around y = 0 so the
    initial flat brine/water interface starts in open pore space.
Permeability cell (--perm):
    fully periodic 4 x 4 mm box, same disorder statistics, no walls.
--small (fingering only): Lx = 4 mm, height 6 mm — pipeline shake-down.

The grain list (centres, radii) is stored in a sidecar .npz next to the GSD
for reproducibility and post-processing.

Usage:
    python3 create_input_geometry.py <num_diameter> [seed] [--perm] [--small]

    num_diameter : fluid particles across the MEAN grain diameter (e.g. 10)
    seed         : RNG seed (default 42)
"""

import sys, math
import numpy as np
import gsd.hoomd

import hoomd
from hoomd import sph

# ─── Parameters ──────────────────────────────────────────────────────────────
num_diameter = int(sys.argv[1])
flags = sys.argv[2:]
seed  = next((int(a) for a in flags if not a.startswith('--')), 42)
perm_mode  = '--perm'  in flags
small_mode = '--small' in flags

d_mean  = 0.4e-3               # mean grain diameter                  [m]
R_lo    = 0.16e-3              # min grain radius                     [m]
R_hi    = 0.24e-3              # max grain radius                     [m]
rho0    = 1000.0               # rest density                         [kg/m³]
c_target = 0.25                # target solid (grain) fraction (upper bound;
                               # RSA + 4dx gap saturates below this at
                               # coarse resolutions — see NOTE above)   [–]

dx   = d_mean / num_diameter
mass = rho0 * dx**3

kernel  = 'WendlandC4'
slength = hoomd.sph.kernel.OptimalH[kernel] * dx
rcut    = hoomd.sph.kernel.Kappa[kernel] * slength

n_solid    = math.ceil(rcut / dx)
part_depth = math.ceil(2.5 * hoomd.sph.kernel.Kappa[kernel] * rcut / dx)

if perm_mode:
    lx, ly_fluid = 4.0e-3, 4.0e-3
elif small_mode:
    lx, ly_fluid = 4.0e-3, 6.0e-3
else:
    lx, ly_fluid = 8.0e-3, 12.0e-3

nx = round(lx / dx)
ny_fluid = round(ly_fluid / dx)
ny = ny_fluid if perm_mode else ny_fluid + 2 * n_solid
nz = part_depth

lx = nx * dx          # snap to lattice
ly_fluid = ny_fluid * dx
ly = ny * dx
lz = nz * dx

# ─── Random sequential addition of grains ────────────────────────────────────
rng = np.random.default_rng(seed)
min_gap = 4.0 * dx
area = lx * ly_fluid
centres, radii = [], []
solid_area = 0.0
attempts, max_attempts = 0, 200000

def dist_x(x1, x2):
    d = x1 - x2
    return d - lx * np.round(d / lx)

def dist_y(y1, y2):
    d = y1 - y2
    if perm_mode:
        d = d - ly_fluid * np.round(d / ly_fluid)
    return d

while solid_area < c_target * area and attempts < max_attempts:
    attempts += 1
    R = rng.uniform(R_lo, R_hi)
    xc = rng.uniform(-lx / 2, lx / 2)
    if perm_mode:
        yc = rng.uniform(-ly_fluid / 2, ly_fluid / 2)
    else:
        # keep a grain-free band of one mean diameter around the initial
        # interface (y = 0), and full grains inside the fluid region
        yc = rng.uniform(-ly_fluid / 2 + R, ly_fluid / 2 - R)
        if abs(yc) < d_mean / 2 + R:
            continue
    ok = True
    for (xo, yo), Ro in zip(centres, radii):
        if dist_x(xc, xo)**2 + dist_y(yc, yo)**2 < (R + Ro + min_gap)**2:
            ok = False
            break
    if ok:
        centres.append((xc, yc))
        radii.append(R)
        solid_area += np.pi * R**2

centres = np.array(centres)
radii   = np.array(radii)
c_actual = solid_area / area

# ─── Particle lattice & solid mask ───────────────────────────────────────────
xs = np.linspace(-lx / 2 + dx / 2, lx / 2 - dx / 2, nx)
ys = np.linspace(-ly / 2 + dx / 2, ly / 2 - dx / 2, ny)
zs = np.linspace(-lz / 2 + dx / 2, lz / 2 - dx / 2, nz)

xg, yg, zg = np.meshgrid(xs, ys, zs, indexing='ij')
positions   = np.column_stack([xg.ravel(), yg.ravel(), zg.ravel()])
n_particles = positions.shape[0]

x_pos = positions[:, 0]
y_pos = positions[:, 1]

is_solid = np.zeros(n_particles, dtype=bool)
if not perm_mode:
    is_solid |= (y_pos < -ly_fluid / 2) | (y_pos > ly_fluid / 2)

for (xc, yc), R in zip(centres, radii):
    dxr = dist_x(x_pos, xc)
    dyr = dist_y(y_pos, yc)
    is_solid |= (dxr**2 + dyr**2) < R**2

typeid = np.where(is_solid, 1, 0).astype(np.int32)

n_fluid   = int(np.sum(typeid == 0))
n_solid_p = int(np.sum(typeid == 1))

in_lattice = np.abs(y_pos) < ly_fluid / 2 if not perm_mode else np.ones(n_particles, bool)
phi = float(np.sum((typeid == 0) & in_lattice)) / float(np.sum(in_lattice))

# ─── Write GSD + sidecar grain list ──────────────────────────────────────────
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

tag = 'perm' if perm_mode else ('small' if small_mode else 'pack')
init_filename = f'porous_ra_{tag}_{num_diameter}_seed{seed}_init.gsd'
with gsd.hoomd.open(name=init_filename, mode='w') as f:
    f.append(snapshot)

np.savez(init_filename.replace('_init.gsd', '_grains.npz'),
         centres=centres, radii=radii, lx=lx, ly_fluid=ly_fluid,
         seed=seed, dx=dx, phi=phi)

print(f'Written {init_filename}: {n_fluid} fluid, {n_solid_p} solid ({n_particles} total)')
print(f'  Mode: {"perm cell (fully periodic)" if perm_mode else ("small shake-down" if small_mode else "full pack")},  seed = {seed}')
print(f'  Domain:  {lx*1e3:.2f} × {ly*1e3:.2f} × {lz*1e3:.3f} mm  (nx={nx}, ny={ny}, nz={nz})')
print(f'  dx = {dx*1e6:.1f} µm,  {len(radii)} grains,  R ∈ [{R_lo*1e3:.2f}, {R_hi*1e3:.2f}] mm')
print(f'  Solid fraction: target {c_target:.2f}, achieved {c_actual:.3f}  '
      f'(porosity φ = {phi:.3f})')
print(f'  Min throat gap: {min_gap/dx:.0f} dx')
if not perm_mode:
    print(f'  Fluid height H = {ly_fluid*1e3:.1f} mm, grain-free interface band ±{d_mean/2*1e3:.1f} mm')
    print(f'  Concentration assigned in run script (flat interface at y = 0).')
