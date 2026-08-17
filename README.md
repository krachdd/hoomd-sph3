# hoomd-sph3 — SPH Implementation in HOOMD-Blue

[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](./LICENSE)
[![Release V1.0](https://img.shields.io/badge/release-V1.0-green)](https://github.com/krachdd/hoomd-sph3/releases/tag/v1.0)
[![HOOMD-Blue](https://img.shields.io/badge/HOOMD--Blue-5.2.0-orange)](https://hoomd-blue.readthedocs.io/en/latest/)

**Language:** C++, Python

---

## Overview

**hoomd-sph3** is an open-source Smoothed Particle Hydrodynamics (SPH) solver built on top of [HOOMD-Blue 5.2.0](https://hoomd-blue.readthedocs.io/en/latest/). It provides modular, extensible C++ components for particle-based fluid simulation, together with a comprehensive suite of benchmark cases covering single-phase, two-phase, free-surface, non-Newtonian, and density-gradient-driven flows.

The solver runs on CPUs (MPI domain decomposition supported; GPU kernels are not yet implemented) via HOOMD-Blue's plugin architecture, and uses two custom I/O libraries — **gsd-sph** (serial) and **pgsd-sph** (MPI-parallel) — for efficient snapshot reading and writing.

---

## Features

- **Multiple time integrators** — Velocity Verlet, Leap Frog, KickDriftKick (KDK), KickDriftKick with Transport Velocity (KDK-TV)
- **Modular physical model system** — filter-based particle selection; clean separation of EOS, viscosity model, and integrator
- **Density computation** — summation-based and continuity-equation-based density update routines
- **Extensible kernel library** — WendlandC2, WendlandC4, and additional SPH kernels
- **Equation-of-state (EOS)** — weakly-compressible (Tait) EOS and user-extensible interface
- **Free-surface flows** — Shepard completeness ratio, contact-angle enforcement (Huber et al. 2016), mean-curvature estimation, and CSF surface tension (Colagrossi & Landrini 2003)
- **Non-Newtonian rheology** — power-law and Bingham viscosity models
- **Two-phase flows** — surface-tension and Transport-Velocity formulations
- **Density-gradient-driven flows** — Boussinesq buoyancy approximation, thermal Couette, heated cavity, and Rayleigh-Taylor instability
- **Parallel I/O** — MPI-IO via pgsd-sph eliminates the serial gather-on-root bottleneck
- **Density diffusion** — Molteni & Colagrossi (2009) diffusion and Antuono/Marrone δ-SPH with L-matrix renormalized density gradients (`model.activateDeltaSPH`)
- **Adaptive time stepping** — `hoomd.sph.update.AdaptiveTimestep` recomputes Δt from the CFL, viscous, and body-force conditions using the MPI-reduced maximum fluid speed
- **Helper modules** — GSD-to-VTU conversion, input-geometry readers, and diagnostic utilities
- **Reproducible benchmarks** — standardized create/run scripts for all benchmark cases
- **Unit tests** — `hoomd/sph/pytest/` covers kernel normalization/derivatives, EOS consistency, lattice completeness, and small end-to-end solver runs (`python -m pytest hoomd/sph/pytest` in the build tree)

---

## Repository Structure

```
hoomd-sph3/
├── hoomd-blue/              # HOOMD-Blue 5.2.0 (submodule)
├── dependencies/
│   ├── gsd-sph/             # GSD 3.4.2 extended with SPH particle fields
│   └── pgsd-sph/            # PGSD 3.2.0 — MPI-parallel GSD fork
├── helper_modules/
│   ├── gsd2vtu/             # GSD → VTU/PVD conversion for ParaView
│   ├── sph_helper.py        # SPH utility functions
│   ├── sph_info.py          # Simulation info / diagnostics
│   ├── read_input_fromtxt.py
│   └── delete_solids_initial_timestep.py
├── sph-simulations/
│   ├── 00_singlephaseflow_benchmarks/
│   ├── 01_twophaseflow_benchmarks/
│   ├── 02_densitygradient_driven_flow/
│   ├── 03_nonnewtonian_benchmarks/
│   ├── 04_free_surface_flow/
│   └── 05_io_benchmarks/
├── buildall.sh              # Full build script (pgsd → gsd → HOOMD-Blue)
├── fast_buildall.sh         # Incremental build (skips clean)
├── addsphpath.sh            # Sets PYTHONPATH for all modules
├── requirements.txt         # Python dependencies (pip)
└── LICENSE
```

---

## Dependencies

### System / compiled

| Dependency | Version | Purpose |
|---|---|---|
| [HOOMD-Blue](https://hoomd-blue.readthedocs.io/) | 5.2.0 | Simulation engine |
| [cereal](https://uscilab.github.io/cereal/) | any | C++11 serialization |
| MPI (OpenMPI ≥ 4 or MPICH ≥ 3) | — | Parallel execution and I/O |
| CMake | ≥ 3.14 | Build system |

**Linux (Debian/Ubuntu):**
```bash
sudo apt install libcereal-dev python3-dev libbz2-dev libopenmpi-dev cmake
```

### Python

| Package | Version | Purpose |
|---|---|---|
| Python | ≥ 3.8 | Scripting and post-processing |
| numpy | 2.2.1 | Array operations |
| scipy | 1.15.1 | Scientific utilities |
| mpi4py | 4.0.1 | Python MPI bindings |
| vtk / pyevtk | 9.3.1 / 1.6.0 | Visualization output |
| pybind11 | 2.13.6 | C++ Python bindings |
| Cython | 3.0.11 | PGSD Python wrapper compilation |
| matplotlib | 3.10.0 | Plotting |

Install all Python requirements:
```bash
pip install -r requirements.txt
```

---

## Getting Started

### 1. Clone with submodules

```bash
git clone --recurse-submodules https://github.com/krachdd/hoomd-sph3.git
cd hoomd-sph3
```

### 2. Build all components

The `buildall.sh` script builds pgsd-sph, gsd-sph, and HOOMD-Blue in the correct order:

```bash
bash buildall.sh
```

For incremental rebuilds (without `rm -rf build`):
```bash
bash fast_buildall.sh
```

### 3. Set up the Python path

Source `addsphpath.sh` to add all build directories and helper modules to `PYTHONPATH`:

```bash
source addsphpath.sh
```

Add this to your shell profile or job script to make it permanent.

### 4. Run a benchmark

Each benchmark follows the same two-step pattern:

```bash
# Step 1 — create the initial GSD geometry
python3 sph-simulations/00_singlephaseflow_benchmarks/01_poiseuille/create_input_geometries_tube.py <num_length>

# Step 2 — run the simulation
python3 sph-simulations/00_singlephaseflow_benchmarks/01_poiseuille/run_tube.py <num_length> <init_gsd_file> [steps]
```

See the `README.md` inside each benchmark folder for case-specific arguments.

---

## Benchmark Suites

### 00 — Single-Phase Flow

Classical incompressible/weakly-compressible SPH benchmarks:

| Case | Description |
|---|---|
| 01_poiseuille | Poiseuille flow in a tube |
| 02_parallel_plates | Couette/Poiseuille between parallel plates |
| 03_couette_flow | Lid-driven Couette flow |
| 04_channel_flow | Pressure-driven channel flow |
| 05_liddriven_cavity | Lid-driven cavity (Re = 1 – 10 000) |
| 06_spherepackings | Flow through random sphere packings |
| 07_adami_cylinder1 | Flow around a cylinder — Adami BC type 1 |
| 08_adami_cylinder2 | Flow around a cylinder — Adami BC type 2 |
| 09_drag_cylinder | Cylinder drag coefficient |
| 10_hydrostatic | Hydrostatic pressure test |
| 11_taylor_green | Taylor–Green vortex decay |
| 12_stokes_second | Stokes second problem (oscillating plate) |
| 13_bcc_spherepacking_permeability | BCC sphere-packing permeability vs. Re |

### 01 — Two-Phase Flow

Surface-tension and interfacial flow benchmarks:

| Case | Description |
|---|---|
| 01_layered_couette | Layered two-fluid Couette flow |
| 02_static_droplet | Static droplet — Laplace pressure jump |
| 03_capillary_rise | Capillary rise in a tube |
| 04_rising_bubble | Rising bubble in a heavier fluid |
| 05_rayleigh_taylor | Rayleigh–Taylor instability |
| 06_sessile_droplet_shear | Sessile droplet under shear |
| 07_sessile_droplet_snapback | Sessile droplet snap-back |
| 08_h2_bubble_brine_shear | H₂ bubble in brine under shear |

### 02 — Density-Gradient-Driven Flow

Buoyancy and Boussinesq approximation benchmarks:

| Case | Description |
|---|---|
| 01_thermal_couette | Thermal Couette flow |
| 02_heated_cavity | Differentially heated cavity |
| 03_boussinesq_rayleigh_taylor | Boussinesq Rayleigh–Taylor instability |

### 03 — Non-Newtonian Benchmarks

| Case | Description |
|---|---|
| 01_powerlaw_poiseuille | Power-law fluid in Poiseuille flow |
| 02_bingham_poiseuille | Bingham fluid in Poiseuille flow |
| 03_twophase_powerlaw_couette | Two-phase power-law Couette flow |
| 04_bingham_slump_test | Bingham slump test |

### 04 — Free-Surface Flow

Benchmarks for the `SinglePhaseFlowFS` solver (free-surface + contact-angle + CSF):

| Case | Description |
|---|---|
| 01_hydrostatic_fs | Hydrostatic free surface |
| 02_dam_break | Dam-break (Martin & Moyce 1952) |
| 03_sessile_droplet | Sessile droplet on a substrate |

### 05 — I/O Benchmarks

Performance and correctness tests for serial (GSD) and parallel (PGSD) I/O:

| Case | Description |
|---|---|
| 00_create_domains | Domain generation utilities |
| 01_gsd | Serial GSD read/write benchmarks |
| 02_pgsd | MPI-parallel PGSD read/write benchmarks |

---

## I/O Libraries

### gsd-sph (v3.4.2)

An extension of the upstream [GSD library](https://github.com/glotzerlab/gsd) that adds the SPH-specific particle fields (`slength`, `density`, `pressure`, `energy`, `auxiliary1–4`) to `gsd.hoomd.ParticleData`. Serial POSIX I/O; intended for single-process runs and post-processing.

See [`dependencies/gsd-sph/README.md`](dependencies/gsd-sph/README.md).

### pgsd-sph (v3.2.0)

A fork of GSD in which every POSIX `read`/`write` call is replaced with `MPI_File_*` collective I/O. All MPI ranks write their own particle partitions simultaneously; no data is gathered to rank 0. This eliminates the serial I/O bottleneck for large-scale parallel runs.

See [`dependencies/pgsd-sph/README.md`](dependencies/pgsd-sph/README.md).

---

## Helper Modules

| Module | Description |
|---|---|
| `gsd2vtu/` | Converts GSD trajectory files to VTU/PVD format for ParaView |
| `sph_helper.py` | SPH kernel functions, particle spacing utilities |
| `sph_info.py` | Print simulation metadata and particle field statistics |
| `read_input_fromtxt.py` | Read raw particle data from text-format input files |
| `delete_solids_initial_timestep.py` | Remove solid-wall particles from the first GSD frame |

---

## Contributing

- Follow the filter-based (not group-based) particle selection pattern throughout.
- Document new integrators, physical models, and benchmark cases with a `README.md`.
- Keep `create_*` and `run_*` scripts self-contained and parameterized by resolution.
- Maintain test cases and update this README when adding new benchmark categories.

---

## License

BSD 3-Clause "New" or "Revised" License — see [LICENSE](./LICENSE) for full details.

---

## Contact

**Developer:**
[David Krach](https://www.mib.uni-stuttgart.de/institute/team/Krach/)
Institute of Applied Mechanics (MIB), University of Stuttgart
E-mail: [david.krach@mib.uni-stuttgart.de](mailto:david.krach@mib.uni-stuttgart.de)
