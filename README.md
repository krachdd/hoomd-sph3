# hoomd-sph3 - A SPH Implementation in HOOMD-Blue

**License:** BSD 3-Clause  
**Language:** C++, python

## Overview

hoomd-sph3 is an open-source implementation of Smoothed Particle Hydrodynamics (SPH) in [HOOMD-Blue](https://hoomd-blue.readthedocs.io/en/latest/), focused on flexible particle-based fluid simulations in C++. The code targets HOOMD-Blue v5.2.0 and provides modular components for various physical models, integrators, and logging mechanisms.

## Features

- Multiple time integration methods (Velocity Verlet, Leap Frog, etc.)
- Modular logger and physical model system using filters
- Density computation and summation routines
- Parallel I/O routines for snapshot reading and writing
- Extensible kernel and equation-of-state (EOS) 
- Standardized input and test files for reproducibility

## Dependencies

- [HOOMD-Blue 5.2.0](https://hoomd-blue.readthedocs.io/en/latest/)
- [cereal](https://uscilab.github.io/cereal/) (C++11 serialization)
- GSD (General Simulation Data) version 3.4.2
- PGSD version 3.2.0
- Python 3 (for auxiliary scripting)
- Standard C++ build tools

**Linux installation example:**
```bash
sudo apt install libcereal-dev python3-dev libbz2-dev
```
Refer to the HOOMD-Blue webpage for more details.

## Getting Started

1. **Clone the repository:**
   ```bash
   git clone --recurse-submodules https://github.com/krachdd/hoomd-sph3.git
   ```
2. **Install dependencies:**
   Ensure HOOMD-Blue, cereal, GSD, PGSD, and other requirements are installed as described above.

3. **Build the project:**
   Follow HOOMD-Blue's build instructions, integrating the `hoomd-sph3` modules as needed.

4. **Run Simulation:**
   Prepare standard input files for your scenario and run using the provided scripts or binary.
   > Example usage scripts and standardized input templates are available in the repository.

## Contributing

- Document new modules, integrators, and changes thoroughly.
- Use filters rather than groups for selecting particle sets.
- Maintain test suites and README documentation with all essential updates.

## License

BSD 3-Clause "New" or "Revised" License. See [LICENSE](./LICENSE) for full details.

## Contact

**Developer:**  
[David Krach](https://www.mib.uni-stuttgart.de/institute/team/Krach/)  
E-mail: [david.krach@mib.uni-stuttgart.de](mailto:david.krach@mib.uni-stuttgart.de)







