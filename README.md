# hoomd-sph3 - A SPH Implementation in HOOMD-Blue
SPH Implementation in HOOMD-Blue 3.3.0 (as of 01.08.2022)

## Libaries to install check HOOMD Webpage
additional: cereal - A C++11 library for serialization
```bash
sudo apt install libcereal-dev
```

## Main Modifications
### Technical aspects
- **pybind11** instead of old BOOST Verison
- can be build not only on NVIDIA cards
- use Type-Parameter dictionaries for class metaparameter handling
- version later than python 3.10
- better internal separation of model, logger, and integrator
- access to particle fields direct in python interface easily possible
- generally more modularized (in all aspects of the model)
- **HOOMD-Blue v3.3.0** instead of **HOOMD-Blue v1.8.2/v2.2.0**. Try to keep that up to date in the future

### Organizational apsects
- differnt Loggers, physical models, integrators can be applied added additivly (use filters instead of groups)
- differnet time integration methods (Velocity Verlet, Basic Velocity Verlet, Leap Frog)

### TODOs
- suggestion: suspensions based on already existing rigid body implementation of HOOMD-Blue
- DENSITYSUMMATION Method benchmarken (this is important) not correct in old SPH, can be seen if mean density is checked
- test limit implementation in Velocity-Verlet (easier in **HOOMD-Blue v3.5.0**)
- timer for different parts of program flow
- density computation dependent on number of ranks (is this now fixed)
- use extern gsd, DK migrated it, but there is still an
- initialize inheritence in subclasses (kernel an eos) with
```python
super().__init__(arg1, ...)
```
- add particle number density method (C++ code already implemented)
- define variables outside the particle loops (see e.g. TODO in compute_normalization_constant_solid on SinglePhaseFlow). Low hanging fruit with impact on performance. Reduce access times and - - number of requests.
- rename Tait Eq. to Cole see Paper: Cole 1948 Underwater Explosions
- get rid of not needed request on functoins. Low hanging fruit with impact on performance. Reduce access times and number of requests.
- Literature Review on kernels: which on to use when.
- plug compute_normalization_constant_solid into compute noslip
- understand why 
```python 
# in sphmodel.py 
def _attach(self):
    """
    """
    self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.half)
``` 
makes wierd errors?
- standardized input files for parameter input
- write standardized tests at least for the essentials 


### Fundamental Errors in old Code
- fictitious pressure computation, specifically the hydrostatic contribution. See Adami2012!
- Velocity Verlet not correct implemented
- Density dependent on discretisation/ this might be kernel related