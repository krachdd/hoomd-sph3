#!/home/david/anaconda3/envs/sph3/bin/python
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
# import gsd.hoomd

cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu)

integrator = hoomd.sph.Integrator(dt=0.005)

# sim.operations.integrator = hoomd.sph.Integrator()
