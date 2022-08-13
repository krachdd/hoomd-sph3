#!/home/david/anaconda3/envs/sph3/bin/python
import hoomd
from hoomd import *
from hoomd import sph
import numpy as np
import itertools
import gsd.hoomd


device = hoomd.device.CPU(notice_level=10)
# device = hoomd.device.CPU(notice_level=2)
sim = hoomd.Simulation(device=device)
# if device.communicator.rank == 0:
#     print("hoomd.version.mpi_enabled: {0}".format(hoomd.version.mpi_enabled))


filename = "test_parallelplates.gsd"
dumpname = "run_{0}".format(filename)

sim.create_state_from_gsd(filename = filename, domain_decomposition=(2, 2, 1))

print(sim.state.box)

# Print the domain decomposition.
domain_decomposition = sim.state.domain_decomposition
if device.communicator.rank == 0:
    print(domain_decomposition)

# Print the location of the split planes.
split_fractions = sim.state.domain_decomposition_split_fractions
if device.communicator.rank == 0:
    print(split_fractions)

# Print the number of particles on each rank.
with sim.state.cpu_local_snapshot as snap:
    N = len(snap.particles.position)
    print(f'{N} particles on rank {device.communicator.rank}')


# hoomd.write.GSD.write(filename = dumpname, state = sim.state, mode = 'wb')