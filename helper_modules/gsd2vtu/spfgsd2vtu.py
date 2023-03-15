#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""


# --- HEADER ---------------------------------------------------
import gsd.fl
import gsd.hoomd
import gsd.pygsd
import numpy as np
import sys
import os
from pyevtk.hl import pointsToVTK as vtk
#--------------------------------------------------------------


# Input GSD file
f = gsd.fl.GSDFile(name = sys.argv[1], mode = 'rb', application = "HOOMD-SPH", schema = "hoomd", schema_version = [1,0])
# f = gsd.pygsd.GSDFile(open('log.gsd', 'rb'))

# Parse GSD file into a trajectory object
t = gsd.hoomd.HOOMDTrajectory(f)

# print(t[0].particles.position)

# Run loop over all snapshots
count = 0
for snapshot in t:
   count += 1
   
   pname = sys.argv[1].replace('.gsd','')
   
   if not os.path.exists(pname):
      os.makedirs(pname)
        
   # Define VTU export filename
   filename = pname+'/'+pname+'_'+str(snapshot.configuration.step)
   
   vtk(filename, np.array(snapshot.particles.position.T[0]),
                 np.array(snapshot.particles.position.T[1]),
                 np.array(snapshot.particles.position.T[2]),
       data = {'Velocity x'      :np.array(snapshot.particles.velocity.T[0]),
           'Velocity y'          :np.array(snapshot.particles.velocity.T[1]),
           'Velocity z'          :np.array(snapshot.particles.velocity.T[2]),
           'TypeId'              :np.array(snapshot.particles.typeid),
           'Slength'             :np.array(snapshot.particles.slength),
           'Mass'                :np.array(snapshot.particles.mass),
           'Density'             :np.array(snapshot.particles.density),
           'Pressure'            :np.array(snapshot.particles.pressure),
           'Energy'              :np.array(snapshot.particles.energy),
           'Ficticious Velx'     :np.array(snapshot.particles.auxiliary1.T[0]),
           'Ficticious Vely'     :np.array(snapshot.particles.auxiliary1.T[1]),
           'Ficticious Velz'     :np.array(snapshot.particles.auxiliary1.T[2]),
           },
       )
    

