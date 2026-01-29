#!/usr/bin/env python3

"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

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
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

"""


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
f = gsd.fl.GSDFile(name = sys.argv[1], mode = 'r', application = "HOOMD-SPH", schema = "hoomd", schema_version = [1,0])
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
   filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
   
   vtk(filename, np.array(snapshot.particles.position.T[0]),
                 np.array(snapshot.particles.position.T[1]),
                 np.array(snapshot.particles.position.T[2]),
       data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
           'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
           'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
           'TypeId'     :np.array(snapshot.particles.typeid),
           'Slength'    :np.array(snapshot.particles.slength),
           'Mass'       :np.array(snapshot.particles.mass),
           'Density'    :np.array(snapshot.particles.density),
           'Pressure'   :np.array(snapshot.particles.pressure),
           'Energy'     :np.array(snapshot.particles.energy),
               #'Aux1x'      :np.array(snapshot.particles.auxiliary1.T[0]),
               #'Aux1y'      :np.array(snapshot.particles.auxiliary1.T[1]),
               #'Aux1z'      :np.array(snapshot.particles.auxiliary1.T[2]),
               #'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
               #'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
               #'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
               #'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
               #'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
               #'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
               #'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
               #'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
               #'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
           },
       )
    

