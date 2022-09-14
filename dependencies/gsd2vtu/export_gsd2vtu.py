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

def import_trajectory(GSDfilename):
    """
    """
    # Input GSD file
    f = gsd.fl.GSDFile(name = GSDfilename, mode = 'rb', application = "HOOMD-SPH", schema = "hoomd", schema_version = [1,0])
    # f = gsd.pygsd.GSDFile(open('log.gsd', 'rb'))

    # Parse GSD file into a trajectory object and return
    return gsd.hoomd.HOOMDTrajectory(f)


def export_basic(GSDfilename):
    """
    filename: str
    export basic fields

    """
    print(f'Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # print(t[0].particles.position)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(snapshot.configuration.step)
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                    np.array(snapshot.particles.position.T[1]),
                    np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'     :np.array(snapshot.particles.typeid),
                    'Slength'          :np.array(snapshot.particles.slength),
                    'Mass'       :np.array(snapshot.particles.mass),
                    'Density'    :np.array(snapshot.particles.dpe.T[0]),
                    'Pressure'   :np.array(snapshot.particles.dpe.T[1]),
                    'Energy'     :np.array(snapshot.particles.dpe.T[2]),
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
    

def export_spf(GSDfilename):
    """
    filename: str
    export basic fields plus aux1 where fictious velocity is stored

    """
    print(f'Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(snapshot.configuration.step)
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                    np.array(snapshot.particles.position.T[1]),
                    np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'     :np.array(snapshot.particles.typeid),
                    'Slength'          :np.array(snapshot.particles.slength),
                    'Mass'       :np.array(snapshot.particles.mass),
                    'Density'    :np.array(snapshot.particles.dpe.T[0]),
                    'Pressure'   :np.array(snapshot.particles.dpe.T[1]),
                    'Energy'     :np.array(snapshot.particles.dpe.T[2]),
                    'Aux1x'      :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Aux1y'      :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Aux1z'      :np.array(snapshot.particles.auxiliary1.T[2]),
                    # 'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    # 'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    # 'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    # 'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    # 'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    # 'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    # 'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    # 'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    # 'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )


def export_all(GSDfilename):
    """
    filename: str
    export all fields

    """
    print(f'Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(snapshot.configuration.step)
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                    np.array(snapshot.particles.position.T[1]),
                    np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'     :np.array(snapshot.particles.typeid),
                    'Slength'          :np.array(snapshot.particles.slength),
                    'Mass'       :np.array(snapshot.particles.mass),
                    'Density'    :np.array(snapshot.particles.dpe.T[0]),
                    'Pressure'   :np.array(snapshot.particles.dpe.T[1]),
                    'Energy'     :np.array(snapshot.particles.dpe.T[2]),
                    'Aux1x'      :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Aux1y'      :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Aux1z'      :np.array(snapshot.particles.auxiliary1.T[2]),
                    'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )