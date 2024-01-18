#!/usr/bin/env python3
import numpy as np
import os, sys, glob
import matplotlib.pyplot as plt

import vtk
from   tvtk import array_handler as ah

# Global Path
# global_path = ''
global_path = '/home/USADR/ac130084/CIFS/sph-simulations-backup/'

# define colorbar
cm = plt.cm.get_cmap('plasma')

simulations_res = [
                 30,
                 50,
                 100
 ]

path_couette = '03_couette_flow'
logfiles_couette = [
    #'couette_flow_20_28_17_vs_5e-05_run.log',
    'couette_flow_30_38_17_vs_3.3333333333333335e-05_run.log',
    'couette_flow_50_58_17_vs_2e-05_run.log',
    'couette_flow_100_108_17_vs_1e-05_run.log'
]

for i in range(len(simulations_res)):
    conv_data = np.genfromtxt(f'{global_path}{path_couette}/{logfiles_couette[i]}', skip_header = 1)[:, -2]
    number_fparticles = np.genfromtxt(f'{global_path}{path_couette}/{logfiles_couette[i]}', skip_header = 1)[-1, -3]
    norm_conv_data = np.divide(conv_data, number_fparticles)
    plt.plot(norm_conv_data, label = f'{simulations_res[i]}')
    plt.legend()
plt.show()

vtufiles =  [
    #'couette_flow_20_28_17_vs_5e-05_run/couette_flow_20_28_17_vs_5e-05_run_20000.vtu',
    'couette_flow_30_38_17_vs_3.3333333333333335e-05_run/couette_flow_30_38_17_vs_3.3333333333333335e-05_run_200000.vtu',
    'couette_flow_50_58_17_vs_2e-05_run/couette_flow_50_58_17_vs_2e-05_run_200000.vtu',
    'couette_flow_100_108_17_vs_1e-05_run/couette_flow_100_108_17_vs_1e-05_run_200000.vtu'
]

for i in range(len(vtufiles)):

    # Initialize VTK reader
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(f'{global_path}{path_couette}/{vtufiles[i]}')
    reader.Update()

    # Get data
    data = reader.GetOutput()

    # SPH particle positions
    vtk_points = data.GetPoints().GetData()
    points = ah.vtk2array(vtk_points)
    #print(points)

    # SPH particle type ID
    vtk_typeID = data.GetPointData().GetArray('TypeId')
    typeid = ah.vtk2array(vtk_typeID)

    # SPH particle type ID
    vtk_pressure = data.GetPointData().GetArray('Pressure')
    pressure = ah.vtk2array(vtk_pressure)

    # SPH particle type ID
    vtk_velocityx = data.GetPointData().GetArray('Velocity x')
    velocityx = ah.vtk2array(vtk_velocityx)

    # SPH particle type ID
    vtk_velocityy = data.GetPointData().GetArray('Velocity y')
    velocityy = ah.vtk2array(vtk_velocityy)

    tdata_points_y = []
    tdata_points_z = []
    tdata_vels = []
    for k in range(len(points)):
        if points[k, 0] < 0.0001  and points[k, 0] > -0.0001:
            tdata_points_y.append(points[k, 0])
            tdata_points_z.append(points[k, 1])
            tdata_vels.append(velocityx[k])

    im = plt.scatter(tdata_points_y, tdata_points_z, c = tdata_vels, s = 20, cmap = cm)
    plt.colorbar(im)
    plt.show()
    plt.close()
    
    plt.scatter(tdata_points_z, tdata_vels, c = tdata_vels, s = 20, cmap = cm)
    plt.show()
    plt.close()