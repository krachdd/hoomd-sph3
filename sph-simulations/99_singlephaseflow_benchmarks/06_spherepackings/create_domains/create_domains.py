#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:14:37 2023

@author: ac126015
"""

## HEADER ## 
import numpy as np
import os, sys, glob
from multiprocessing import Process

# Import py modules to handle .vtu files
sys.path.append('/scratch/local1/krach/data/00_processing/daves_py_modules/')
sys.path.append('/home/david/Arbeit/data/00_processing/daves_py_modules/')
import create_raw_files
import prepare_simulations
from mpi4py import MPI

# EXECUTE ON 8 CORES
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


types = ['bcc', 'fcc', 'sc']

porosities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# porosities = [0.1, 0.9]
resolutions = [100, 60]
# resolutions = [60]

lporosity = porosities[rank]


for i in range(len(types)):
    for k in range(len(resolutions)):
        lref = 1e-03
        print(rank, types[i], resolutions[k], lporosity)
        fn = create_raw_files.create_spherepacking( types[i], lporosity, lref, resolutions[k] , nx = 1, ny = 1, nz = 1)
        ifn = prepare_simulations.create_sph_input_file(fn, 'WendlandC4',  
                                lref/resolutions[k], resolutions[k], resolutions[k], resolutions[k], 
                                lporosity)
        if resolutions[k] == 60:
            cores = 64
        else:
            cores = 128
        prepare_simulations.create_slurm_file_sph(f'sp_{types[i]}_{lporosity}_{resolutions[k]}', cores, 1, './run_spherepacking_TV.py', ifn, ifn, partition = 'cpu', time = '24:00:00')
        