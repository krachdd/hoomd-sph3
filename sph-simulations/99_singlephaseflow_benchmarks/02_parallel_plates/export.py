#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""
# ----- HEADER -----------------------------------------------
import hoomd
from hoomd import *
from hoomd import sph
from hoomd.sph import _sph
import numpy as np
import math
# import itertools
from datetime import datetime
import export_gsd2vtu 
import read_input_fromtxt
import delete_solids_initial_timestep
import sys, os

import gsd.hoomd

filename = 'poiseuille_flow_20_28_28_vs_0.005_run.gsd'

export_gsd2vtu.export_spf(filename)