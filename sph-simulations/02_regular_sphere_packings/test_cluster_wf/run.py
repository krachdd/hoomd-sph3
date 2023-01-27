#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import os, glob

for file in glob.glob("input*N60*.txt"):
    os.system('mpirun -np 4 run_spherepacking.py {0}'.format(file))
