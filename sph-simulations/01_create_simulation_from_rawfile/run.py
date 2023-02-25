#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import os, glob

for file in glob.glob("input*100*.txt"):
    os.system('mpirun -np 1 create_gsd_from_raw_singlecore.py {0}'.format(file))
