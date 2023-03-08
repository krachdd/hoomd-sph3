#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import os, glob

for file in glob.glob("input*.txt"):
    os.system('mpirun -np 1 create_sc.py {0}'.format(file))
