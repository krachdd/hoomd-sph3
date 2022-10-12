#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import os, glob
fs = [0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 50.0, 100.]
for f in range(len(fs)):
    # print(fs[f])
    os.system('mpirun -np 8 sphere_in_box_run.py {0}'.format(fs[f]))
