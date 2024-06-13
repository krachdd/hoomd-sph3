#!/usr/bin/env python3

"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

import os, glob

for file in glob.glob("input*.txt"):
    os.system('python3 create_gsd_from_raw_singlecore.py {0}'.format(file))
