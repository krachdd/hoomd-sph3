# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""SPH models."""

import math
import numpy
import hoomd
import os
import inspect
import hoomd.sph
from hoomd     import _hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch
from hoomd.sph.force import Force

class _ModelBaseClass(Force):
    r"""
    Base class for models.
    
    Warning:
        This class should not be instantiated by users. The class can be used
        for `isinstance` or `issubclass` checks.
    """

    def __init__(self):
        super().__init__()

    def _attach(self):
        # check that some angles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No angles are defined.\n")

        # create the c++ mirror class
        if isinstance(self._simulation.device, hoomd.device.CPU):
            cpp_cls = getattr(_sph, self._cpp_class_name)
        else:
            cpp_cls = getattr(_sph, self._cpp_class_name + "GPU")

        self._cpp_obj = cpp_cls(self._simulation.state._cpp_sys_def)

        super()._attach()

    