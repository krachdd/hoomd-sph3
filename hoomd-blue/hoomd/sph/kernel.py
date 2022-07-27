R""" SPH kernel classes
"""

import math
import numpy
import hoomd
import hoomd.sph
from hoomd     import _hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch

## \internal
# \brief Base class for smoothing kernel function classes
#
class _SmoothingKernel(hoomd.meta._metadata):
    ## \internal
    # \brief Constructs the kernel class
    #
    # Initializes the cpp_smoothingkernel.
    def __init__(self, name=None):
        # check if initialization has occurred
        if not hoomd.init.is_initialized():
            hoomd.context.msg.error("Cannot create a smoothing kernel before initialization\n");
            raise RuntimeError('Error creating smoothing kernel');

        self.kappa = 0;
        self.cpp_smoothingkernel = None;
        self.nlist = None;
        self.cpp_nlist = None;

        # Allow kernel class to store a name.
        if name is None:
            self.name = "";
        else:
            self.name="_" + name;

        # base class constructor
        hoomd.meta._metadata.__init__(self)

    def check_initialization(self):
        # check that we have been initialized properly
        if self.cpp_smoothingkernel is None:
            hoomd.context.msg.error('Bug in hoomd_script: cpp_smoothingkernel not set, please report\n');
            raise RuntimeError();

    def getKernelKappa(self):
        return self.cpp_smoothingkernel.getKernelKappa()

    def EvalKernel(self,h,rij):
        return self.cpp_smoothingkernel.EvalKernel(h,rij)

    def EvalKernelDerivative(self,h,rij):
        return self.cpp_smoothingkernel.EvalKernelDerivative(h,rij)

    def setNeighborList(self,nlist):
        # Neighbor list
        self.nlist = nlist

        # Set kernel scaling factor in neighbor list class
        self.nlist.cpp_nlist.setKernelFactor(self.kappa)

        # Set nlist in SmoothingKernel class
        #self.cpp_smoothingkernel.setNeighborList(self.nlist.cpp_nlist)
        
class WendlandC2(_SmoothingKernel):
    R""" Wendland C2 Kernel
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # Initialize base class
        _SmoothingKernel.__init__(self, "WendlandC2");

        # Kernel scaling parameter
        self.kappa = 2.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.WendlandC2();

    def OptimalH(self):
        return 1.7

    def Kappa(self):
        return self.kappa

class WendlandC4(_SmoothingKernel):
    R""" Wendland C4 Kernel
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # Initialize base class
        _SmoothingKernel.__init__(self, "WendlandC4");

        # Kernel scaling parameter
        self.kappa = 2.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.WendlandC4();

    def OptimalH(self):
        return 1.7

    def Kappa(self):
        return self.kappa
        
class WendlandC6(_SmoothingKernel):
    R""" Wendland C6 Kernel
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # Initialize base class
        _SmoothingKernel.__init__(self, "WendlandC6");

        # Kernel scaling parameter
        self.kappa = 2.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.WendlandC6();

    def OptimalH(self):
        return 1.7

    def Kappa(self):
        return self.kappa
               
class Quintic(_SmoothingKernel):
    R""" Quintic Kernel
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # Initialize base class
        _SmoothingKernel.__init__(self, "Quintic");

        # Kernel scaling parameter
        self.kappa = 3.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.Quintic();

    def OptimalH(self):
        return 1.45

    def Kappa(self):
        return self.kappa
        
class CubicSpline(_SmoothingKernel):
    R""" Cubic Spline Kernel
    """
    def __init__(self):
        hoomd.util.print_status_line();
        # Initialize base class
        _SmoothingKernel.__init__(self, "Quintic");

        # Kernel scaling parameter
        self.kappa = 2.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.Quintic();

    def OptimalH(self):
        return 1.7

    def Kappa(self):
        return self.kappa
        
            
# Dicts
Kernels  = {'WendlandC2':WendlandC2,'WendlandC4':WendlandC4,'WendlandC6':WendlandC6,'Quintic':Quintic,'CubicSpline':CubicSpline}
OptimalH = {'WendlandC2':1.7,'WendlandC4':1.7,'WendlandC6':1.7,'Quintic':1.45,'CubicSpline':1.7}
Kappa    = {'WendlandC2':2.0,'WendlandC4':2.0,'WendlandC6':2.0,'Quintic':3.0,'CubicSpline':2.0}
