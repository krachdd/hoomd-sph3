# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""SPH kernel classes."""

import hoomd 
import hoomd.sph
import hoomd.nsearch 
from hoomd import _hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch
import numpy
import math


class _SmoothingKernel(_HOOMDBaseObject):
    r"""
    Base class for smoothing kernel function classes
    """

    def __init__(self):
        self._in_context_manager = False

        self.kappa = 0;
        self.cpp_smoothingkernel = None;
        self.nlist = None;
        self.cpp_nlist = None;

        # Allow kernel class to store a name.
        if name is None:
            self.name = "";
        else:
            self.name="_" + name;

    def check_initialization(self):
        # check that we have been initialized properly
        if self.cpp_smoothingkernel is None:
            # hoomd.context.msg.error('Bug in hoomd_script: cpp_smoothingkernel not set, please report\n');
            raise RuntimeError('Bug in hoomd_script: cpp_smoothingkernel not set, please report\n');

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



class WendlandC2(_SmoothingKernel):
    R""" Wendland C2 Kernel
    """
    def __init__(self):
        # hoomd.util.print_status_line();
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
        # hoomd.util.print_status_line();
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
        # hoomd.util.print_status_line();
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
        # hoomd.util.print_status_line();
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
        # hoomd.util.print_status_line();
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
