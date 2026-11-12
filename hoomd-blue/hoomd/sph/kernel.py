# Copyright (c) 2025-2026 David Krach, Daniel Rostan.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de

r"""Smoothing kernels.

Smoothing kernels weight the contribution of neighboring particles in every SPH
summation. A kernel has compact support of radius :math:`\kappa h`, where
:math:`h` is the smoothing length.

The module also provides lookup tables keyed by the kernel class name:

* ``Kernels``: the kernel class.
* ``Kappa``: the support factor :math:`\kappa`.
* ``OptimalH``: the recommended ratio :math:`h / \Delta x` of smoothing length
  to initial particle spacing.

Example::

    kernel = hoomd.sph.kernel.WendlandC4()
    slength = hoomd.sph.kernel.OptimalH["WendlandC4"] * dx
    r_cut = hoomd.sph.kernel.Kappa["WendlandC4"] * slength
"""

import hoomd 
import hoomd.sph
import hoomd.nsearch 
from hoomd import _hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch
from hoomd.operation import _HOOMDBaseObject

import numpy
import math


class _SmoothingKernel(_HOOMDBaseObject):
    r"""
    Base class for smoothing kernel function classes
    """

    def __init__(self, name = None):
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
            # hoomd.context.msg.error("Bug in hoomd_script: cpp_smoothingkernel not set, please report\n");
            raise RuntimeError("Bug in hoomd_script: cpp_smoothingkernel not set, please report\n");

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
        self.nlist._cpp_obj.setKernelFactor(self.kappa)



class WendlandC2(_SmoothingKernel):
    R"""Wendland :math:`C^2` smoothing kernel.

    .. math::

        W(q) = \frac{\alpha_d}{h^d}
               \left(1 - \frac{q}{2}\right)^4 (1 + 2q),
        \quad q = r_{ij}/h,\; q \in [0, 2]

    where :math:`\alpha_d` is the normalisation constant in :math:`d` spatial
    dimensions and :math:`h` is the smoothing length.  The support radius
    is :math:`r_\mathrm{cut} = \kappa h` with :math:`\kappa = 2`.
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
    R"""Wendland :math:`C^4` smoothing kernel.

    .. math::

        W(q) = \frac{\alpha_d}{h^d}
               \left(1 - \frac{q}{2}\right)^6
               \left(1 + 3q + \frac{35}{12}q^2\right),
        \quad q = r_{ij}/h,\; q \in [0, 2]

    The support radius is :math:`r_\mathrm{cut} = \kappa h` with
    :math:`\kappa = 2`.
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
    R"""Wendland :math:`C^6` smoothing kernel.

    .. math::

        W(q) = \frac{\alpha_d}{h^d}
               \left(1 - \frac{q}{2}\right)^8
               \left(1 + 4q + \frac{25}{4}q^2 + 4q^3\right),
        \quad q = r_{ij}/h,\; q \in [0, 2]

    The support radius is :math:`r_\mathrm{cut} = \kappa h` with
    :math:`\kappa = 2`.
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
    R"""Quintic (Morris) smoothing kernel.

    .. math::

        W(q) = \frac{\alpha_d}{h^d}
               \begin{cases}
               (3-q)^5 - 6(2-q)^5 + 15(1-q)^5, & 0 \leq q < 1 \\
               (3-q)^5 - 6(2-q)^5,              & 1 \leq q < 2 \\
               (3-q)^5,                          & 2 \leq q < 3 \\
               0,                                & q \geq 3
               \end{cases}

    where :math:`q = r_{ij}/h`.  The support radius is
    :math:`r_\mathrm{cut} = \kappa h` with :math:`\kappa = 3`.
    """
    def __init__(self):
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
    R"""Cubic spline (M4) smoothing kernel.

    .. math::

        W(q) = \frac{\alpha_d}{h^d}
               \begin{cases}
               \frac{2}{3} - q^2 + \frac{1}{2}q^3, & 0 \leq q < 1 \\
               \frac{1}{6}(2 - q)^3,                & 1 \leq q < 2 \\
               0,                                    & q \geq 2
               \end{cases}

    where :math:`q = r_{ij}/h`.  The support radius is
    :math:`r_\mathrm{cut} = \kappa h` with :math:`\kappa = 2`.
    """
    def __init__(self):
        # Initialize base class
        _SmoothingKernel.__init__(self, "CubicSpline");

        # Kernel scaling parameter
        self.kappa = 2.0

        # create the c++ mirror class
        self.cpp_smoothingkernel = _sph.CubicSpline();

    def OptimalH(self):
        return 1.7

    def Kappa(self):
        return self.kappa
        
            
# Dicts
Kernels  = {"WendlandC2":WendlandC2,"WendlandC4":WendlandC4,"WendlandC6":WendlandC6,"Quintic":Quintic,"CubicSpline":CubicSpline}
OptimalH = {"WendlandC2":1.7,"WendlandC4":1.7,"WendlandC6":1.7,"Quintic":1.45,"CubicSpline":1.7}
Kappa    = {"WendlandC2":2.0,"WendlandC4":2.0,"WendlandC6":2.0,"Quintic":3.0,"CubicSpline":2.0}


__all__ = [
    "WendlandC2",
    "WendlandC4",
    "WendlandC6",
    "Quintic",
    "CubicSpline",
]