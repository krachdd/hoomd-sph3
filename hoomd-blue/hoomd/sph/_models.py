# Copyright (c) 2009-2016 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

R""" sph - Smoothed-Particle Hydrodynamics models.
"""

import math
import numpy
import hoomd
import os
import inspect
import hoomd.sph
from hoomd     import _hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch
from hoomd.sph.force import _force

class _SPHBaseClass(_force):
    R""""
    """
    def __init__(self, kernel, eos, nlist):
        hoomd.util.print_status_line();

        # Initialize the base class
        hoomd.sph.force._force.__init__(self);

        # Private attributes
        self.kernel     = kernel
        self.eos        = eos
        self.nlist      = nlist
        self.accel_set  = False
        self.params_set = False
        self.gx         = 0.0
        self.gy         = 0.0
        self.gz         = 0.0

        # Detect maximum cut-off radius for nlist parameterization
        kappa     = self.kernel.Kappa()
        pdata     = hoomd.context.current.system_definition.getParticleData()
        globalN   = pdata.getNGlobal()

        maxh      = pdata.getMaxSmoothingLength()
        self.rcut = kappa * maxh

        # Neighbor list
        self.nlist.subscribe(lambda:self.get_rcut())
        self.nlist.update_rcut()
        self.nlist.cpp_nlist.setStorageMode(_nsearch.NeighborList.storageMode.full);

        # Set neighbor list in kernel class
        self.kernel.setNeighborList(self.nlist)

        # C++ mirror classes
        self.cpp_force = None

    def update_coeffs(self):
        """Noop for this compute"""
        pass

    def get_rcut(self):

        # Go through the list of only the active particle types in the simulation
        ntypes = hoomd.context.current.system_definition.getParticleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(hoomd.context.current.system_definition.getParticleData().getNameByType(i));

        # update the rcut by pair type
        r_cut_dict = hoomd.nsearch.nlist.rcut();
        for i in range(0,ntypes):
            for j in range(i,ntypes):
                r_cut_dict.set_pair(type_list[i],type_list[j],self.rcut);

        return r_cut_dict;

    def setBodyAcceleration(self,gx,gy,gz,damp=0):
        self.accel_set = True
        self.check_initialization();
        self.gx   = gx.item() if isinstance(gx, numpy.generic) else gx
        self.gy   = gy.item() if isinstance(gy, numpy.generic) else gy
        self.gz   = gz.item() if isinstance(gz, numpy.generic) else gz
        self.damp = int(damp.item()) if isinstance(damp,numpy.generic) else int(damp)
        self.damp = abs(self.damp)
        self.cpp_force.setAcceleration(self.gx,self.gy,self.gz,self.damp)

class SinglePhaseFlow(_SPHBaseClass):
    R""" SinglePhaseFlow solver
    """
    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup=None,
                 solidgroup=None,
                 densitymethod='SUMMATION',
                 viscositymethod='HARMONICAVERAGE'):
        hoomd.util.print_status_line();

        # Initialize base class
        _SPHBaseClass.__init__(self, kernel, eos, nlist);

        # Set force name
        self.force_name = 'SinglePhaseFlow'

        # Check if group is associated
        if fluidgroup == None or solidgroup == None:
            raise ValueError("SinglePhaseFlow: fluidgroup and solidgroup can not be None")
        else:
            self.fluidgroup = fluidgroup
            self.solidgroup = solidgroup
            self.cpp_fluidgroup = fluidgroup.cpp_group
            self.cpp_solidgroup = solidgroup.cpp_group

        # C++ mirror class
        func_name = '_sph.SinglePF'
        if hoomd.context.exec_conf.isCUDAEnabled():
            func_name += 'GPU'
        func_name += '_' + EOS[self.eos.name] + '_' + Kernel[self.kernel.name]
        func_name += """(hoomd.context.current.system_definition,
                         self.kernel.cpp_smoothingkernel,
                         self.eos.cpp_stateequation,
                         self.nlist.cpp_nlist,
                         self.cpp_fluidgroup,
                         self.cpp_solidgroup,
                         self.DENSITYMETHODS[densitymethod],
                         self.VISCOSITYMETHODS[viscositymethod]
                        )"""
        self.cpp_force = eval(func_name)

        # Check if given smoothing lengths are equal, if so
        # set constant smoothing length flag in model class
        pdata     = hoomd.context.current.system_definition.getParticleData()
        globalN   = pdata.getNGlobal()

        self.consth = pdata.constSmoothingLength()

        if self.consth:
            self.maxh = pdata.getSmoothingLength(0)
            self.cpp_force.setConstSmoothingLength(self.maxh)
        else:
            self.maxh = pdata.getMaxSmoothingLength()

        # Add to computes so that communicator gets set
        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

    def set_params(self,mu):
        self.check_initialization();
        self.mu   = mu.item()   if isinstance(mu, numpy.generic)   else mu
        self.cpp_force.setParams(self.mu)
        self.params_set = True

    @property
    def densitymethod(self):
        self.check_initialization();
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.iteritems())
        return invD[self.cpp_force.getDensityMethod()]

    @densitymethod.setter
    def densitymethod(self, method):
        self.check_initialization();
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self.cpp_force.setDensityMethod(self.DENSITYMETHODS[method])

    @property
    def viscositymethod(self):
        self.check_initialization();
        # Invert key mapping
        invD = dict((v,k) for k, v in self.VISCOSITYMETHODS.iteritems())
        return invD[self.cpp_force.getViscosityMethod()]

    @viscositymethod.setter
    def viscositymethod(self, method):
        self.check_initialization();
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self.cpp_force.setViscosityMethod(self.VISCOSITYMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.check_initialization();
        self.alpha   = alpha.item()  if isinstance(alpha, numpy.generic)   else alpha
        self.beta    = beta.item()   if isinstance(beta, numpy.generic)   else beta
        self.cpp_force.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self.check_initialization();
        self.cpp_force.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.check_initialization();
        self.ddiff   = ddiff.item()   if isinstance(ddiff, numpy.generic)   else ddiff
        self.cpp_force.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self.check_initialization();
        self.cpp_force.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.check_initialization();
        self.shepardfreq   = shepardfreq.item()   if isinstance(shepardfreq, numpy.generic)   else shepardfreq
        self.cpp_force.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self.check_initialization();
        self.cpp_force.deactivateShepardRenormalization()

    def computeSolidForces(self):
        self.check_initialization();
        self.cpp_force.computeSolidForces()

    def compute_dt(self,LREF,UREF,DRHO=0.01,COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        UREF = numpy.abs(UREF)

        # Compute required quantities
        # Magnitude of body force
        if not self.accel_set:
            GMAG = 0.0
        else:
            GMAG = numpy.sqrt(self.gx**2+self.gy**2+self.gz**2)
        # Smoothing length
        H   = self.maxh
        # Viscosity
        MU  = self.mu
        # Rest density
        RHO0 = self.eos.RestDensity

        # Speed of sound
        # CFL condition
        C2_1 = UREF*UREF/DRHO
        # Gravity waves condition
        C2_2 = GMAG*LREF/DRHO
        # Fourier condition
        C2_3 = (MU*UREF)/(RHO0*LREF*DRHO)
        # Maximum speed of sound
        C = numpy.sqrt(numpy.max([C2_1,C2_2,C2_3]))
        self.eos.set_speedofsound(C)

        # CFL condition
        DT_1 = 0.25*H/C
        # Fourier condition
        DT_2 = (H*H*RHO0)/(8.0*MU)
        if GMAG > 0.0:
            # Gravity waves condition
            DT_3 = numpy.sqrt(H/(16.0*GMAG))
            return COURANT*numpy.min([DT_1,DT_2,DT_3])
        else:
            return COURANT*numpy.min([DT_1,DT_2])


# class TwoPhaseFlow(_SPHBaseClass):
#     R""" TwoPhaseFlow solver
#     """
#     DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
#                       'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

#     VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

#     def __init__(self,
#                  kernel,
#                  eos1,
#                  eos2,
#                  nlist,
#                  fluidgroup1=None,
#                  fluidgroup2=None,
#                  solidgroup=None,
#                  densitymethod='SUMMATION',
#                  viscositymethod='HARMONICAVERAGE'):
#         hoomd.util.print_status_line();

#         # Initialize base class
#         _SPHBaseClass.__init__(self, kernel, eos1, nlist);

#         # Set force name
#         self.force_name = 'TwoPhaseFlow'

#         # Set eos instances
#         self.eos1 = eos1
#         self.eos2 = eos2

#         # Check if group is associated
#         if fluidgroup1 == None or fluidgroup2 == None or solidgroup == None:
#             raise ValueError("TwoPhaseFlow: fluidgroups and solidgroup can not be None")
#         else:
#             self.fluidgroup1 = fluidgroup1
#             self.fluidgroup2 = fluidgroup2
#             self.solidgroup  = solidgroup
#             self.cpp_fluidgroup1 = fluidgroup1.cpp_group
#             self.cpp_fluidgroup2 = fluidgroup2.cpp_group
#             self.cpp_solidgroup  = solidgroup.cpp_group

#         # C++ mirror class
#         func_name = '_sph.TwoPF'
#         if hoomd.context.exec_conf.isCUDAEnabled():
#             func_name += 'GPU'
#         func_name += '_' + EOS[self.eos1.name] + EOS[self.eos2.name] + '_' + Kernel[self.kernel.name]
#         func_name += """(hoomd.context.current.system_definition,
#                          self.kernel.cpp_smoothingkernel,
#                          self.eos1.cpp_stateequation,
#                          self.eos2.cpp_stateequation,
#                          self.nlist.cpp_nlist,
#                          self.cpp_fluidgroup1,
#                          self.cpp_fluidgroup2,
#                          self.cpp_solidgroup,
#                          self.DENSITYMETHODS[densitymethod],
#                          self.VISCOSITYMETHODS[viscositymethod]
#                         )"""
#         self.cpp_force = eval(func_name)

#         # Check if given smoothing lengths are equal, if so
#         # set constant smoothing length flag in model class
#         pdata     = hoomd.context.current.system_definition.getParticleData()
#         globalN   = pdata.getNGlobal()

#         self.consth = pdata.constSmoothingLength()

#         if self.consth:
#             self.maxh = pdata.getSmoothingLength(0)
#             self.cpp_force.setConstSmoothingLength(self.maxh)
#         else:
#             self.maxh = pdata.getMaxSmoothingLength()

#         # Add to computes so that communicator gets set
#         hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

#     def set_params(self,mu1,mu2,sigma12,omega):
#         self.check_initialization();
#         self.mu1     = mu1.item()     if isinstance(mu1, numpy.generic)     else mu1
#         self.mu2     = mu2.item()     if isinstance(mu2, numpy.generic)     else mu2
#         self.sigma12 = sigma12.item() if isinstance(sigma12, numpy.generic) else sigma12
#         self.omega   = omega.item()   if isinstance(omega, numpy.generic)   else omega
#         self.cpp_force.setParams(self.mu1,self.mu2,self.sigma12,self.omega)
#         self.params_set = True

#     @property
#     def densitymethod(self):
#         self.check_initialization();
#         # Invert key mapping
#         invD = dict((v,k) for k, v in self.DENSITYMETHODS.iteritems())
#         return invD[self.cpp_force.getDensityMethod()]

#     @densitymethod.setter
#     def densitymethod(self, method):
#         self.check_initialization();
#         if method not in self.DENSITYMETHODS:
#             raise ValueError("Undefined DensityMethod.")
#         self.cpp_force.setDensityMethod(self.DENSITYMETHODS[method])

#     @property
#     def viscositymethod(self):
#         self.check_initialization();
#         # Invert key mapping
#         invD = dict((v,k) for k, v in self.VISCOSITYMETHODS.iteritems())
#         return invD[self.cpp_force.getViscosityMethod()]

#     @viscositymethod.setter
#     def viscositymethod(self, method):
#         self.check_initialization();
#         if method not in self.VISCOSITYMETHODS:
#             raise ValueError("Undefined ViscosityMethod.")
#         self.cpp_force.setViscosityMethod(self.VISCOSITYMETHODS[method])

#     def activateArtificialViscosity(self, alpha, beta):
#         self.check_initialization();
#         self.alpha   = alpha.item()  if isinstance(alpha, numpy.generic)   else alpha
#         self.beta    = beta.item()   if isinstance(beta, numpy.generic)   else beta
#         self.cpp_force.activateArtificialViscosity(alpha, beta)

#     def deactivateArtificialViscosity(self):
#         self.check_initialization();
#         self.cpp_force.deactivateArtificialViscosity()

#     def activateDensityDiffusion(self, ddiff):
#         self.check_initialization();
#         self.ddiff   = ddiff.item()   if isinstance(ddiff, numpy.generic)   else ddiff
#         self.cpp_force.activateDensityDiffusion(ddiff)

#     def deactivateDensityDiffusion(self):
#         self.check_initialization();
#         self.cpp_force.deactivateDensityDiffusion()

#     def compute_dt(self,LREF,UREF,DRHO=0.01,COURANT=0.25):
#         # Input sanity
#         if LREF == 0.0:
#             raise ValueError('Reference length LREF may not be zero.')
#         if DRHO == 0.0:
#             raise ValueError('Maximum density variation DRHO may not be zero.')
#         UREF = numpy.abs(UREF)

#         # Compute required quantities
#         # Magnitude of body force
#         if not self.accel_set:
#             GMAG = 0.0
#         else:
#             GMAG = numpy.sqrt(self.gx**2+self.gy**2+self.gz**2)
#         # Smoothing length
#         H    = self.maxh
#         # Viscosity
#         MU1  = self.mu1
#         MU2  = self.mu2
#         # Rest density
#         RHO01 = self.eos1.RestDensity
#         RHO02 = self.eos2.RestDensity
#         # Surface tension
#         SIGMA = self.sigma12

#         # Speed of sound
#         # CFL condition
#         C2_1 = UREF*UREF/DRHO
#         # Gravity waves condition
#         C2_2 = GMAG*LREF/DRHO
#         # Surface wave condition
#         C2_31 = SIGMA/(RHO01*LREF*DRHO)
#         C2_32 = SIGMA/(RHO02*LREF*DRHO)
#         # Fourier condition
#         C2_41 = (MU1*UREF)/(RHO01*LREF*DRHO)
#         C2_42 = (MU2*UREF)/(RHO02*LREF*DRHO)

#         # Maximum speed of sound
#         C_1 = numpy.sqrt(numpy.max([C2_1,C2_2,C2_31,C2_41]))
#         C_2 = numpy.sqrt(numpy.max([C2_1,C2_2,C2_32,C2_42]))
#         self.eos1.set_speedofsound(C_1)
#         self.eos2.set_speedofsound(C_2)

#         # CFL condition
#         DT_1 = 0.25*H/numpy.max([C_1,C_2])
#         # Fourier condition
#         DT_2 = (H*H*numpy.min([RHO01,RHO02]))/(8.0*numpy.max([MU1,MU2]))
#         # Surface wave condition
#         if SIGMA > 0.0:
#             DT_3 = numpy.sqrt((H*H*H*numpy.min([RHO01,RHO02]))/(32.0*numpy.pi*SIGMA))
#         else:
#             DT_3 = DT_1*2
#         # Gravity waves condition
#         if GMAG > 0.0:
#             DT_4 = numpy.sqrt(H/(16.0*GMAG))
#         else:
#             DT_4 = DT_1*3

#         return COURANT*numpy.min([DT_1,DT_2,DT_3,DT_4])


# Dicts
Kernel = {'_WendlandC2':'WC2','_WendlandC4':'WC4','_WendlandC6':'WC6','_Quintic':'Q','_CubicSplibe':'CS'}
EOS = {'_Linear':'L','_Tait':'T'}
