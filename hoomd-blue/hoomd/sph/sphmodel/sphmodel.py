"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de


SPH Momentum interaction forces."""

import copy
import warnings

import hoomd
from hoomd.sph import _sph
from hoomd.nsearch import _nsearch
from hoomd.sph import force
from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
import numpy as np
from hoomd.data.typeconverter import OnlyFrom, nonnegative_real
from itertools import combinations_with_replacement

class SPHModel(force.Force):
    r""" Base class for all SPH Models

    """

    # Module where the C++ class is defined. Reassign this when developing an
    # external plugin.
    _ext_module = _sph
    _ext_module_nlist = _nsearch

    def __init__(self, kernel, eos, nlist):
        super().__init__()



        # default exclusions
        params = ParameterDict(accel_set = bool(False),
                               params_set = bool(False),
                               gx=float(0.0),
                               gy=float(0.0),
                               gz=float(0.0),
                               damp=int(0))
        self._param_dict.update(params)

        self.kernel     = kernel
        self.eos        = eos

        type_params = []
        self._extend_typeparam(type_params)
        self._param_dict.update(
            ParameterDict(nlist=hoomd.nsearch.nlist.NeighborList))
        self.nlist = nlist

    def _attach_hook(self):
        # check that some Particles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")

        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        self._cpp_baseclass_name = "SPHBaseClass" + "_" + Kernel[self.kernel.name] + "_" + EOS[self.eos.name]
        base_cls = getattr(_sph, self._cpp_baseclass_name)
        self._cpp_base_obj = base_cls(self._simulation.state._cpp_sys_def, self.kernel.cpp_smoothingkernel,
                                 self.eos.cpp_stateequation, self.nlist._cpp_obj)



    def _setattr_param(self, attr, value):
        if attr == "nlist":
            self._nlist_setter(value)
            return
        super()._setattr_param(attr, value)

    def _nlist_setter(self, new_nlist):
        if new_nlist is self.nlist:
            return
        if self._attached:
            raise RuntimeError("nlist cannot be set after scheduling.")
        self._param_dict._dict["nlist"] = new_nlist

    def get_rcut(self):

        # Go through the list of only the active particle types in the simulation
        # ntypes = hoomd.context.current.system_definition.getParticleData().getNTypes();
        ntypes = self._simulation.state._cpp_sys_def.getParticleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(self._simulation.state._cpp_sys_def.getParticleData().getNameByType(i));

        # update the rcut by pair type
        r_cut_dict = hoomd.nsearch.nlist.rcut();
        for i in range(0,ntypes):
            for j in range(i,ntypes):
                r_cut_dict.set_pair(type_list[i],type_list[j],self.rcut);

        return r_cut_dict;

    def get_typelist(self):
        # Go through the list of only the active particle types in the simulation
        ntypes = self._simulation.state._cpp_sys_def.getParticleData().getNTypes();
        type_list = [];
        for i in range(0,ntypes):
            type_list.append(self._simulation.state._cpp_sys_def.getParticleData().getNameByType(i));
        return type_list





class SinglePhaseFlow(SPHModel):
    R""" SinglePhaseFlow solver
    """
    DENSITYMETHODS = {"SUMMATION":_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      "CONTINUITY":_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {"HARMONICAVERAGE":_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter = None,
                 solidgroup_filter = None,
                 densitymethod=None,
                 viscositymethod="HARMONICAVERAGE"):

        super().__init__(kernel, eos, nlist)

        self._param_dict.update(ParameterDict(
                          densitymethod = densitymethod,
                          viscositymethod = viscositymethod, 
                          mu = float(0.0), 
                          artificialviscosity = bool(True), 
                          alpha = float(0.2),
                          beta = float(0.0),
                          densitydiffusion = bool(False),
                          ddiff = float(0.0),
                          shepardrenormanlization = bool(False),
                          densityreinitialization = bool(False),
                          shepardfreq = int(0),
                          densityreinitfreq = int(0),
                          compute_solid_forces = bool(False),
                          max_sl = float(0.0)
                          ))



        # self._state = self._simulation.state
        self._cpp_SPFclass_name = "SinglePF_" + Kernel[self.kernel.name] + "_" + EOS[self.eos.name]
        self.fluidgroup_filter = fluidgroup_filter
        self.solidgroup_filter = solidgroup_filter
        self.str_densitymethod = densitymethod
        self.str_viscositymethod = viscositymethod
        self.accel_set = False
        self.params_set = False

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")


    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spf_cls = getattr(_sph, self._cpp_SPFclass_name)
        else:
            print("GPU not implemented")

        # check that some Particles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")
        
        # This should never happen, but leaving it in case the logic for adding
        # missed some edge case.
        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        cpp_sys_def = self._simulation.state._cpp_sys_def
        cpp_fluidgroup  = self._simulation.state._get_group(self.fluidgroup_filter)
        cpp_solidgroup  = self._simulation.state._get_group(self.solidgroup_filter)
        cpp_kernel = self.kernel.cpp_smoothingkernel
        cpp_eos = self.eos.cpp_stateequation
        cpp_nlist =  self.nlist._cpp_obj

        # Set Kernel specific Kappa in cpp-Nlist
        self.kernel.setNeighborList(self.nlist)

        self._cpp_obj = spf_cls(cpp_sys_def, cpp_kernel, cpp_eos, cpp_nlist, cpp_fluidgroup, 
                                cpp_solidgroup, self.cpp_densitymethod, self.cpp_viscositymethod)

        # Set kernel parameters
        kappa = self.kernel.Kappa()
        mycpp_kappa = self.kernel.cpp_smoothingkernel.getKernelKappa()

        pdata = self._simulation.state._cpp_sys_def.getParticleData()
        globalN = pdata.getNGlobal()

        self.consth = pdata.constSmoothingLength()
        if self.consth:
            self.maxh = pdata.getSlength(0)
            if (self._simulation.device.communicator.rank == 0):
                print(f"Using constant Smoothing Length: {self.maxh}")

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print("Non-Constant Smoothing length")
        self.rcut = kappa * self.maxh


        # Set rcut in neigbour list
        self._param_dict.update(ParameterDict(
                          rcut = self.rcut, 
                          max_sl = self.maxh
                          ))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        # get all params in line
        self.mu = self._param_dict["mu"]
        self.artificialviscosity = self._param_dict["artificialviscosity"]
        self.alpha = self._param_dict["alpha"]
        self.beta = self._param_dict["beta"]
        self.densitydiffusion = self._param_dict["densitydiffusion"]
        self.ddiff = self._param_dict["ddiff"]
        self.shepardrenormanlization = self._param_dict["shepardrenormanlization"]
        self.densityreinitialization = self._param_dict["densityreinitialization"]
        self.shepardfreq = self._param_dict["shepardfreq"]
        self.densityreinitfreq = self._param_dict["densityreinitfreq"]
        self.compute_solid_forces = self._param_dict["compute_solid_forces"]

        self.set_params(self.mu)
        self.setdensitymethod(self.str_densitymethod)
        self.setviscositymethod(self.str_viscositymethod)
        
        if (self.artificialviscosity == True):
            self.activateArtificialViscosity(self.alpha, self.beta)
        else:
            self.deactivateArtificialViscosity()
        
        if (self.densitydiffusion == True):
            self.activateDensityDiffusion(self.ddiff)
        else:
            self.deactivateDensityDiffusion()
        
        if (self.shepardrenormanlization == True):
            self.activateShepardRenormalization(self.shepardfreq)
        else:
            self.deactivateShepardRenormalization()

        if (self.densityreinitialization == True):
            self.activateDensityReinitialization(self.densityreinitfreq)
        else:
            self.deactivateDensityReinitialization()
        
        if (self.compute_solid_forces == True):
            self.computeSolidForces()

        self.setrcut(self.rcut, self.get_typelist())

        self.setBodyAcceleration(self.gx, self.gy, self.gz, self.damp)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    def _detach_hook(self):
        self.nlist._detach()

    def set_params(self,mu):
        # self.mu   = mu.item()   if isinstance(mu, np.generic)   else mu
        self._cpp_obj.setParams(self.mu)
        self.params_set = True
        self._param_dict.__setattr__("params_set", True)

    # @rcut.setter
    def setrcut(self, rcut, types):
        if rcut <= 0.0:
            raise ValueError("Rcut has to be > 0.0.")
        for p in combinations_with_replacement(types, 2):
            self._cpp_obj.setRCut(p, rcut)

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    # @property
    def viscositymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.VISCOSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getViscosityMethod()]

    # @viscositymethod.setter
    def setviscositymethod(self, method):
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self._cpp_obj.setViscosityMethod(self.VISCOSITYMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.alpha   = alpha.item()  if isinstance(alpha, np.generic)   else alpha
        self.beta    = beta.item()   if isinstance(beta, np.generic)   else beta
        self._cpp_obj.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self._cpp_obj.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.ddiff   = ddiff.item()   if isinstance(ddiff, np.generic)   else ddiff
        self._cpp_obj.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self._cpp_obj.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.shepardfreq   = shepardfreq.item()   if isinstance(shepardfreq, np.generic)   else shepardfreq
        self._cpp_obj.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self._cpp_obj.deactivateShepardRenormalization()

    def activateDensityReinitialization(self, densityreinitfreq=20):
        self.densityreinitfreq   = densityreinitfreq.item()   if isinstance(densityreinitfreq, np.generic)   else densityreinitfreq
        self._cpp_obj.activateDensityReinitialization(int(densityreinitfreq))

    def deactivateDensityReinitialization(self):
        self._cpp_obj.deactivateDensityReinitialization()

    def activatePowerLaw(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw(float(K), float(n), float(mu_min))

    def activateCarreau(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian(self):
        self._cpp_obj.deactivateNonNewtonian()

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

    def setBodyAcceleration(self,gx,gy,gz,damp=0):
        self.accel_set = True
        self._param_dict.__setattr__("accel_set", True)
        # self.check_initialization();
        # self.gx   = gx.item() if isinstance(gx, np.generic) else gx
        # self.gy   = gy.item() if isinstance(gy, np.generic) else gy
        # self.gz   = gz.item() if isinstance(gz, np.generic) else gz
        # self.damp = int(damp.item()) if isinstance(damp,np.generic) else int(damp)
        self.damp = abs(self.damp)

        if ( self.gx == 0 and self.gy == 0 and self.gz == 0):
            if ( self._simulation.device.communicator.rank == 0 ):
                print(f"{self._cpp_SPFclass_name} does NOT use a body force!" )

        self._cpp_obj.setAcceleration(self.gx,self.gy,self.gz,self.damp)

    def get_speedofsound(self):
        return self.eos.SpeedOfSound

    def set_speedofsound(self, c):
        self.eos.set_speedofsound(c)

    def get_GMAG(self):
        # Magnitude of body force
        if (abs(self.gx) > 0.0 or abs(self.gy) > 0.0 or abs(self.gz) > 0.0):
            return  np.sqrt(self.gx**2+self.gy**2+self.gz**2)
        else:
            return 0.0

    def compute_speedofsound(self, LREF, UREF, DX, DRHO, H, MU, RHO0):
        # Input sanity
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")

        UREF = np.abs(UREF)

        C_a = []
        # $c_0^2 \geq U_\mathrm{ref}^2 / \delta\rho$  (CFL)
        C_a.append(UREF*UREF/DRHO)
        # $c_0^2 \geq g\,l_\mathrm{ref} / \delta\rho$  (gravity-wave)
        C_a.append(self.get_GMAG()*LREF/DRHO)
        # $c_0^2 \geq \mu U_\mathrm{ref} / (\rho_0 l_\mathrm{ref} \delta\rho)$  (Fourier)
        C_a.append((MU*UREF)/(RHO0*LREF*DRHO))
        # $c_0 = \sqrt{\max(c_0^2)}$

        C_a = np.asarray(C_a)
        conditions = ["CFL-condition", "Gravity_waves-condition", "Fourier-condition"]
        condition = [conditions[i] for i in np.where(C_a == C_a.max())[0]]
        C = np.sqrt(np.max(C_a))

        # Set speed of sound
        self.eos.set_speedofsound(C)

        return C, condition


    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")
        if H != self._param_dict["max_sl"]:
            raise ValueError("Given H not equal to stored H self._param_dict[max_sl]!")
        if MU != self._param_dict["mu"]:
            raise ValueError("Given MU not equal to stored MU self._param_dict[mu]!")
        if RHO0 != self.eos.RestDensity:
            raise ValueError("Given RHO0 not equal to stored RHO0 self.eos.RestDensity!")
        
        UREF = np.abs(UREF)

        C = self.get_speedofsound()

        DT_a = []
        # $\Delta t \leq \Delta x / c_0$  (CFL)
        DT_a.append(DX/C)
        # $\Delta t \leq \rho_0 \Delta x^2 / (8\mu)$  (Fourier/viscous)
        DT_a.append((DX*DX*RHO0)/(8.0*MU))

        if self.get_GMAG() > 0.0:
            # $\Delta t \leq \sqrt{h / (16 g)}$  (gravity-wave)
            DT_a.append(np.sqrt(H/(16.0*self.get_GMAG())))
        DT_a = np.asarray(DT_a)
        conditions = ["CFL-condition", "Fourier-condition", "Gravity_waves-condition"]
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        DT = COURANT * np.min(DT_a)

        return DT, condition


class SinglePhaseFlowGDGD(SPHModel):
    R"""General Density Gradient Driven Flow solver.

    Extends SinglePhaseFlow with a transported scalar field T (temperature or
    concentration) that drives buoyancy via one of two models:

    * **VRD** (``boussinesq=False``, default):
      Per-particle rest density $\rho_{0,i} = \rho_0 (1 - \beta (T_i - T_\mathrm{ref}))$.
      Buoyancy emerges implicitly from the pressure gradient.
      For SUMMATION density method, VRD pressures are computed on-the-fly
      in the pair loop.  For CONTINUITY, the VRD $\partial P/\partial\rho$ derivative is used
      in the $\dot{p}$ chain rule.

    * **Boussinesq** (``boussinesq=True``):
      Standard EOS with global $\rho_0$; explicit per-particle buoyancy correction
      $\Delta F_b = m \, g \, (-\beta (T_i - T_\mathrm{ref}))$ is added to the momentum equation.

    The scalar $T$ is stored in ``aux4.x``.  Its rate of change $\dot{T}$ is
    accumulated into ``ratedpe.z`` (by the scalar diffusion operator) and
    time-marched by the integrator in the same half-step scheme as density.

    Wall boundary temperatures are set by assigning ``aux4.x`` directly to
    solid particles (e.g. via ``sim.state.cpu_local_snapshot``).

    Parameters
    ----------
    kernel : SmoothingKernel
        Smoothing kernel object.
    eos : StateEquation
        Equation of state object (Tait or Linear).
    nlist : NeighborList
        Neighbour list object.
    fluidgroup_filter : hoomd.filter
        Filter selecting fluid particles.
    solidgroup_filter : hoomd.filter
        Filter selecting solid (wall) particles.
    densitymethod : str
        ``'SUMMATION'`` or ``'CONTINUITY'``.
    viscositymethod : str
        ``'HARMONICAVERAGE'`` (only option currently).
    kappa_s : float
        Scalar diffusivity [m²/s] (thermal diffusivity α = λ/(ρ·cₚ) or
        mass diffusivity D).  Default 0.
    beta_s : float
        Expansion coefficient [1/K] (thermal β or solutal βc).  Default 0.
    scalar_ref : float
        Reference scalar value T_ref (or c_ref).  Default 0.
    boussinesq : bool
        If True, Boussinesq approximation; if False, Variable Reference
        Density (VRD).  Default False.
    """

    DENSITYMETHODS = {"SUMMATION": _sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      "CONTINUITY": _sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {"HARMONICAVERAGE": _sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter=None,
                 solidgroup_filter=None,
                 densitymethod=None,
                 viscositymethod="HARMONICAVERAGE",
                 kappa_s=0.0,
                 beta_s=0.0,
                 scalar_ref=0.0,
                 boussinesq=False):

        super().__init__(kernel, eos, nlist)

        self._param_dict.update(ParameterDict(
            densitymethod=densitymethod,
            viscositymethod=viscositymethod,
            mu=float(0.0),
            artificialviscosity=bool(True),
            alpha=float(0.2),
            beta=float(0.0),
            densitydiffusion=bool(False),
            ddiff=float(0.0),
            shepardrenormanlization=bool(False),
            densityreinitialization=bool(False),
            shepardfreq=int(0),
            densityreinitfreq=int(0),
            compute_solid_forces=bool(False),
            max_sl=float(0.0),
            # GDGD-specific parameters
            kappa_s=float(kappa_s),
            beta_s=float(beta_s),
            scalar_ref=float(scalar_ref),
            boussinesq=bool(boussinesq),
        ))

        self._cpp_SPFclass_name = ("SinglePFGDGD"
                                   "_" + Kernel[self.kernel.name]
                                   + "_" + EOS[self.eos.name])
        self.fluidgroup_filter  = fluidgroup_filter
        self.solidgroup_filter  = solidgroup_filter
        self.str_densitymethod  = densitymethod
        self.str_viscositymethod= viscositymethod
        self.accel_set          = False
        self.params_set         = False
        self.gdgd_params_set    = False

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spf_cls = getattr(_sph, self._cpp_SPFclass_name)
        else:
            print("GPU not implemented")

        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")

        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        else:
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        cpp_sys_def    = self._simulation.state._cpp_sys_def
        cpp_fluidgroup = self._simulation.state._get_group(self.fluidgroup_filter)
        cpp_solidgroup = self._simulation.state._get_group(self.solidgroup_filter)
        cpp_kernel     = self.kernel.cpp_smoothingkernel
        cpp_eos        = self.eos.cpp_stateequation
        cpp_nlist      = self.nlist._cpp_obj

        self.kernel.setNeighborList(self.nlist)

        self._cpp_obj = spf_cls(cpp_sys_def, cpp_kernel, cpp_eos, cpp_nlist,
                                cpp_fluidgroup, cpp_solidgroup,
                                self.cpp_densitymethod, self.cpp_viscositymethod)

        kappa = self.kernel.Kappa()

        pdata   = self._simulation.state._cpp_sys_def.getParticleData()
        self.consth = pdata.constSmoothingLength()
        if self.consth:
            self.maxh = pdata.getSlength(0)
            if self._simulation.device.communicator.rank == 0:
                print(f"Using constant Smoothing Length: {self.maxh}")
            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else:
            self.maxh = pdata.getMaxSmoothingLength()
            if self._simulation.device.communicator.rank == 0:
                print("Non-Constant Smoothing length")
        self.rcut = kappa * self.maxh

        self._param_dict.update(ParameterDict(rcut=self.rcut, max_sl=self.maxh))

        # Reload parameters (may have been set before attach)
        self.str_densitymethod  = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod= self._param_dict._dict["viscositymethod"]

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        # Apply all stored parameters
        self.mu                   = self._param_dict["mu"]
        self.artificialviscosity  = self._param_dict["artificialviscosity"]
        self.alpha                = self._param_dict["alpha"]
        self.beta                 = self._param_dict["beta"]
        self.densitydiffusion     = self._param_dict["densitydiffusion"]
        self.ddiff                = self._param_dict["ddiff"]
        self.shepardrenormanlization = self._param_dict["shepardrenormanlization"]
        self.densityreinitialization = self._param_dict["densityreinitialization"]
        self.shepardfreq          = self._param_dict["shepardfreq"]
        self.densityreinitfreq    = self._param_dict["densityreinitfreq"]
        self.compute_solid_forces = self._param_dict["compute_solid_forces"]
        kappa_s_val               = self._param_dict["kappa_s"]
        beta_s_val                = self._param_dict["beta_s"]
        scalar_ref_val            = self._param_dict["scalar_ref"]
        boussinesq_val            = self._param_dict["boussinesq"]

        self.set_params(self.mu)
        self.setdensitymethod(self.str_densitymethod)
        self.setviscositymethod(self.str_viscositymethod)

        if self.artificialviscosity:
            self.activateArtificialViscosity(self.alpha, self.beta)
        else:
            self.deactivateArtificialViscosity()

        if self.densitydiffusion:
            self.activateDensityDiffusion(self.ddiff)
        else:
            self.deactivateDensityDiffusion()

        if self.shepardrenormanlization:
            self.activateShepardRenormalization(self.shepardfreq)
        else:
            self.deactivateShepardRenormalization()

        if self.densityreinitialization:
            self.activateDensityReinitialization(self.densityreinitfreq)
        else:
            self.deactivateDensityReinitialization()

        if self.compute_solid_forces:
            self.computeSolidForces()

        # Set GDGD-specific parameters in the C++ object
        self.setGDGDParams(kappa_s_val, beta_s_val, scalar_ref_val, boussinesq_val)

        self.setrcut(self.rcut, self.get_typelist())
        self.setBodyAcceleration(self.gx, self.gy, self.gz, self.damp)

        super()._attach_hook()

    def _detach_hook(self):
        self.nlist._detach()

    def set_params(self, mu):
        self._cpp_obj.setParams(self.mu)
        self.params_set = True
        self._param_dict.__setattr__("params_set", True)

    def setGDGDParams(self, kappa_s, beta_s, scalar_ref, boussinesq=False):
        """Set GDGD scalar-transport parameters.

        Parameters
        ----------
        kappa_s : float
            Scalar diffusivity [m²/s].
        beta_s : float
            Expansion coefficient [1/K].
        scalar_ref : float
            Reference scalar value T_ref.
        boussinesq : bool
            If True, Boussinesq approximation.  If False (default), VRD.
        """
        self.gdgd_params_set = True
        self._cpp_obj.setGDGDParams(float(kappa_s), float(beta_s),
                                     float(scalar_ref), bool(boussinesq))
        self._param_dict.__setattr__("kappa_s",    float(kappa_s))
        self._param_dict.__setattr__("beta_s",     float(beta_s))
        self._param_dict.__setattr__("scalar_ref", float(scalar_ref))
        self._param_dict.__setattr__("boussinesq", bool(boussinesq))

    # ── Methods inherited from SinglePhaseFlowTV / SinglePhaseFlow ────────────
    # These delegate directly to the C++ object via the same pybind11 bindings.

    def setrcut(self, rcut, types):
        if rcut <= 0.0:
            raise ValueError("Rcut has to be > 0.0.")
        for p in list(combinations_with_replacement(types, 2)):
            self._cpp_obj.setRCut(p, rcut)

    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    def setviscositymethod(self, method):
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self._cpp_obj.setViscosityMethod(self.VISCOSITYMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.alpha = float(alpha)
        self.beta  = float(beta)
        self._cpp_obj.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self._cpp_obj.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.ddiff = float(ddiff)
        self._cpp_obj.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self._cpp_obj.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.shepardfreq = int(shepardfreq)
        self._cpp_obj.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self._cpp_obj.deactivateShepardRenormalization()

    def activateDensityReinitialization(self, densityreinitfreq=20):
        self.densityreinitfreq = int(densityreinitfreq)
        self._cpp_obj.activateDensityReinitialization(int(densityreinitfreq))

    def deactivateDensityReinitialization(self):
        self._cpp_obj.deactivateDensityReinitialization()

    def activatePowerLaw(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw(float(K), float(n), float(mu_min))

    def activateCarreau(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian(self):
        self._cpp_obj.deactivateNonNewtonian()

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

    def setBodyAcceleration(self, gx, gy, gz, damp=0):
        self.accel_set = True
        self._param_dict.__setattr__("accel_set", True)
        self.damp = abs(self.damp)
        if self.gx == 0 and self.gy == 0 and self.gz == 0:
            if self._simulation.device.communicator.rank == 0:
                print(f"{self._cpp_SPFclass_name} does NOT use a body force!")
        self._cpp_obj.setAcceleration(self.gx, self.gy, self.gz, self.damp)

    def get_speedofsound(self):
        return self.eos.SpeedOfSound

    def set_speedofsound(self, c):
        self.eos.set_speedofsound(c)

    def get_GMAG(self):
        if abs(self.gx) > 0.0 or abs(self.gy) > 0.0 or abs(self.gz) > 0.0:
            return np.sqrt(self.gx**2 + self.gy**2 + self.gz**2)
        return 0.0

    def compute_speedofsound(self, LREF, UREF, DX, DRHO, H, MU, RHO0):
        """Delegate to SinglePhaseFlow compute_speedofsound logic."""
        # Reuse computation from parent — instantiate a temporary SPF for the
        # same EOS/kernel, or duplicate the formula directly.
        # For simplicity, call the equivalent formula (same as SinglePhaseFlow).
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            DRHO = 0.01
        c_list = []
        cond_list = []
        c_list.append(10.0 * UREF)
        cond_list.append("10*UREF")
        if MU > 0.0:
            c_list.append(np.sqrt(10.0 * DRHO * MU * UREF / (RHO0 * H * H)))
            cond_list.append("viscous-condition")
        c_array = np.asarray(c_list)
        idx = np.argmax(c_array)
        return c_array[idx], cond_list[idx]

    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.2):
        """Estimate maximum stable timestep."""
        C = self.eos.SpeedOfSound
        DT_a = []
        DT_a.append(COURANT * H / C)
        DT_a.append((DX * DX * RHO0) / (8.0 * MU))
        if self.get_GMAG() > 0.0:
            DT_a.append(np.sqrt(H / (16.0 * self.get_GMAG())))
        DT_a     = np.asarray(DT_a)
        conds    = ["CFL-condition", "Fourier-condition", "Gravity_waves-condition"]
        cond     = [conds[i] for i in np.where(DT_a == DT_a.min())[0]]
        return COURANT * np.min(DT_a), cond


class SinglePhaseFlowTV(SPHModel):
    R""" SinglePhaseFlow solver
    """
    DENSITYMETHODS = {"SUMMATION":_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      "CONTINUITY":_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {"HARMONICAVERAGE":_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter = None,
                 solidgroup_filter = None,
                 densitymethod=None,
                 viscositymethod="HARMONICAVERAGE"):

        super().__init__(kernel, eos, nlist)

        self._param_dict.update(ParameterDict(
                          densitymethod = densitymethod,
                          viscositymethod = viscositymethod, 
                          mu = float(0.0), 
                          artificialviscosity = bool(True), 
                          alpha = float(0.2),
                          beta = float(0.0),
                          densitydiffusion = bool(False),
                          ddiff = float(0.0),
                          shepardrenormanlization = bool(False),
                          densityreinitialization = bool(False),
                          shepardfreq = int(0),
                          densityreinitfreq = int(0),
                          compute_solid_forces = bool(False),
                          max_sl = float(0.0)
                          ))




        # self._state = self._simulation.state
        self._cpp_SPFclass_name = "SinglePFTV" "_" + Kernel[self.kernel.name] + "_" + EOS[self.eos.name]
        self.fluidgroup_filter = fluidgroup_filter
        self.solidgroup_filter = solidgroup_filter
        self.str_densitymethod = densitymethod
        self.str_viscositymethod = viscositymethod
        self.accel_set = False
        self.params_set = False

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spf_cls = getattr(_sph, self._cpp_SPFclass_name)
        else:
            print("GPU not implemented")

        # check that some Particles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")
        
        # This should never happen, but leaving it in case the logic for adding
        # missed some edge case.
        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        cpp_sys_def = self._simulation.state._cpp_sys_def
        cpp_fluidgroup  = self._simulation.state._get_group(self.fluidgroup_filter)
        cpp_solidgroup  = self._simulation.state._get_group(self.solidgroup_filter)
        cpp_kernel = self.kernel.cpp_smoothingkernel
        cpp_eos = self.eos.cpp_stateequation
        cpp_nlist =  self.nlist._cpp_obj

        # Set Kernel specific Kappa in cpp-Nlist
        self.kernel.setNeighborList(self.nlist)

        self._cpp_obj = spf_cls(cpp_sys_def, cpp_kernel, cpp_eos, cpp_nlist, cpp_fluidgroup, 
                                cpp_solidgroup, self.cpp_densitymethod, self.cpp_viscositymethod)

        # Set kernel parameters
        kappa = self.kernel.Kappa()
        mycpp_kappa = self.kernel.cpp_smoothingkernel.getKernelKappa()

        pdata = self._simulation.state._cpp_sys_def.getParticleData()
        globalN = pdata.getNGlobal()

        self.consth = pdata.constSmoothingLength()
        if self.consth:
            self.maxh = pdata.getSlength(0)
            if (self._simulation.device.communicator.rank == 0):
                print(f"Using constant Smoothing Length: {self.maxh}")

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print("Non-Constant Smoothing length")
        self.rcut = kappa * self.maxh


        # Set rcut in neigbour list
        self._param_dict.update(ParameterDict(
                          rcut = self.rcut, 
                          max_sl = self.maxh
                          ))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        # get all params in line
        self.mu = self._param_dict["mu"]
        self.artificialviscosity = self._param_dict["artificialviscosity"]
        self.alpha = self._param_dict["alpha"]
        self.beta = self._param_dict["beta"]
        self.densitydiffusion = self._param_dict["densitydiffusion"]
        self.ddiff = self._param_dict["ddiff"]
        self.shepardrenormanlization = self._param_dict["shepardrenormanlization"]
        self.densityreinitialization = self._param_dict["densityreinitialization"]
        self.shepardfreq = self._param_dict["shepardfreq"]
        self.densityreinitfreq = self._param_dict["densityreinitfreq"]
        self.compute_solid_forces = self._param_dict["compute_solid_forces"]

        self.set_params(self.mu)
        self.setdensitymethod(self.str_densitymethod)
        self.setviscositymethod(self.str_viscositymethod)
        
        if (self.artificialviscosity == True):
            self.activateArtificialViscosity(self.alpha, self.beta)
        else:
            self.deactivateArtificialViscosity()
        
        if (self.densitydiffusion == True):
            self.activateDensityDiffusion(self.ddiff)
        else:
            self.deactivateDensityDiffusion()
        
        if (self.shepardrenormanlization == True):
            self.activateShepardRenormalization(self.shepardfreq)
        else:
            self.deactivateShepardRenormalization()
        
        if (self.densityreinitialization == True):
            self.activateDensityReinitialization(self.densityreinitfreq)
        else:
            self.deactivateDensityReinitialization()

        if (self.compute_solid_forces == True):
            self.computeSolidForces()

        self.setrcut(self.rcut, self.get_typelist())

        self.setBodyAcceleration(self.gx, self.gy, self.gz, self.damp)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    def _detach_hook(self):
        self.nlist._detach()

    def set_params(self,mu):
        # self.mu   = mu.item()   if isinstance(mu, np.generic)   else mu
        self._cpp_obj.setParams(self.mu)
        self.params_set = True
        self._param_dict.__setattr__("params_set", True)

    # @rcut.setter
    def setrcut(self, rcut, types):
        if rcut <= 0.0:
            raise ValueError("Rcut has to be > 0.0.")
        for p in list(combinations_with_replacement(types, 2)):
            self._cpp_obj.setRCut(p, rcut)

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    # @property
    def viscositymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.VISCOSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getViscosityMethod()]

    # @viscositymethod.setter
    def setviscositymethod(self, method):
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self._cpp_obj.setViscosityMethod(self.VISCOSITYMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.alpha   = alpha.item()  if isinstance(alpha, np.generic)   else alpha
        self.beta    = beta.item()   if isinstance(beta, np.generic)   else beta
        self._cpp_obj.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self._cpp_obj.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.ddiff   = ddiff.item()   if isinstance(ddiff, np.generic)   else ddiff
        self._cpp_obj.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self._cpp_obj.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.shepardfreq   = shepardfreq.item()   if isinstance(shepardfreq, np.generic)   else shepardfreq
        self._cpp_obj.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self._cpp_obj.deactivateShepardRenormalization()

    def activateDensityReinitialization(self, densityreinitfreq=20):
        self.densityreinitfreq   = densityreinitfreq.item()   if isinstance(densityreinitfreq, np.generic)   else densityreinitfreq
        self._cpp_obj.activateDensityReinitialization(int(densityreinitfreq))

    def deactivateDensityReinitialization(self):
        self._cpp_obj.deactivateDensityReinitialization()

    def activatePowerLaw(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw(float(K), float(n), float(mu_min))

    def activateCarreau(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian(self):
        self._cpp_obj.deactivateNonNewtonian()

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

    def setBodyAcceleration(self,gx,gy,gz,damp=0):
        self.accel_set = True
        self._param_dict.__setattr__("accel_set", True)
        # self.check_initialization();
        # self.gx   = gx.item() if isinstance(gx, np.generic) else gx
        # self.gy   = gy.item() if isinstance(gy, np.generic) else gy
        # self.gz   = gz.item() if isinstance(gz, np.generic) else gz
        # self.damp = int(damp.item()) if isinstance(damp,np.generic) else int(damp)
        self.damp = abs(self.damp)

        if ( self.gx == 0 and self.gy == 0 and self.gz == 0):
            if ( self._simulation.device.communicator.rank == 0 ):
                print(f"{self._cpp_SPFclass_name} does NOT use a body force!" )

        self._cpp_obj.setAcceleration(self.gx,self.gy,self.gz,self.damp)

    def get_speedofsound(self):
        return self.eos.SpeedOfSound

    def set_speedofsound(self, c):
        self.eos.set_speedofsound(c)

    def get_GMAG(self):
        # Magnitude of body force
        if (abs(self.gx) > 0.0 or abs(self.gy) > 0.0 or abs(self.gz) > 0.0):
            return  np.sqrt(self.gx**2+self.gy**2+self.gz**2)
        else:
            return 0.0

    def compute_speedofsound(self, LREF, UREF, DX, DRHO, H, MU, RHO0):
        # Input sanity
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")

        UREF = np.abs(UREF)

        C_a = []
        # Speed of sound
        # CFL condition
        C_a.append(UREF*UREF/DRHO)
        # Gravity waves condition
        C_a.append(self.get_GMAG()*LREF/DRHO)
        # Fourier condition
        C_a.append((MU*UREF)/(RHO0*LREF*DRHO))
        # Adami type 
        C_a.append(0.01 * self.get_GMAG() * LREF)

        # Maximum speed of sound
        C_a = np.asarray(C_a)
        conditions = ["CFL-condition", "Gravity_waves-condition", "Fourier-condition", "Adami-condition"]
        condition = [conditions[i] for i in np.where(C_a == C_a.max())[0]]
        C = np.sqrt(np.max(C_a))

        # Set speed of sound
        self.eos.set_speedofsound(C)

        return C, condition


    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")
        if H != self._param_dict["max_sl"]:
            raise ValueError("Given H not equal to stored H self._param_dict[max_sl]!")
        if MU != self._param_dict["mu"]:
            raise ValueError("Given MU not equal to stored MU self._param_dict[mu]!")
        if RHO0 != self.eos.RestDensity:
            raise ValueError("Given RHO0 not equal to stored RHO0 self.eos.RestDensity!")
        
        UREF = np.abs(UREF)

        C = self.get_speedofsound()

        DT_a = []
        # CFL condition
        # DT_1 = 0.25*H/C
        DT_a.append(DX/C)
        # Fourier condition
        DT_a.append((DX*DX*RHO0)/(8.0*MU))
        # Adami max flow
        DT_a.append(H/(C+abs(UREF)))
        # Adami viscous condition
        DT_a.append(H**2/(MU/RHO0))

        if self.get_GMAG() > 0.0:
            # Gravity waves condition
            DT_a.append(np.sqrt(H/(16.0*self.get_GMAG())))
        
        DT_a = np.asarray(DT_a)
        conditions = ["CFL-condition", "Fourier-condition", "Adami_max_flow-condition", "Adami_viscous-condition" "Gravity_waves-condition"]
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        
        DT = COURANT * np.min(DT_a)

        return DT, condition


class SinglePhaseFlowFS(SPHModel):
    R"""Free-surface SPH solver for single-phase flows (thin films, waves, jets).

    Extends the transport-velocity (TV) formulation (Adami et al. 2013) with:

    1. **Free-surface detection** via the Shepard kernel-completeness ratio
       λ_i = V_i·W₀(h) + Σ_{j≠i} V_j·W(r_ij, h).
       Particles with λ_i < ``fs_threshold`` are flagged as free-surface
       particles; their outward unit normal is n̂_i = −∇λ_i / |∇λ_i|.

    2. **Contact-angle enforcement** (Huber et al. 2016): near solid walls the
       free-surface normal is blended with the wall normal so that the liquid
       meets the wall at the prescribed equilibrium contact angle ``contact_angle``.

    3. **Curvature estimation**:
       κ_i = (1/V_i) Σ_j V_j (n̂_j − n̂_i)·∇W_ij
       Only computed for surface particles.

    4. **Free-surface pressure clamping**: P ← max(0, P) for surface particles.

    5. **CSF surface tension force**:
       F_{σ,i} = −σ · κ_i · n̂_i · (m_i/ρ_i).

    Parameters
    ----------
    kernel : SmoothingKernel
        Smoothing kernel object.
    eos : StateEquation
        Equation of state object (Tait or Linear).
    nlist : NeighborList
        Neighbour list object.
    fluidgroup_filter : hoomd.filter
        Filter selecting fluid particles.
    solidgroup_filter : hoomd.filter
        Filter selecting solid (wall) particles.
    densitymethod : str
        ``'SUMMATION'`` or ``'CONTINUITY'``.
    viscositymethod : str
        ``'HARMONICAVERAGE'`` (only option currently).
    sigma : float
        Surface tension coefficient σ [N/m].  Set to 0 to disable surface
        tension while keeping free-surface detection active.  Default 0.
    fs_threshold : float
        Kernel-completeness threshold λ ∈ (0,1).  Particles with
        λ < ``fs_threshold`` are treated as free-surface particles.
        Typical value: 0.75.  Default 0.75.
    contact_angle : float
        Equilibrium contact angle θ [rad] measured inside the liquid from
        the solid wall.  π/2 = neutral wetting (no correction),
        0 = complete wetting, π = complete non-wetting.  Default π/2.
    """

    DENSITYMETHODS = {"SUMMATION": _sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      "CONTINUITY": _sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {"HARMONICAVERAGE": _sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter=None,
                 solidgroup_filter=None,
                 densitymethod=None,
                 viscositymethod="HARMONICAVERAGE",
                 sigma=0.0,
                 fs_threshold=0.75,
                 contact_angle=np.pi / 2):

        super().__init__(kernel, eos, nlist)

        self._param_dict.update(ParameterDict(
            densitymethod=densitymethod,
            viscositymethod=viscositymethod,
            mu=float(0.0),
            artificialviscosity=bool(True),
            alpha=float(0.2),
            beta=float(0.0),
            densitydiffusion=bool(False),
            ddiff=float(0.0),
            shepardrenormanlization=bool(False),
            densityreinitialization=bool(False),
            shepardfreq=int(0),
            densityreinitfreq=int(0),
            compute_solid_forces=bool(False),
            max_sl=float(0.0),
            # FS-specific parameters
            sigma=float(sigma),
            fs_threshold=float(fs_threshold),
            contact_angle=float(contact_angle),
        ))

        self._cpp_SPFclass_name = ("SinglePFFS"
                                   "_" + Kernel[self.kernel.name]
                                   + "_" + EOS[self.eos.name])
        self.fluidgroup_filter   = fluidgroup_filter
        self.solidgroup_filter   = solidgroup_filter
        self.str_densitymethod   = densitymethod
        self.str_viscositymethod = viscositymethod
        self.accel_set           = False
        self.params_set          = False
        self.fs_params_set       = False

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            spf_cls = getattr(_sph, self._cpp_SPFclass_name)
        else:
            print("GPU not implemented")

        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")

        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        else:
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        cpp_sys_def    = self._simulation.state._cpp_sys_def
        cpp_fluidgroup = self._simulation.state._get_group(self.fluidgroup_filter)
        cpp_solidgroup = self._simulation.state._get_group(self.solidgroup_filter)
        cpp_kernel     = self.kernel.cpp_smoothingkernel
        cpp_eos        = self.eos.cpp_stateequation
        cpp_nlist      = self.nlist._cpp_obj

        self.kernel.setNeighborList(self.nlist)

        self._cpp_obj = spf_cls(cpp_sys_def, cpp_kernel, cpp_eos, cpp_nlist,
                                cpp_fluidgroup, cpp_solidgroup,
                                self.cpp_densitymethod, self.cpp_viscositymethod)

        kappa = self.kernel.Kappa()

        pdata = self._simulation.state._cpp_sys_def.getParticleData()
        self.consth = pdata.constSmoothingLength()
        if self.consth:
            self.maxh = pdata.getSlength(0)
            if self._simulation.device.communicator.rank == 0:
                print(f"Using constant Smoothing Length: {self.maxh}")
            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else:
            self.maxh = pdata.getMaxSmoothingLength()
            if self._simulation.device.communicator.rank == 0:
                print("Non-Constant Smoothing length")
        self.rcut = kappa * self.maxh

        self._param_dict.update(ParameterDict(rcut=self.rcut, max_sl=self.maxh))

        # Reload parameters (may have been set before attach)
        self.str_densitymethod   = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]

        if self.str_densitymethod == str("SUMMATION"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str("CONTINUITY"):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str("HARMONICAVERAGE"):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        # Apply all stored parameters
        self.mu                      = self._param_dict["mu"]
        self.artificialviscosity     = self._param_dict["artificialviscosity"]
        self.alpha                   = self._param_dict["alpha"]
        self.beta                    = self._param_dict["beta"]
        self.densitydiffusion        = self._param_dict["densitydiffusion"]
        self.ddiff                   = self._param_dict["ddiff"]
        self.shepardrenormanlization = self._param_dict["shepardrenormanlization"]
        self.densityreinitialization = self._param_dict["densityreinitialization"]
        self.shepardfreq             = self._param_dict["shepardfreq"]
        self.densityreinitfreq       = self._param_dict["densityreinitfreq"]
        self.compute_solid_forces    = self._param_dict["compute_solid_forces"]
        sigma_val                    = self._param_dict["sigma"]
        fs_threshold_val             = self._param_dict["fs_threshold"]
        contact_angle_val            = self._param_dict["contact_angle"]

        self.set_params(self.mu)
        self.setdensitymethod(self.str_densitymethod)
        self.setviscositymethod(self.str_viscositymethod)

        if self.artificialviscosity:
            self.activateArtificialViscosity(self.alpha, self.beta)
        else:
            self.deactivateArtificialViscosity()

        if self.densitydiffusion:
            self.activateDensityDiffusion(self.ddiff)
        else:
            self.deactivateDensityDiffusion()

        if self.shepardrenormanlization:
            self.activateShepardRenormalization(self.shepardfreq)
        else:
            self.deactivateShepardRenormalization()

        if self.densityreinitialization:
            self.activateDensityReinitialization(self.densityreinitfreq)
        else:
            self.deactivateDensityReinitialization()

        if self.compute_solid_forces:
            self.computeSolidForces()

        # Set FS-specific parameters in the C++ object
        self.setFSParams(sigma_val, fs_threshold_val, contact_angle_val)

        self.setrcut(self.rcut, self.get_typelist())
        self.setBodyAcceleration(self.gx, self.gy, self.gz, self.damp)

        super()._attach_hook()

    def _detach_hook(self):
        self.nlist._detach()

    def set_params(self, mu):
        self._cpp_obj.setParams(self.mu)
        self.params_set = True
        self._param_dict.__setattr__("params_set", True)

    def setFSParams(self, sigma, fs_threshold, contact_angle):
        """Set free-surface parameters.

        Parameters
        ----------
        sigma : float
            Surface tension coefficient σ [N/m].  Set to 0 to disable surface
            tension while keeping free-surface detection active.
        fs_threshold : float
            Kernel-completeness threshold λ ∈ (0,1).  Particles with
            λ < ``fs_threshold`` are treated as free-surface particles.
            Typical value: 0.75.
        contact_angle : float
            Equilibrium contact angle θ [rad] measured inside the liquid from
            the solid wall.  π/2 = neutral wetting (no correction applied),
            0 = complete wetting, π = complete non-wetting.
        """
        self.fs_params_set = True
        self._cpp_obj.setFSParams(float(sigma), float(fs_threshold),
                                  float(contact_angle))
        self._param_dict.__setattr__("sigma",         float(sigma))
        self._param_dict.__setattr__("fs_threshold",  float(fs_threshold))
        self._param_dict.__setattr__("contact_angle", float(contact_angle))

    # ── Methods delegating to C++ object ─────────────────────────────────────

    def setrcut(self, rcut, types):
        if rcut <= 0.0:
            raise ValueError("Rcut has to be > 0.0.")
        for p in list(combinations_with_replacement(types, 2)):
            self._cpp_obj.setRCut(p, rcut)

    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    def setviscositymethod(self, method):
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self._cpp_obj.setViscosityMethod(self.VISCOSITYMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.alpha = float(alpha)
        self.beta  = float(beta)
        self._cpp_obj.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self._cpp_obj.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.ddiff = float(ddiff)
        self._cpp_obj.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self._cpp_obj.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.shepardfreq = int(shepardfreq)
        self._cpp_obj.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self._cpp_obj.deactivateShepardRenormalization()

    def activateDensityReinitialization(self, densityreinitfreq=20):
        self.densityreinitfreq = int(densityreinitfreq)
        self._cpp_obj.activateDensityReinitialization(int(densityreinitfreq))

    def deactivateDensityReinitialization(self):
        self._cpp_obj.deactivateDensityReinitialization()

    def activatePowerLaw(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw(float(K), float(n), float(mu_min))

    def activateCarreau(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian(self):
        self._cpp_obj.deactivateNonNewtonian()

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

    def setBodyAcceleration(self, gx, gy, gz, damp=0):
        self.accel_set = True
        self._param_dict.__setattr__("accel_set", True)
        self.damp = abs(self.damp)
        if self.gx == 0 and self.gy == 0 and self.gz == 0:
            if self._simulation.device.communicator.rank == 0:
                print(f"{self._cpp_SPFclass_name} does NOT use a body force!")
        self._cpp_obj.setAcceleration(self.gx, self.gy, self.gz, self.damp)

    def get_speedofsound(self):
        return self.eos.SpeedOfSound

    def set_speedofsound(self, c):
        self.eos.set_speedofsound(c)

    def get_GMAG(self):
        if abs(self.gx) > 0.0 or abs(self.gy) > 0.0 or abs(self.gz) > 0.0:
            return np.sqrt(self.gx**2 + self.gy**2 + self.gz**2)
        return 0.0

    def compute_speedofsound(self, LREF, UREF, DX, DRHO, H, MU, RHO0):
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")
        UREF  = np.abs(UREF)
        C_a   = []
        C_a.append(UREF * UREF / DRHO)
        C_a.append(self.get_GMAG() * LREF / DRHO)
        C_a.append((MU * UREF) / (RHO0 * LREF * DRHO))
        C_a.append(0.01 * self.get_GMAG() * LREF)
        C_a = np.asarray(C_a)
        conditions = ["CFL-condition", "Gravity_waves-condition",
                      "Fourier-condition", "Adami-condition"]
        condition = [conditions[i] for i in np.where(C_a == C_a.max())[0]]
        C = np.sqrt(np.max(C_a))
        self.eos.set_speedofsound(C)
        return C, condition

    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.25):
        if LREF == 0.0:
            raise ValueError("Reference length LREF may not be zero.")
        if DRHO == 0.0:
            raise ValueError("Maximum density variation DRHO may not be zero.")
        if DX <= 0.0:
            raise ValueError("DX may not be zero or negative.")
        UREF  = np.abs(UREF)
        C     = self.get_speedofsound()
        DT_a  = []
        DT_a.append(DX / C)
        DT_a.append((DX * DX * RHO0) / (8.0 * MU))
        DT_a.append(H / (C + abs(UREF)))
        DT_a.append(H**2 / (MU / RHO0))
        if self.get_GMAG() > 0.0:
            DT_a.append(np.sqrt(H / (16.0 * self.get_GMAG())))
        DT_a      = np.asarray(DT_a)
        conditions = ["CFL-condition", "Fourier-condition",
                      "Adami_max_flow-condition", "Adami_viscous-condition",
                      "Gravity_waves-condition"]
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        return COURANT * np.min(DT_a), condition


class TwoPhaseFlow(SPHModel):
    R""" TwoPhaseFlow solver
    """
    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    COLORGRADIENTMETHODS = {'DENSITYRATIO':_sph.PhaseFlowColorGradientMethod.DENSITYRATIO,
                            'NUMBERDENSITY':_sph.PhaseFlowColorGradientMethod.NUMBERDENSITY}

    def __init__(self,
                 kernel,
                 eos1,
                 eos2,
                 nlist,
                 fluidgroup1_filter = None,
                 fluidgroup2_filter = None,
                 solidgroup_filter = None,
                 densitymethod = None,
                 viscositymethod = 'HARMONICAVERAGE',
                 colorgradientmethod = 'DENSITYRATIO'):

        super().__init__(kernel, eos1, nlist)

        self._param_dict.update(ParameterDict(
                          densitymethod = densitymethod,
                          viscositymethod = viscositymethod,
                          colorgradientmethod = colorgradientmethod, 
                          mu1 = float(0.0), 
                          mu2 = float(0.0),
                          sigma12 = float(0.0), 
                          omega = float(0.0),
                          omega_adv = float(180.0),
                          omega_rec = float(0.0),
                          hysteresis = bool(False),
                          artificialviscosity = bool(True),
                          alpha = float(0.2),
                          beta = float(0.0),
                          densitydiffusion = bool(False),
                          ddiff = float(0.0),
                          shepardrenormanlization = bool(False),
                          shepardfreq = int(30),
                          compute_solid_forces = bool(False),
                          fickian_shifting = bool(False),
                          max_sl = float(0.0)
                          ))




        # self._state = self._simulation.state
        self.eos1 = eos1
        self.eos2 = eos2
        self._cpp_TPFclass_name = 'TwoPF' '_' + Kernel[self.kernel.name] + '_' + EOS[self.eos1.name] + EOS[self.eos2.name]
        self.fluidgroup1_filter = fluidgroup1_filter
        self.fluidgroup2_filter = fluidgroup2_filter
        self.solidgroup_filter = solidgroup_filter
        self.str_densitymethod = densitymethod
        self.str_viscositymethod = viscositymethod
        self.str_colorgradientmethod = colorgradientmethod
        self.accel_set = False
        self.params_set = False

        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str('HARMONICAVERAGE'):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        if self.str_colorgradientmethod == str('DENSITYRATIO'):
            self.cpp_colorgradientmethod = hoomd.sph._sph.PhaseFlowColorGradientMethod.DENSITYRATIO
        elif self.str_colorgradientmethod == str('NUMBERDENSITY'):
            self.cpp_colorgradientmethod = hoomd.sph._sph.PhaseFlowColorGradientMethod.NUMBERDENSITY
        else:
            raise ValueError("Using undefined ColorGradientMethod.")

    def _attach_hook(self):
        if isinstance(self._simulation.device, hoomd.device.CPU):
            tpf_cls = getattr(_sph, self._cpp_TPFclass_name)
        else:
            print("GPU not implemented")

        # check that some Particles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")
        
        # This should never happen, but leaving it in case the logic for adding
        # missed some edge case.
        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happening since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)

        cpp_sys_def = self._simulation.state._cpp_sys_def
        cpp_fluidgroup1  = self._simulation.state._get_group(self.fluidgroup1_filter)
        cpp_fluidgroup2  = self._simulation.state._get_group(self.fluidgroup2_filter)
        cpp_solidgroup  = self._simulation.state._get_group(self.solidgroup_filter)
        cpp_kernel = self.kernel.cpp_smoothingkernel
        cpp_eos1 = self.eos1.cpp_stateequation
        cpp_eos2 = self.eos2.cpp_stateequation
        cpp_nlist =  self.nlist._cpp_obj

        # Set Kernel specific Kappa in cpp-Nlist
        self.kernel.setNeighborList(self.nlist)

        self._cpp_obj = tpf_cls(cpp_sys_def, cpp_kernel, 
                                cpp_eos1, cpp_eos2, 
                                cpp_nlist, 
                                cpp_fluidgroup1, cpp_fluidgroup2, 
                                cpp_solidgroup, 
                                self.cpp_densitymethod, 
                                self.cpp_viscositymethod, 
                                self.cpp_colorgradientmethod)

        # Set kernel parameters
        kappa = self.kernel.Kappa()
        mycpp_kappa = self.kernel.cpp_smoothingkernel.getKernelKappa()

        pdata = self._simulation.state._cpp_sys_def.getParticleData()
        globalN = pdata.getNGlobal()

        self.consth = pdata.constSmoothingLength()
        if self.consth:
            self.maxh = pdata.getSlength(0)
            if (self._simulation.device.communicator.rank == 0):
                print(f'Using constant Smoothing Length: {self.maxh}')

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print('Non-Constant Smoothing length')
        self.rcut = kappa * self.maxh


        # Set rcut in neigbour list
        self._param_dict.update(ParameterDict(
                          rcut = self.rcut, 
                          max_sl = self.maxh
                          ))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]
        self.str_colorgradientmethod = self._param_dict._dict["colorgradientmethod"]

        if self.str_densitymethod == str('SUMMATION'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYSUMMATION
        elif self.str_densitymethod == str('CONTINUITY'):
            self.cpp_densitymethod = hoomd.sph._sph.PhaseFlowDensityMethod.DENSITYCONTINUITY
        else:
            raise ValueError("Using undefined DensityMethod.")

        if self.str_viscositymethod == str('HARMONICAVERAGE'):
            self.cpp_viscositymethod = hoomd.sph._sph.PhaseFlowViscosityMethod.HARMONICAVERAGE
        else:
            raise ValueError("Using undefined ViscosityMethod.")

        if self.str_colorgradientmethod == str('DENSITYRATIO'):
            self.cpp_colorgradientmethod = hoomd.sph._sph.PhaseFlowColorGradientMethod.DENSITYRATIO
        elif self.str_colorgradientmethod == str('NUMBERDENSITY'):
            self.cpp_colorgradientmethod = hoomd.sph._sph.PhaseFlowColorGradientMethod.NUMBERDENSITY
        else:
            raise ValueError("Using undefined ColorGradientMethod.")

        # get all params in line
        self.mu1 = self._param_dict['mu1']
        self.mu2 = self._param_dict['mu2']
        self.sigma12 = self._param_dict['sigma12']
        self.omega = self._param_dict['omega']
        self.artificialviscosity = self._param_dict['artificialviscosity']
        self.alpha = self._param_dict['alpha']
        self.beta = self._param_dict['beta']
        self.densitydiffusion = self._param_dict['densitydiffusion']
        self.ddiff = self._param_dict['ddiff']
        self.shepardrenormanlization = self._param_dict['shepardrenormanlization']
        self.shepardfreq = self._param_dict['shepardfreq']
        self.compute_solid_forces = self._param_dict['compute_solid_forces']
        self.fickian_shifting = self._param_dict['fickian_shifting']

        self.set_params(self.mu1, self.mu2, self.sigma12, self.omega)
        self.hysteresis = self._param_dict['hysteresis']
        self.omega_adv  = self._param_dict['omega_adv']
        self.omega_rec  = self._param_dict['omega_rec']
        if self.hysteresis:
            self._cpp_obj.setHysteresis(self.omega_rec, self.omega_adv)
        self.setdensitymethod(self.str_densitymethod)
        self.setviscositymethod(self.str_viscositymethod)
        self.setcolorgradientmethod(self.str_colorgradientmethod)
        
        if (self.artificialviscosity == True):
            self.activateArtificialViscosity(self.alpha, self.beta)
        else:
            self.deactivateArtificialViscosity()
        
        if (self.densitydiffusion == True):
            self.activateDensityDiffusion(self.ddiff)
        else:
            self.deactivateDensityDiffusion()
        
        if (self.shepardrenormanlization == True):
            self.activateShepardRenormalization(self.shepardfreq)
        else:
            self.deactivateShepardRenormalization()
        
        if (self.compute_solid_forces == True):
            self.computeSolidForces()

        if (self.fickian_shifting == True):
            self.activateFickianShifting()

        self.setrcut(self.rcut, self.get_typelist())

        self.setBodyAcceleration(self.gx, self.gy, self.gz, self.damp)

        # Attach param_dict and typeparam_dict
        super()._attach_hook()

    def _detach_hook(self):
        self.nlist._detach()

    def set_params(self, mu1, mu2, sigma12, omega):
        self._cpp_obj.setParams(self.mu1, self.mu2, self.sigma12, self.omega)
        self.params_set = True
        self._param_dict.__setattr__('params_set', True)

    # @rcut.setter
    def setrcut(self, rcut, types):
        if rcut <= 0.0:
            raise ValueError("Rcut has to be > 0.0.")
        for p in combinations_with_replacement(types, 2):
            self._cpp_obj.setRCut(p, rcut)

    # @property
    def densitymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.DENSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getDensityMethod()]

    # @densitymethod.setter
    def setdensitymethod(self, method):
        if method not in self.DENSITYMETHODS:
            raise ValueError("Undefined DensityMethod.")
        self._cpp_obj.setDensityMethod(self.DENSITYMETHODS[method])

    # @property
    def viscositymethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.VISCOSITYMETHODS.iteritems())
        return invD[self._cpp_obj.getViscosityMethod()]

    # @viscositymethod.setter
    def setviscositymethod(self, method):
        if method not in self.VISCOSITYMETHODS:
            raise ValueError("Undefined ViscosityMethod.")
        self._cpp_obj.setViscosityMethod(self.VISCOSITYMETHODS[method])

    # @property
    def colorgradientmethod(self):
        # Invert key mapping
        invD = dict((v,k) for k, v in self.COLORGRADIENTMETHODS.iteritems())
        return invD[self._cpp_obj.getColorGradientMethod()]

    # @colorgradientmethod.setter
    def setcolorgradientmethod(self, method):
        if method not in self.COLORGRADIENTMETHODS:
            raise ValueError("Undefined ColorGradientMethod.")
        self._cpp_obj.setColorGradientMethod(self.COLORGRADIENTMETHODS[method])

    def activateArtificialViscosity(self, alpha, beta):
        self.alpha   = alpha.item()  if isinstance(alpha, np.generic)   else alpha
        self.beta    = beta.item()   if isinstance(beta, np.generic)   else beta
        self._cpp_obj.activateArtificialViscosity(alpha, beta)

    def deactivateArtificialViscosity(self):
        self._cpp_obj.deactivateArtificialViscosity()

    def activateDensityDiffusion(self, ddiff):
        self.ddiff   = ddiff.item()   if isinstance(ddiff, np.generic)   else ddiff
        self._cpp_obj.activateDensityDiffusion(ddiff)

    def deactivateDensityDiffusion(self):
        self._cpp_obj.deactivateDensityDiffusion()

    def activateShepardRenormalization(self, shepardfreq=30):
        self.shepardfreq   = shepardfreq.item()   if isinstance(shepardfreq, np.generic)   else shepardfreq
        self._cpp_obj.activateShepardRenormalization(int(shepardfreq))

    def deactivateShepardRenormalization(self):
        self._cpp_obj.deactivateShepardRenormalization()

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

    def activatePowerLaw1(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw1(float(K), float(n), float(mu_min))

    def activateCarreau1(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau1(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham1(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham1(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley1(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley1(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian1(self):
        self._cpp_obj.deactivateNonNewtonian1()

    def activatePowerLaw2(self, K, n, mu_min=0.0):
        self._cpp_obj.activatePowerLaw2(float(K), float(n), float(mu_min))

    def activateCarreau2(self, mu0, mu_inf, lam, n):
        self._cpp_obj.activateCarreau2(float(mu0), float(mu_inf), float(lam), float(n))

    def activateBingham2(self, mu_p, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateBingham2(float(mu_p), float(tau_y), float(m_reg), float(mu_min))

    def activateHerschelBulkley2(self, K, n, tau_y, m_reg, mu_min=0.0):
        self._cpp_obj.activateHerschelBulkley2(float(K), float(n), float(tau_y), float(m_reg), float(mu_min))

    def deactivateNonNewtonian2(self):
        self._cpp_obj.deactivateNonNewtonian2()

    def activateFickianShifting(self):
        self._cpp_obj.activateFickianShifting()

    def setBodyAcceleration(self,gx,gy,gz,damp=0):
        self.accel_set = True
        self._param_dict.__setattr__('accel_set', True)
        # self.check_initialization();
        # self.gx   = gx.item() if isinstance(gx, np.generic) else gx
        # self.gy   = gy.item() if isinstance(gy, np.generic) else gy
        # self.gz   = gz.item() if isinstance(gz, np.generic) else gz
        # self.damp = int(damp.item()) if isinstance(damp,np.generic) else int(damp)
        self.damp = abs(self.damp)

        if ( self.gx == 0 and self.gy == 0 and self.gz == 0):
            if ( self._simulation.device.communicator.rank == 0 ):
                print(f'{self._cpp_TPFclass_name} does NOT use a body force!' )

        self._cpp_obj.setAcceleration(self.gx,self.gy,self.gz,self.damp)

    def get_speedofsound(self):
        return self.eos1.SpeedOfSound, self.eos2.SpeedOfSound

    def set_speedofsound(self, c1, c2):
        self.eos1.set_speedofsound(c1)
        self.eos2.set_speedofsound(c2)

    def get_GMAG(self):
        # Magnitude of body force
        if (abs(self.gx) > 0.0 or abs(self.gy) > 0.0 or abs(self.gz) > 0.0):
            return  np.sqrt(self.gx**2+self.gy**2+self.gz**2)
        else:
            return 0.0

    def compute_speedofsound(self, LREF, UREF, DX, DRHO, H, MU1, MU2, RHO01, RHO02, SIGMA12):
        # Input sanity
        if LREF == 0.0:
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')

        UREF = np.abs(UREF)

        C_a1 = []
        C_a2 = []
        # Speed of sound
        # CFL condition
        C_a1.append(UREF*UREF/DRHO)
        C_a2.append(UREF*UREF/DRHO)
        # Gravity waves condition
        C_a1.append(self.get_GMAG()*LREF/DRHO)
        C_a2.append(self.get_GMAG()*LREF/DRHO)
        # Surface Wave condition
        C_a1.append( SIGMA12/(RHO01 * LREF * DRHO) )
        C_a2.append( SIGMA12/(RHO02 * LREF * DRHO) )
        # Fourier condition
        C_a1.append((MU1*UREF)/(RHO01*LREF*DRHO))
        C_a2.append((MU2*UREF)/(RHO02*LREF*DRHO))
        # Maximum speed of sound
        
        C_a1 = np.asarray(C_a1)
        C_a2 = np.asarray(C_a2)
        conditions = ['CFL-condition', 'Gravity_waves-condition', 'Surface_waves-condition', 'Fourier-condition']
        condition1 = [conditions[i] for i in np.where(C_a1 == C_a1.max())[0]]
        condition2 = [conditions[i] for i in np.where(C_a2 == C_a2.max())[0]]
        C1 = np.sqrt(np.max(C_a1))
        C2 = np.sqrt(np.max(C_a2))

        # Set speed of sound
        self.eos1.set_speedofsound(C1)
        self.eos2.set_speedofsound(C2)

        return C1, condition1, C2, condition2


    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU1, MU2, RHO01, RHO02, SIGMA12, COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')
        if H != self._param_dict['max_sl']:
            raise ValueError('Given H not equal to stored H self._param_dict[max_sl]!')
        if MU1 != self._param_dict['mu1'] or MU2 != self._param_dict['mu2'] :
            raise ValueError('Given MU not equal to stored MU self._param_dict[mu]!')
        if RHO01 != self.eos1.RestDensity or RHO02 != self.eos2.RestDensity:
            raise ValueError('Given RHO0 not equal to stored RHO0 self.eos.RestDensity!')
        
        UREF = np.abs(UREF)

        C1, C2 = self.get_speedofsound()

        DT_a = []
        # CFL condition
        # DT_1 = 0.25*H/C
        DT_a.append( DX/np.max( [C1, C2] ) )
        # Fourier condition
        DT_a.append((DX*DX*np.max([RHO01, RHO02]))/(8.0 * np.max([MU1, MU2]) ))
        # Surface Waves condition
        if SIGMA12 > 0.0:
            DT_a.append( np.sqrt((DX * DX * DX * np.min([RHO01, RHO02]))/(32.0 * np.pi * SIGMA12)) )

        if self.get_GMAG() > 0.0:
            # Gravity waves condition
            DT_a.append(np.sqrt(DX/(16.0*self.get_GMAG())))
        DT_a = np.asarray(DT_a)
        conditions = ['CFL-condition', 'Fourier-condition', 'Surface_waves-condition', 'Gravity_waves-condition']
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        DT = COURANT * np.min(DT_a)

        return DT, condition




class TwoPhaseFlowTV(TwoPhaseFlow):
    R"""Two-phase SPH with transport-velocity (TV) formulation (Adami et al. 2013).

    Extends TwoPhaseFlow with two additional momentum terms:

    1. Artificial-stress correction (Adami 2013, Eq. 11):
         F_i^AS += (1/2) Σ_j (Vi²+Vj²) (A_i+A_j)·∇W_ij
         A_k = ρ_k v_k ⊗ (tv_k − v_k)
       Counteracts the tensile instability (particle clustering) by penalising
       deviations of the transport velocity from the physical velocity.

    2. Background-pressure contribution (BPC, written to aux2):
         bpc_i = −Σ_j (Vi²+Vj²) P_bg/m_i · (∂W/∂r)/r · r_ij
       P_bg = eos.TransportVelocityPressure (set via eos.set_params(..., tvp=P_bg)).
       KickDriftKickTV reads bpc_i on the next half-step to advance particle
       positions along the smooth transport velocity field.

    All TwoPhaseFlow options remain available (CIP, Riemann dissipation,
    density diffusion, surface tension, Fickian shifting).

    IMPORTANT: Use KickDriftKickTV as the integration method — NOT VelocityVerletBasic.
    The EOS transport velocity pressure must be set:
        eos1.set_params(rho01, backpressure, tvp=P_bg1)
        eos2.set_params(rho02, backpressure, tvp=P_bg2)
    Typical choice: P_bg = backpressure_coefficient * rho0 * c^2.

    References:
        Adami, Hu & Adams (2013) J. Comput. Phys. 241, 292–307
    """

    def __init__(self,
                 kernel,
                 eos1,
                 eos2,
                 nlist,
                 fluidgroup1_filter=None,
                 fluidgroup2_filter=None,
                 solidgroup_filter=None,
                 densitymethod=None,
                 viscositymethod='HARMONICAVERAGE',
                 colorgradientmethod='DENSITYRATIO'):

        super().__init__(kernel, eos1, eos2, nlist,
                         fluidgroup1_filter, fluidgroup2_filter, solidgroup_filter,
                         densitymethod, viscositymethod, colorgradientmethod)

        # Override C++ class name — uses TwoPFTV_* instead of TwoPF_*
        self._cpp_TPFclass_name = ('TwoPFTV' '_'
                                   + Kernel[self.kernel.name] + '_'
                                   + EOS[self.eos1.name] + EOS[self.eos2.name])

        # Extra parameters for TV-specific optional features
        self._param_dict.update(ParameterDict(
            riemann_dissipation=bool(False),
            riemann_beta=float(1.0),
            consistent_interface_pressure=bool(False),
        ))

        self.riemann_dissipation = False
        self.riemann_beta = 1.0
        self.consistent_interface_pressure = False

    def _attach_hook(self):
        # Parent creates self._cpp_obj using self._cpp_TPFclass_name (TwoPFTV_*)
        super()._attach_hook()

        # Apply TV-specific settings after C++ object is created
        self.riemann_dissipation = self._param_dict['riemann_dissipation']
        self.riemann_beta = self._param_dict['riemann_beta']
        self.consistent_interface_pressure = self._param_dict['consistent_interface_pressure']

        if self.riemann_dissipation:
            self._cpp_obj.activateRiemannDissipation(float(self.riemann_beta))
        else:
            self._cpp_obj.deactivateRiemannDissipation()

        if self.consistent_interface_pressure:
            self._cpp_obj.activateConsistentInterfacePressure()
        else:
            self._cpp_obj.deactivateConsistentInterfacePressure()

    def activateRiemannDissipation(self, beta=1.0):
        """Activate Riemann-based dissipation (Zhang, Hu & Adams 2017).

        Replaces Monaghan AV with impedance-weighted upwind dissipation.
        Mutually exclusive with activateArtificialViscosity.

        Args:
            beta: Scaling coefficient β_R (default 1.0).
        """
        self.riemann_beta = float(beta)
        self._cpp_obj.activateRiemannDissipation(self.riemann_beta)

    def deactivateRiemannDissipation(self):
        """Deactivate Riemann-based dissipation."""
        self._cpp_obj.deactivateRiemannDissipation()

    def activateConsistentInterfacePressure(self):
        """Activate consistent interface pressure (Hu & Adams 2009).

        Uses rest-density weighting plus hydrostatic correction for
        cross-phase pressure averages.  Reduces parasitic currents at
        large density-ratio interfaces.  Requires SUMMATION density method.
        """
        self._cpp_obj.activateConsistentInterfacePressure()

    def deactivateConsistentInterfacePressure(self):
        """Deactivate consistent interface pressure."""
        self._cpp_obj.deactivateConsistentInterfacePressure()


# Dicts
Kernel = {"_WendlandC2":"WC2","_WendlandC4":"WC4","_WendlandC6":"WC6","_Quintic":"Q","_CubicSpline":"CS"}
EOS = {"_Linear":"L","_Tait":"T"}


