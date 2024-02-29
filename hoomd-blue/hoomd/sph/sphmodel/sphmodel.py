f"""----------------------------------------------------------
maintainer: dkrach, david.krach@mib.uni-stuttgart.de
-----------------------------------------------------------"""

"""SPH Momentum interaction forces."""

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
        # params["exclusions"] = exclusions
        self._param_dict.update(params)

        self.kernel     = kernel
        self.eos        = eos

        # # Detect maximum cut-off radius for nlist parameterization
        # kappa = self.kernel.Kappa()
        # print(dir(self))
        # pdata = self.state._cpp_sys_def.getParticleData()
        # globalN   = pdata.getNGlobal()

        # maxh      = pdata.getMaxSmoothingLength()
        # self.rcut = kappa * maxh

        # # Neighbor list
        # self.nlist.subscribe(lambda:self.get_rcut())
        # self.nlist.update_rcut()
        # # self.nlist.cpp_nlist.setStorageMode(_nsearch.NeighborList.storageMode.full);

        # # Set neighbor list in kernel class
        # self.kernel.setNeighborList(self.nlist)

        # print(self._simulation.device)

        type_params = []
        self._extend_typeparam(type_params)
        self._param_dict.update(
            ParameterDict(nlist=hoomd.nsearch.nlist.NeighborList))
        self.nlist = nlist


    # def _add(self, simulation):
    #     super()._add(simulation)
    #     self._add_nlist()

    # def _add_nlist(self):
    #     nlist = self.nlist
    #     deepcopy = False
    #     if not isinstance(self._simulation, hoomd.Simulation):
    #         if nlist._added:
    #             deepcopy = True
    #         else:
    #             return
    #     # We need to check if the force is added since if it is not then this is
    #     # being called by a SyncedList object and a disagreement between the
    #     # simulation and nlist._simulation is an error. If the force is added
    #     # then the nlist is compatible. We cannot just check the nlist's _added
    #     # property because _add is also called when the SyncedList is synced.
    #     if deepcopy or nlist._added and nlist._simulation != self._simulation:
    #         warnings.warn(
    #             f"{self} object is creating a new equivalent neighbor list."
    #             f" This is happending since the force is moving to a new "
    #             f"simulation. Set a new nlist to suppress this warning.",
    #             RuntimeWarning)
    #         self.nlist = copy.deepcopy(nlist)
    #     self.nlist._add(self._simulation)
    #     # This is ideopotent, but we need to ensure that if we change
    #     # neighbor list when not attached we handle correctly.
    #     self._add_dependency(self.nlist)

    # def _attach_hook(self):
    #     # check that some Particles are defined
    #     if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
    #         self._simulation.device._cpp_msg.warning("No particles are defined.\n")
        
    #     # This should never happen, but leaving it in case the logic for adding
    #     # missed some edge case.
    #     if self.nlist._attached and self._simulation != self.nlist._simulation:
    #         raise RuntimeError("{} object's neighbor list is used in a "
    #                            "different simulation.".format(type(self)))
    #     self.nlist._attach()
    #     if isinstance(self._simulation.device, hoomd.device.CPU):
    #         self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
    #         # self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.half)
    #     else:
    #         self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
    #     # self._cpp_obj = cls(self._simulation.state._cpp_sys_def,
    #     #                      self.nlist._cpp_obj)

    #     self._cpp_baseclass_name = 'SPHBaseClass' + '_' + Kernel[self.kernel.name] + '_' + EOS[self.eos.name]
    #     base_cls = getattr(_sph, self._cpp_baseclass_name)
    #     self._cpp_base_obj = base_cls(self._simulation.state._cpp_sys_def, self.kernel.cpp_smoothingkernel,
    #                              self.eos.cpp_stateequation, self.nlist._cpp_obj)
    #     super()._attach_hook()

    def _attach_hook(self):
        # check that some Particles are defined
        if self._simulation.state._cpp_sys_def.getParticleData().getNGlobal() == 0:
            self._simulation.device._cpp_msg.warning("No particles are defined.\n")

        if self.nlist._attached and self._simulation != self.nlist._simulation:
            warnings.warn(
                f"{self} object is creating a new equivalent neighbor list."
                f" This is happending since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            # cls = getattr(self._ext_module_nlist, self._cpp_class_name)
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        else:
            # cls = getattr(self._ext_module_nlist, self._cpp_class_name + "GPU")
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        # self._cpp_obj = cls(self._simulation.state._cpp_sys_def,
        #                     self.nlist._cpp_obj)

        self._cpp_baseclass_name = 'SPHBaseClass' + '_' + Kernel[self.kernel.name] + '_' + EOS[self.eos.name]
        base_cls = getattr(_sph, self._cpp_baseclass_name)
        self._cpp_base_obj = base_cls(self._simulation.state._cpp_sys_def, self.kernel.cpp_smoothingkernel,
                                 self.eos.cpp_stateequation, self.nlist._cpp_obj)



    # def _detach_hook(self):
    #     self.nlist._detach()

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
        # old_nlist = self.nlist
        self._param_dict._dict["nlist"] = new_nlist
        # if self._added:
        #     self._add_nlist()
        #     old_nlist._remove_dependent(self)

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
    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter = None,
                 solidgroup_filter = None,
                 densitymethod=None,
                 viscositymethod='HARMONICAVERAGE'):

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
                          shepardfreq = int(0),
                          compute_solid_forces = bool(False),
                          max_sl = float(0.0)
                          ))




        # self._state = self._simulation.state
        self._cpp_SPFclass_name = 'SinglePF' '_' + Kernel[self.kernel.name] + '_' + EOS[self.eos.name]
        self.fluidgroup_filter = fluidgroup_filter
        self.solidgroup_filter = solidgroup_filter
        self.str_densitymethod = densitymethod
        self.str_viscositymethod = viscositymethod
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


        # IMPORTANT TO ADD CUDA! 
        # create the c++ mirror class
        # if isinstance(self._simulation.device, hoomd.device.CPU):
        #     _cpp_class_name += 'GPU'


        # print(self._cpp_class_name)
        # cls = getattr(_sph, self._cpp_class_name)
        # self._cpp_obj = ()

        # self._attach()
        # Check if group is associated
        # if fluidgroup == None or solidgroup == None:
        #     raise ValueError("SinglePhaseFlow: fluidgroup and solidgroup can not be None")
        # else:
        #     self.fluidgroup = fluidgroup
        #     self.solidgroup = solidgroup
        #     self.cpp_fluidgroup = fluidgroup.cpp_group
        #     self.cpp_solidgroup = solidgroup.cpp_group

        # # Check if given smoothing lengths are equal, if so
        #     # set constant smoothing length flag in model class
        # pdata = self._simulation.state._cpp_sys_def.getParticleData()
        # globalN   = pdata.getNGlobal()

        # self.consth = pdata.constSmoothingLength()

        # if self.consth:
        #     self.maxh = pdata.getSmoothingLength(0)
        #     self._cpp_obj.setConstSmoothingLength(self.maxh)
        # else:
        #     self.maxh = pdata.getMaxSmoothingLength()

    # def _add(self, simulation):
    #     super()._add(simulation)

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
                f" This is happending since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            # self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.half)
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        # TODO: understand why _nsearch.NeighborList.storageMode.half makes wierd errors!

        else:
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
                print(f'Using constant Smooting Length: {self.maxh}')

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print('Non-Constant Smooting length')
        self.rcut = kappa * self.maxh

        # print(f'self.rcut {self.rcut}')

        # Set rcut in neigbour list
        self._param_dict.update(ParameterDict(
                          rcut = self.rcut, 
                          max_sl = self.maxh
                          ))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]

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

        # get all params in line
        self.mu = self._param_dict['mu']
        self.artificialviscosity = self._param_dict['artificialviscosity']
        self.alpha = self._param_dict['alpha']
        self.beta = self._param_dict['beta']
        self.densitydiffusion = self._param_dict['densitydiffusion']
        self.ddiff = self._param_dict['ddiff']
        self.shepardrenormanlization = self._param_dict['shepardrenormanlization']
        self.shepardfreq = self._param_dict['shepardfreq']
        print(self.shepardfreq)
        print(self.mu)
        self.compute_solid_forces = self._param_dict['compute_solid_forces']

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
                print(f'{self._cpp_SPFclass_name} does NOT use a body force!' )

        # self.cpp_force.setAcceleration(self.gx,self.gy,self.gz,self.damp)
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
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')

        UREF = np.abs(UREF)

        C_a = []
        # Speed of sound
        # CFL condition
        C_a.append(UREF*UREF/DRHO)
        # Gravity waves condition
        C_a.append(self.get_GMAG()*LREF/DRHO)
        # Fourier condition
        C_a.append((MU*UREF)/(RHO0*LREF*DRHO))
        # Maximum speed of sound
        
        C_a = np.asarray(C_a)
        conditions = ['CFL-condition', 'Gravity_waves-condition', 'Fourier-condition']
        condition = [conditions[i] for i in np.where(C_a == C_a.max())[0]]
        C = np.sqrt(np.max(C_a))

        # Set speed of sound
        self.eos.set_speedofsound(C)

        return C, condition


    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')
        if H != self._param_dict['max_sl']:
            raise ValueError('Given H not equal to stored H self._param_dict[max_sl]!')
        if MU != self._param_dict['mu']:
            raise ValueError('Given MU not equal to stored MU self._param_dict[mu]!')
        if RHO0 != self.eos.RestDensity:
            raise ValueError('Given RHO0 not equal to stored RHO0 self.eos.RestDensity!')
        
        UREF = np.abs(UREF)

        C = self.get_speedofsound()

        DT_a = []
        # CFL condition
        # DT_1 = 0.25*H/C
        DT_a.append(DX/C)
        # Fourier condition
        DT_a.append((DX*DX*RHO0)/(8.0*MU))
        
        if self.get_GMAG() > 0.0:
            # Gravity waves condition
            DT_a.append(np.sqrt(H/(16.0*self.get_GMAG())))
        DT_a = np.asarray(DT_a)
        conditions = ['CFL-condition', 'Fourier-condition', 'Gravity_waves-condition']
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        DT = COURANT * np.min(DT_a)

        return DT, condition


class SinglePhaseFlowTV(SPHModel):
    R""" SinglePhaseFlow solver
    """
    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    def __init__(self,
                 kernel,
                 eos,
                 nlist,
                 fluidgroup_filter = None,
                 solidgroup_filter = None,
                 densitymethod=None,
                 viscositymethod='HARMONICAVERAGE'):

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
                          shepardfreq = int(0),
                          compute_solid_forces = bool(False),
                          max_sl = float(0.0)
                          ))




        # self._state = self._simulation.state
        self._cpp_SPFclass_name = 'SinglePFTV' '_' + Kernel[self.kernel.name] + '_' + EOS[self.eos.name]
        self.fluidgroup_filter = fluidgroup_filter
        self.solidgroup_filter = solidgroup_filter
        self.str_densitymethod = densitymethod
        self.str_viscositymethod = viscositymethod
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
                f" This is happending since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            # self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.half)
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        # TODO: understand why _nsearch.NeighborList.storageMode.half makes wierd errors!

        else:
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
                print(f'Using constant Smooting Length: {self.maxh}')

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print('Non-Constant Smooting length')
        self.rcut = kappa * self.maxh

        # print(f'self.rcut {self.rcut}')

        # Set rcut in neigbour list
        self._param_dict.update(ParameterDict(
                          rcut = self.rcut, 
                          max_sl = self.maxh
                          ))

        # Reload density and viscosity methods from __dict__
        self.str_densitymethod = self._param_dict._dict["densitymethod"]
        self.str_viscositymethod = self._param_dict._dict["viscositymethod"]

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

        # get all params in line
        self.mu = self._param_dict['mu']
        self.artificialviscosity = self._param_dict['artificialviscosity']
        self.alpha = self._param_dict['alpha']
        self.beta = self._param_dict['beta']
        self.densitydiffusion = self._param_dict['densitydiffusion']
        self.ddiff = self._param_dict['ddiff']
        self.shepardrenormanlization = self._param_dict['shepardrenormanlization']
        self.shepardfreq = self._param_dict['shepardfreq']
        self.compute_solid_forces = self._param_dict['compute_solid_forces']

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
        self._param_dict.__setattr__('params_set', True)

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

    def computeSolidForces(self):
        self._cpp_obj.computeSolidForces()

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
                print(f'{self._cpp_SPFclass_name} does NOT use a body force!' )

        # self.cpp_force.setAcceleration(self.gx,self.gy,self.gz,self.damp)
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
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')

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
        conditions = ['CFL-condition', 'Gravity_waves-condition', 'Fourier-condition', 'Adami-condition']
        condition = [conditions[i] for i in np.where(C_a == C_a.max())[0]]
        C = np.sqrt(np.max(C_a))

        # Set speed of sound
        self.eos.set_speedofsound(C)

        return C, condition


    def compute_dt(self, LREF, UREF, DX, DRHO, H, MU, RHO0, COURANT=0.25):
        # Input sanity
        if LREF == 0.0:
            raise ValueError('Reference length LREF may not be zero.')
        if DRHO == 0.0:
            raise ValueError('Maximum density variation DRHO may not be zero.')
        if DX <= 0.0:
            raise ValueError('DX may not be zero or negative.')
        if H != self._param_dict['max_sl']:
            raise ValueError('Given H not equal to stored H self._param_dict[max_sl]!')
        if MU != self._param_dict['mu']:
            raise ValueError('Given MU not equal to stored MU self._param_dict[mu]!')
        if RHO0 != self.eos.RestDensity:
            raise ValueError('Given RHO0 not equal to stored RHO0 self.eos.RestDensity!')
        
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
        conditions = ['CFL-condition', 'Fourier-condition', 'Adami_max_flow-condition', 'Adami_viscous-condition' 'Gravity_waves-condition']
        condition = [conditions[i] for i in np.where(DT_a == DT_a.min())[0]]
        
        DT = COURANT * np.min(DT_a)

        return DT, condition


class TwoPhaseFlow(SPHModel):
    R""" TwoPhaseFlow solver
    """
    DENSITYMETHODS = {'SUMMATION':_sph.PhaseFlowDensityMethod.DENSITYSUMMATION,
                      'CONTINUITY':_sph.PhaseFlowDensityMethod.DENSITYCONTINUITY}

    VISCOSITYMETHODS = {'HARMONICAVERAGE':_sph.PhaseFlowViscosityMethod.HARMONICAVERAGE}

    COLORGRADIENTMETHODS = {'DENSITYRATIO':_sph.PhaseFlowColorGradientMethod.DENSITYRATIO}

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
                          artificialviscosity = bool(True), 
                          alpha = float(0.2),
                          beta = float(0.0),
                          densitydiffusion = bool(False),
                          ddiff = float(0.0),
                          shepardrenormanlization = bool(False),
                          shepardfreq = int(30),
                          compute_solid_forces = bool(False),
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
                f" This is happending since the force is moving to a new "
                f"simulation. Set a new nlist to suppress this warning.",
                RuntimeWarning)
            self.nlist = copy.deepcopy(self.nlist)
        self.nlist._attach(self._simulation)
        if isinstance(self._simulation.device, hoomd.device.CPU):
            # self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.half)
            self.nlist._cpp_obj.setStorageMode(_nsearch.NeighborList.storageMode.full)
        # TODO: understand why _nsearch.NeighborList.storageMode.half makes wierd errors!

        else:
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
                print(f'Using constant Smooting Length: {self.maxh}')

            self._cpp_obj.setConstSmoothingLength(self.maxh)
        else: 
            self.maxh      = pdata.getMaxSmoothingLength()
            if (self._simulation.device.communicator.rank == 0):
                print('Non-Constant Smooting length')
        self.rcut = kappa * self.maxh

        # print(f'self.rcut {self.rcut}')

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

        self.set_params(self.mu1, self.mu2, self.sigma12, self.omega)
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

        # self.cpp_force.setAcceleration(self.gx,self.gy,self.gz,self.damp)
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




# Dicts
Kernel = {'_WendlandC2':'WC2','_WendlandC4':'WC4','_WendlandC6':'WC6','_Quintic':'Q','_CubicSpline':'CS'}
EOS = {'_Linear':'L','_Tait':'T'}


