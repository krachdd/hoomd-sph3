# coding: utf-8

# Copyright (c) 2009-2016 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

R""" sph - Integrators.
"""

import copy
import sys
import hoomd
from hoomd import _hoomd
from hoomd.sph import _sph
from hoomd.integrate import _integrator, _integration_method

class BaseIntegrator(_integrator):
    R""" Enables a variety of standard integration methods.

    Args:
        dt (float): Each time step of the simulation :py:func:`hoomd.run()` will advance the real time of the system forward by *dt* (in time units).

    :py:class:`BaseIntegrator` performs a standard time step integration technique to move the system forward. At each time
    step, all of the specified forces are evaluated and used in moving the system forward to the next step.

    By itself, :py:class:`BaseIntegrator` does nothing. You must specify one or more integration methods to apply to the
    system. Each integration method can be applied to only a specific group of particles enabling advanced simulation
    techniques.

    The following commands can be used to specify the integration methods used by integrate.BaseIntegrator.

    - :py:class:`VelocityVerlet`

    There can only be one integration mode active at a time. If there are more than one ``integrate.BaseIntegrator`` commands in
    a hoomd script, only the most recent before a given :py:func:`run()` will take effect.

    Examples::

        integrate.BaseIntegrator(dt=0.005)
        integrator_mode = integrate.BaseIntegrator(dt=0.001)
    """
    def __init__(self, dt):
        hoomd.util.print_status_line();

        # initialize base class
        _integrator.__init__(self);

        # Store metadata
        self.dt = dt
        self.metadata_fields = ['dt']

        # initialize the reflected c++ class
        self.cpp_integrator = _sph.SPHIntegratorTwoStep(hoomd.context.current.system_definition, dt);
        self.supports_methods = True;

        hoomd.context.current.system.setIntegrator(self.cpp_integrator);

        hoomd.util.quiet_status();
        hoomd.util.unquiet_status();

    def set_params(self, dt=None):
        R""" Changes parameters of an existing integration mode.

        Args:
            dt (float): New time step delta (if set) (in time units).

        Examples::

            integrator_mode.set_params(dt=0.007)
        """
        hoomd.util.print_status_line();
        self.check_initialization();

        # change the parameters
        if dt is not None:
            self.dt = dt
            self.cpp_integrator.setDeltaT(dt);

class VelocityVerlet(_integration_method):
    R""" Velocity-Verlet time integrator.

    Args:
        dt (float): Each time step of the simulation :py:func:`hoomd.run()` will advance the real time of the system forward by *dt* (in time units).

    :py:class:`VelocityVerlet` performs a two-step time step integration technique to move the system forward. At each time
    step, all of the specified forces are evaluated and used in moving the system forward to the next step.
    """
    def __init__(self, group):
        hoomd.util.print_status_line();

        # initialize base class
        _integration_method.__init__(self);

        # create the compute thermo
        hoomd.compute._get_unique_thermo(group=group);

        # initialize the reflected c++ class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_method = _sph.VelocityVerlet(hoomd.context.current.system_definition, group.cpp_group, False);
        else:
            self.cpp_method = _sph.VelocityVerletGPU(hoomd.context.current.system_definition, group.cpp_group, False);

        self.cpp_method.validateGroup()

        # Store metadata
        self.group = group
        self.metadata_fields = ['group']

class SuspendedObjectIntegrator(_integration_method):
    R""" SuspendedObjectIntegrator time integrator.

    Args:
        group (:py:mod:`hoomd.group`): Group of particles on which to apply this method.

    :py:class:`SuspendedObjectIntegrator` performs a two-step time step integration technique on suspended object to
    integrate them forward in time
    """
    def __init__(self, group):
        hoomd.util.print_status_line();

        # initialize base class
        _integration_method.__init__(self);

        # create the compute thermo
        hoomd.compute._get_unique_thermo(group=group);

        # initialize the reflected c++ class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_method = _sph.SuspendedObjectIntegrator(hoomd.context.current.system_definition, group.cpp_group, False);
        else:
            self.cpp_method = _sph.SuspendedObjectIntegratorGPU(hoomd.context.current.system_definition, group.cpp_group, False);


        self.cpp_method.validateGroup()

        # Store metadata
        self.group = group
        self.metadata_fields = ['group']

class RigidBodyIntegrator(_integration_method):
    R""" RigidBodyIntegrator time integrator.

    Args:
        group (:py:mod:`hoomd.group`): Group of particles on which to apply this method.
        transvel (List of 3 :py:mod:`hoomd.variant` or :py:obj:`float`): Translation velocity
        rotvel (:py:mod:`hoomd.variant` or :py:obj:`float`): Rotational velocity
        pivotpnt (List of 3 :py:obj:`float`): Rotational pivot point
        rotaxis (List of 3 :py:obj:`float`): Rotational axis

    :py:class:`RigidBodyIntegrator` performs a two-step time step integration technique on rigid bodies to
    integrate them forward in time
    """
    def __init__(self, group, transvel, rotvel, pivotpnt, rotaxis):
        hoomd.util.print_status_line();

        # initialize base class
        _integration_method.__init__(self);

        # check input sanity
        if len(transvel) != 3:
            raise ValueError('transvel need to be a list of 3 hoomd.variants or floats')
        if len(pivotpnt) != 3:
            raise ValueError('pivotpnt need to be a list of 3 floats')
        if len(rotaxis) != 3:
            raise ValueError('rotaxis need to be a list of 3 floats')

        # setup the variant inputs
        transvel_x = hoomd.variant._setup_variant_input(transvel[0]);
        transvel_y = hoomd.variant._setup_variant_input(transvel[1]);
        transvel_z = hoomd.variant._setup_variant_input(transvel[2]);
        rotvel = hoomd.variant._setup_variant_input(rotvel);

        # create the compute thermo
        hoomd.compute._get_unique_thermo(group=group);

        # Store metadata
        self.group = group
        self.transvel_x = transvel_x
        self.transvel_y = transvel_y
        self.transvel_z = transvel_z
        self.rotvel     = rotvel
        self.pivotpnt_x = pivotpnt[0]
        self.pivotpnt_y = pivotpnt[1]
        self.pivotpnt_z = pivotpnt[2]
        self.rotaxis_x  = rotaxis[0]
        self.rotaxis_y  = rotaxis[1]
        self.rotaxis_z  = rotaxis[2]
        self.metadata_fields = ['group']

        # initialize the reflected c++ class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_method = _sph.RigidBodyIntegrator(hoomd.context.current.system_definition, group.cpp_group,
                                                       transvel_x.cpp_variant, transvel_y.cpp_variant, transvel_z.cpp_variant,
                                                       rotvel.cpp_variant, pivotpnt[0], pivotpnt[1], pivotpnt[2],
                                                       rotaxis[0], rotaxis[1], rotaxis[2], False);
        else:
            self.cpp_method = _sph.RigidBodyIntegratorGPU(hoomd.context.current.system_definition, group.cpp_group,
                                                       transvel_x.cpp_variant, transvel_y.cpp_variant, transvel_z.cpp_variant,
                                                       rotvel.cpp_variant, pivotpnt[0], pivotpnt[1], pivotpnt[2],
                                                       rotaxis[0], rotaxis[1], rotaxis[2], False);


        self.cpp_method.validateGroup()

    def set_params(self, transvel=None, rotvel=None, pivotpnt=None, rotaxis=None):
        R""" Changes parameters of an existing integrator.

        Args:
            transvel (list of floats or variants): New translational velocity
            rotvel (float or variant): New rotational velocity
            pivotpnt (list of floats): New rotational pivot point
            rotaxis (list of floats): New rotation axis
        """
        hoomd.util.print_status_line();
        self.check_initialization();

        # change the parameters
        if transvel is not None:
            if len(transvel) != 3:
               raise ValueError('transvel need to be a list of 3 hoomd.variants or floats')
            transvel_x = hoomd.variant._setup_variant_input(transvel[0]);
            transvel_y = hoomd.variant._setup_variant_input(transvel[1]);
            transvel_y = hoomd.variant._setup_variant_input(transvel[2]);
            self.cpp_method.setTranslationalVelocity(transvel_x.cpp_variant, transvel_y.cpp_variant, transvel_z.cpp_variant);
            self.transvel_x = transvel_x
            self.transvel_y = transvel_y
            self.transvel_z = transvel_z
        if rotvel is not None:
            rotvel = hoomd.variant._setup_variant_input(rotvel);
            self.cpp_method.setRotationSpeed(rotvel.cpp_variant);
            self.rotvel = rotvel
        if pivotpnt is not None:
            if len(pivotpnt) != 3:
                raise ValueError('pivotpnt need to be a list of 3 floats')
            self.cpp_method.setPivotPoint(pivotpnt[0], pivotpnt[1], pivotpnt[2]);
            self.pivotpnt_x = pivotpnt[0]
            self.pivotpnt_y = pivotpnt[1]
            self.pivotpnt_z = pivotpnt[2]
        if rotaxis is not None:
            if len(rotaxis) != 3:
                raise ValueError('rotaxis need to be a list of 3 floats')
            self.cpp_method.setRotationAxis(rotaxis[0], rotaxis[1], rotaxis[2]);
            self.rotaxis_x  = rotaxis[0]
            self.rotaxis_y  = rotaxis[1]
            self.rotaxis_z  = rotaxis[2]
