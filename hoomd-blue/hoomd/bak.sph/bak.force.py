R""" Apply forces to particles.
"""

import sys
import hoomd
from hoomd import _hoomd

## \internal
# \brief Base class for forces
#
# A force in hoomd_script reflects a ForceCompute in c++. It is responsible
# for all high-level management that happens behind the scenes for hoomd_script
# writers. 1) The instance of the c++ analyzer itself is tracked and added to the
# System 2) methods are provided for disabling the force from being added to the
# net force on each particle
class _force(hoomd.meta._metadata):
    ## \internal
    # \brief Constructs the force
    #
    # \param name name of the force instance
    #
    # Initializes the cpp_analyzer to None.
    # If specified, assigns a name to the instance
    # Assigns a name to the force in force_name;
    def __init__(self, name=None):
        # check if initialization has occured
        if not hoomd.init.is_initialized():
            hoomd.context.msg.error("Cannot create force before initialization\n");
            raise RuntimeError('Error creating force');

        # Allow force to store a name.  Used for discombobulation in the logger
        if name is None:
            self.name = "";
        else:
            self.name="_" + name;

        self.cpp_force = None;

        # increment the id counter
        id = _force.cur_id;
        _force.cur_id += 1;

        self.force_name = "force%d" % (id);
        self.enabled = True;
        self.log =True;
        hoomd.context.current.forces.append(self);

        # base class constructor
        hoomd.meta._metadata.__init__(self)

    ## \var enabled
    # \internal
    # \brief True if the force is enabled

    ## \var cpp_force
    # \internal
    # \brief Stores the C++ side ForceCompute managed by this class

    ## \var force_name
    # \internal
    # \brief The Force's name as it is assigned to the System

    ## \internal
    # \brief Checks that proper initialization has completed
    def check_initialization(self):
        # check that we have been initialized properly
        if self.cpp_force is None:
            hoomd.context.msg.error('Bug in hoomd_script: cpp_force not set, please report\n');
            raise RuntimeError();

    def disable(self, log=False):
        R""" Disable the force.

        Args:
            log (bool): Set to True if you plan to continue logging the potential energy associated with this force.

        Examples::

            force.disable()
            force.disable(log=True)

        Executing the disable command will remove the force from the simulation.
        Any :py:func:`hoomd.run()` command executed after disabling a force will not calculate or
        use the force during the simulation. A disabled force can be re-enabled
        with :py:meth:`enable()`.

        By setting *log* to True, the values of the force can be logged even though the forces are not applied
        in the simulation.  For forces that use cutoff radii, setting *log=True* will cause the correct *r_cut* values
        to be used throughout the simulation, and therefore possibly drive the neighbor list size larger than it
        otherwise would be. If *log* is left False, the potential energy associated with this force will not be
        available for logging.

        """
        hoomd.util.print_status_line();
        self.check_initialization();

        # check if we are already disabled
        if not self.enabled:
            hoomd.context.msg.warning("Ignoring command to disable a force that is already disabled");
            return;

        self.enabled = False;
        self.log = log;

        # remove the compute from the system if it is not going to be logged
        if not log:
            hoomd.context.current.system.removeCompute(self.force_name);
            hoomd.context.current.forces.remove(self)

    def enable(self):
        R""" Enable the force.

        Examples::

            force.enable()

        See :py:meth:`disable()`.
        """
        hoomd.util.print_status_line();
        self.check_initialization();

        # check if we are already disabled
        if self.enabled:
            hoomd.context.msg.warning("Ignoring command to enable a force that is already enabled");
            return;

        # add the compute back to the system if it was removed
        if not self.log:
            hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);
            hoomd.context.current.forces.append(self)

        self.enabled = True;
        self.log = True;

    def get_energy(self,group):
        R""" Get the energy of a particle group.

        Args:
            group (:py:mod:`hoomd.group`): The particle group to query the energy for.

        Returns:
            The last computed energy for the members in the group.

        Examples::

            g = group.all()
            energy = force.get_energy(g)
        """
        return self.cpp_force.calcEnergyGroup(group.cpp_group)

    ## \internal
    # \brief updates force coefficients
    def update_coeffs(self):
        pass
        raise RuntimeError("_force.update_coeffs should not be called");
        # does nothing: this is for derived classes to implement

    ## \internal
    # \brief Returns the force data
    #
    def __forces(self):
        return hoomd.data.force_data(self);

    forces = property(__forces);

    ## \internal
    # \brief Get metadata
    def get_metadata(self):
        data = hoomd.meta._metadata.get_metadata(self)
        data['enabled'] = self.enabled
        data['log'] = self.log
        if self.name is not "":
            data['name'] = self.name

        return data

# set default counter
_force.cur_id = 0;

class constant(_force):
    R""" Constant force.

    Args:
        fx (float): x-component of the force (in force units).
        fy (float): y-component of the force (in force units).
        fz (float): z-component of the force (in force units).
        group (:py:mod:`hoomd.group`): Group for which the force will be set.

    :py:class:`constant` specifies that a constant force should be added to every
    particle in the simulation or optionally to all particles in a group.

    Examples::

        force.constant(fx=1.0, fy=0.5, fz=0.25)
        const = force.constant(fx=0.4, fy=1.0, fz=0.5)
        const = force.constant(fx=0.4, fy=1.0, fz=0.5,group=fluid)
    """
    def __init__(self, fx, fy, fz, group=None):
        hoomd.util.print_status_line();

        # initialize the base class
        _force.__init__(self);

        # create the c++ mirror class
        if (group is not None):
            self.cpp_force = _hoomd.ConstForceCompute(hoomd.context.current.system_definition, group.cpp_group, fx, fy, fz);
        else:
            self.cpp_force = _hoomd.ConstForceCompute(hoomd.context.current.system_definition, hoomd.context.current.group_all.cpp_group, fx, fy, fz);

        # store metadata
        self.metadata_fields = ['fx','fy','fz']
        self.fx = fx
        self.fy = fy
        self.fz = fz
        if group is not None:
            self.metadata_fields.append('group')
            self.group = group

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

    ## Change the value of the force
    #
    # \param fx New x-component of the force (in force units)
    # \param fy New y-component of the force (in force units)
    # \param fz New z-component of the force (in force units)
    # \param group Group for which the force will be set
    #
    # Using set_force() requires that you saved the created constant force in a variable. i.e.
    # \code
    # const = force.constant(fx=0.4, fy=1.0, fz=0.5)
    # \endcode
    #
    # \b Example:
    # \code
    # const.set_force(fx=0.2, fy=0.1, fz=-0.5)
    # const.set_force(fx=0.2, fy=0.1, fz=-0.5, group=fluid)
    # \endcode
    def set_force(self, fx, fy, fz, group=None):
        self.check_initialization();
        if (group is not None):
            self.cpp_force.setGroupForce(group.cpp_group,fx,fy,fz);
        else:
            self.cpp_force.setForce(fx, fy, fz);

        self.fx = fx
        self.fy = fy
        self.fz = fz

    # there are no coeffs to update in the constant force compute
    def update_coeffs(self):
        pass
