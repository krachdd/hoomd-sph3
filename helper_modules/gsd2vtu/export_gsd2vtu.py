#!/usr/bin/env python3

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

GSD-to-VTU export library for SPH simulations.

Provides functions that read a HOOMD-SPH GSD trajectory file and write one
``.vtu`` (VTK Unstructured Grid) file per snapshot into a sub-directory
named after the GSD file (without the ``.gsd`` extension).  Each snapshot
is exported as a point cloud with the particle fields listed below.

Available export functions and the fields they include:

===============  ===================================================================
Function         Fields exported
===============  ===================================================================
export_basic     position, velocity, typeid, slength, mass, density, pressure, energy
export_spf       + fictitious solid velocity (aux1)
export_tvspf     + fictitious solid velocity (aux1), back-pressure correction (aux2),
                 transport velocity (aux3)
export_tpf       + fictitious solid velocity (aux1), solid-fluid normal (aux2),
                 fluid-fluid normal (aux3), surface-force density (aux4)
export_gdgd      + fictitious solid velocity (aux1), scalar field T (aux4.x)
export_fs        + fictitious solid velocity (aux1), back-pressure correction (aux2),
                 transport velocity (aux3), kernel completeness :math:`\lambda` (aux4.x),
                 curvature :math:`\kappa` (aux4.y)
export_tpftv     + fictitious solid velocity (aux1), back-pressure correction (aux2),
                 transport velocity (aux3), surface-force density (aux4)
export_all       + all four auxiliary vectors (raw labels)
===============  ===================================================================

Usage example::

    from gsd2vtu import export_gsd2vtu as exp
    exp.export_spf('results/my_simulation.gsd')
"""

# --- HEADER ---------------------------------------------------
import gsd.fl
import gsd.hoomd
import gsd.pygsd
import numpy as np
import sys
import os
from pyevtk.hl import pointsToVTK as vtk
#--------------------------------------------------------------

def import_trajectory(GSDfilename):
    """Open a HOOMD-SPH GSD file and return its trajectory object.

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` file.

    Returns
    -------
    gsd.hoomd.HOOMDTrajectory
        Trajectory object that can be iterated over to access individual
        snapshots.
    """
    f = gsd.fl.GSDFile(name=GSDfilename, mode='r+', application="HOOMD-SPH",
                       schema="hoomd", schema_version=[1, 0])
    return gsd.hoomd.HOOMDTrajectory(f)


def export_basic(GSDfilename):
    """Export basic particle fields from a GSD trajectory to VTU files.

    Writes one ``.vtu`` file per snapshot containing: position, velocity
    (x/y/z), type ID, smoothing length, mass, density, pressure, and energy.

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # print(t[0].particles.position)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                    np.array(snapshot.particles.position.T[1]),
                    np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'     :np.array(snapshot.particles.typeid),
                    'Slength'    :np.array(snapshot.particles.slength),
                    'Mass'       :np.array(snapshot.particles.mass),
                    'Density'    :np.array(snapshot.particles.density),
                    'Pressure'   :np.array(snapshot.particles.pressure),
                    'Energy'     :np.array(snapshot.particles.energy),
                    #'Aux1x'      :np.array(snapshot.particles.auxiliary1.T[0]),
                    #'Aux1y'      :np.array(snapshot.particles.auxiliary1.T[1]),
                    #'Aux1z'      :np.array(snapshot.particles.auxiliary1.T[2]),
                    #'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    #'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    #'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    #'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    #'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    #'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    #'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    #'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    #'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )
    

def export_spf(GSDfilename):
    """Export single-phase flow fields from a GSD trajectory to VTU files.

    Extends :func:`export_basic` with the fictitious solid-boundary velocity
    stored in ``auxiliary1`` (Adami 2012 no-slip boundary condition).

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                      np.array(snapshot.particles.position.T[1]),
                      np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x'          :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y'          :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z'          :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'              :np.array(snapshot.particles.typeid),
                    'Slength'             :np.array(snapshot.particles.slength),
                    'Mass'                :np.array(snapshot.particles.mass),
                    'Density'             :np.array(snapshot.particles.density),
                    'Pressure'            :np.array(snapshot.particles.pressure),
                    'Energy'              :np.array(snapshot.particles.energy),
                    'Ficticious Velx'     :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Ficticious Vely'     :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Ficticious Velz'     :np.array(snapshot.particles.auxiliary1.T[2]),
                    # 'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    # 'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    # 'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    # 'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    # 'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    # 'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    # 'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    # 'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    # 'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )


def export_tvspf(GSDfilename):
    """Export transport-velocity single-phase flow fields to VTU files.

    Extends :func:`export_spf` with the back-pressure correction vector
    (``auxiliary2``) and transport velocity (``auxiliary3``) written by the
    ``SinglePhaseFlowTV`` (Adami 2013) integrator.

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                      np.array(snapshot.particles.position.T[1]),
                      np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x'          :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y'          :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z'          :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'              :np.array(snapshot.particles.typeid),
                    'Slength'             :np.array(snapshot.particles.slength),
                    'Mass'                :np.array(snapshot.particles.mass),
                    'Density'             :np.array(snapshot.particles.density),
                    'Pressure'            :np.array(snapshot.particles.pressure),
                    'Energy'              :np.array(snapshot.particles.energy),
                    'Ficticious Velx'     :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Ficticious Vely'     :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Ficticious Velz'     :np.array(snapshot.particles.auxiliary1.T[2]),
                    'Backpressure corrx'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    'Backpressure corry'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    'Backpressure corrz'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    'Transport Velx'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    'Transport Vely'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    'Transport Velz'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    # 'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    # 'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    # 'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )


def export_tpf(GSDfilename):
    """Export two-phase flow fields from a GSD trajectory to VTU files.

    Exports all four auxiliary vectors with physically meaningful labels:

    - ``auxiliary1`` — fictitious solid-boundary velocity (Adami 2012)
    - ``auxiliary2`` — solid-fluid interface normal vector
    - ``auxiliary3`` — fluid-fluid interface normal vector (colour gradient)
    - ``auxiliary4`` — surface-tension force density

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                      np.array(snapshot.particles.position.T[1]),
                      np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x'                :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y'                :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z'                :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'                    :np.array(snapshot.particles.typeid),
                    'Slength'                   :np.array(snapshot.particles.slength),
                    'Mass'                      :np.array(snapshot.particles.mass),
                    'Density'                   :np.array(snapshot.particles.density),
                    'Pressure'                  :np.array(snapshot.particles.pressure),
                    'Energy'                    :np.array(snapshot.particles.energy),
                    'Ficticious Velocity x'     :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Ficticious Velocity y'     :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Ficticious Velocity z'     :np.array(snapshot.particles.auxiliary1.T[2]),
                    'Solid normal-vector x'     :np.array(snapshot.particles.auxiliary2.T[0]),
                    'Solid normal-vector y'     :np.array(snapshot.particles.auxiliary2.T[1]),
                    'Solid normal-vector z'     :np.array(snapshot.particles.auxiliary2.T[2]),
                    'Fluid-Fluid n-vector x'    :np.array(snapshot.particles.auxiliary3.T[0]),
                    'Fluid-Fluid n-vector y'    :np.array(snapshot.particles.auxiliary3.T[1]),
                    'Fluid-Fluid n-vector z'    :np.array(snapshot.particles.auxiliary3.T[2]),
                    'Surface Force density x'   :np.array(snapshot.particles.auxiliary4.T[0]),
                    'Surface Force density y'   :np.array(snapshot.particles.auxiliary4.T[1]),
                    'Surface Force density z'   :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )

def export_all(GSDfilename):
    """Export all particle fields from a GSD trajectory to VTU files.

    Exports all basic fields and all four auxiliary vectors using generic
    labels (``Aux1x/y/z`` … ``Aux4x/y/z``).  Useful for debugging or when
    the physical meaning of the auxiliary arrays is not yet determined.

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename = GSDfilename)

    # Run loop over all snapshots
    count = 0
    for snapshot in t:
        count += 1
       
        pname = GSDfilename.replace('.gsd','')
        
        if not os.path.exists(pname):
            os.makedirs(pname)
            
        # Define VTU export filename
        filename = pname+'/'+pname+'_'+str(f'{snapshot.configuration.step:09d}')
        
        vtk(filename, np.array(snapshot.particles.position.T[0]),
                    np.array(snapshot.particles.position.T[1]),
                    np.array(snapshot.particles.position.T[2]),
            data = {'Velocity x' :np.array(snapshot.particles.velocity.T[0]),
                    'Velocity y' :np.array(snapshot.particles.velocity.T[1]),
                    'Velocity z' :np.array(snapshot.particles.velocity.T[2]),
                    'TypeId'     :np.array(snapshot.particles.typeid),
                    'Slength'    :np.array(snapshot.particles.slength),
                    'Mass'       :np.array(snapshot.particles.mass),
                    'Density'    :np.array(snapshot.particles.density),
                    'Pressure'   :np.array(snapshot.particles.pressure),
                    'Energy'     :np.array(snapshot.particles.energy),
                    'Aux1x'      :np.array(snapshot.particles.auxiliary1.T[0]),
                    'Aux1y'      :np.array(snapshot.particles.auxiliary1.T[1]),
                    'Aux1z'      :np.array(snapshot.particles.auxiliary1.T[2]),
                    'Aux2x'      :np.array(snapshot.particles.auxiliary2.T[0]),
                    'Aux2y'      :np.array(snapshot.particles.auxiliary2.T[1]),
                    'Aux2z'      :np.array(snapshot.particles.auxiliary2.T[2]),
                    'Aux3x'      :np.array(snapshot.particles.auxiliary3.T[0]),
                    'Aux3y'      :np.array(snapshot.particles.auxiliary3.T[1]),
                    'Aux3z'      :np.array(snapshot.particles.auxiliary3.T[2]),
                    'Aux4x'      :np.array(snapshot.particles.auxiliary4.T[0]),
                    'Aux4y'      :np.array(snapshot.particles.auxiliary4.T[1]),
                    'Aux4z'      :np.array(snapshot.particles.auxiliary4.T[2]),
                      },
                  )


def export_gdgd(GSDfilename):
    """Export density-gradient-driven flow (GDGD) fields to VTU files.

    Extends :func:`export_spf` with the transported scalar field (temperature
    or concentration) stored in ``auxiliary4.x``.  The remaining components
    of ``auxiliary4`` (y, z) are not used by ``SinglePhaseFlowGDGD`` and are
    not exported.

    Auxiliary-field layout (``SinglePhaseFlowGDGD``):

    - ``auxiliary1`` — fictitious solid-boundary velocity (Adami 2012)
    - ``auxiliary4.x`` — scalar field T (temperature / concentration)

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename=GSDfilename)

    for snapshot in t:
        pname = GSDfilename.replace('.gsd', '')

        if not os.path.exists(pname):
            os.makedirs(pname)

        filename = pname + '/' + pname + '_' + str(f'{snapshot.configuration.step:09d}')

        vtk(filename,
            np.array(snapshot.particles.position.T[0]),
            np.array(snapshot.particles.position.T[1]),
            np.array(snapshot.particles.position.T[2]),
            data={
                'Velocity x'       : np.array(snapshot.particles.velocity.T[0]),
                'Velocity y'       : np.array(snapshot.particles.velocity.T[1]),
                'Velocity z'       : np.array(snapshot.particles.velocity.T[2]),
                'TypeId'           : np.array(snapshot.particles.typeid),
                'Slength'          : np.array(snapshot.particles.slength),
                'Mass'             : np.array(snapshot.particles.mass),
                'Density'          : np.array(snapshot.particles.density),
                'Pressure'         : np.array(snapshot.particles.pressure),
                'Energy'           : np.array(snapshot.particles.energy),
                'Ficticious Velx'  : np.array(snapshot.particles.auxiliary1.T[0]),
                'Ficticious Vely'  : np.array(snapshot.particles.auxiliary1.T[1]),
                'Ficticious Velz'  : np.array(snapshot.particles.auxiliary1.T[2]),
                'Scalar T'         : np.array(snapshot.particles.auxiliary4.T[0]),
            })


def export_fs(GSDfilename):
    """Export free-surface flow (FS) fields to VTU files.

    Exports the fields written by ``SinglePhaseFlowFS`` (free-surface variant
    of the transport-velocity single-phase solver).  At GSD write time the
    auxiliary slots contain:

    - ``auxiliary1`` — fictitious solid-boundary velocity (Adami 2012)
    - ``auxiliary2`` — back-pressure correction vector (BPC); the
      free-surface outward normal that temporarily occupies this slot during
      ``detect_freesurface()`` has already been overwritten by ``forcecomputation()``
    - ``auxiliary3`` — transport velocity
    - ``auxiliary4.x`` — kernel completeness :math:`\lambda` (Shepard sum;
      :math:`\lambda < ` ``fs_threshold`` identifies free-surface particles)
    - ``auxiliary4.y`` — mean curvature :math:`\kappa` (computed only for surface particles;
      zero for bulk particles)

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename=GSDfilename)

    for snapshot in t:
        pname = GSDfilename.replace('.gsd', '')

        if not os.path.exists(pname):
            os.makedirs(pname)

        filename = pname + '/' + pname + '_' + str(f'{snapshot.configuration.step:09d}')

        vtk(filename,
            np.array(snapshot.particles.position.T[0]),
            np.array(snapshot.particles.position.T[1]),
            np.array(snapshot.particles.position.T[2]),
            data={
                'Velocity x'         : np.array(snapshot.particles.velocity.T[0]),
                'Velocity y'         : np.array(snapshot.particles.velocity.T[1]),
                'Velocity z'         : np.array(snapshot.particles.velocity.T[2]),
                'TypeId'             : np.array(snapshot.particles.typeid),
                'Slength'            : np.array(snapshot.particles.slength),
                'Mass'               : np.array(snapshot.particles.mass),
                'Density'            : np.array(snapshot.particles.density),
                'Pressure'           : np.array(snapshot.particles.pressure),
                'Energy'             : np.array(snapshot.particles.energy),
                'Ficticious Velx'    : np.array(snapshot.particles.auxiliary1.T[0]),
                'Ficticious Vely'    : np.array(snapshot.particles.auxiliary1.T[1]),
                'Ficticious Velz'    : np.array(snapshot.particles.auxiliary1.T[2]),
                'Backpressure corrx' : np.array(snapshot.particles.auxiliary2.T[0]),
                'Backpressure corry' : np.array(snapshot.particles.auxiliary2.T[1]),
                'Backpressure corrz' : np.array(snapshot.particles.auxiliary2.T[2]),
                'Transport Velx'     : np.array(snapshot.particles.auxiliary3.T[0]),
                'Transport Vely'     : np.array(snapshot.particles.auxiliary3.T[1]),
                'Transport Velz'     : np.array(snapshot.particles.auxiliary3.T[2]),
                'Lambda FS'          : np.array(snapshot.particles.auxiliary4.T[0]),
                'Curvature kappa'    : np.array(snapshot.particles.auxiliary4.T[1]),
            })


def export_tpftv(GSDfilename):
    """Export two-phase flow with transport velocity (TPFTV) fields to VTU files.

    ``TwoPhaseFlowTV`` uses the transport-velocity (Adami 2013) formulation,
    which changes the meaning of ``auxiliary2`` and ``auxiliary3`` relative to
    the standard ``TwoPhaseFlow`` (see :func:`export_tpf`).  At GSD write
    time the auxiliary slots contain:

    - ``auxiliary1`` — fictitious solid-boundary velocity (Adami 2012)
    - ``auxiliary2`` — back-pressure correction vector (BPC); replaces the
      solid-fluid normal vector used in the non-TV two-phase solver
    - ``auxiliary3`` — transport velocity; replaces the fluid-fluid colour-
      gradient normal used in the non-TV solver
    - ``auxiliary4`` — surface-tension force density (same as non-TV TPF)

    Parameters
    ----------
    GSDfilename : str
        Path to the ``.gsd`` input file.

    Returns
    -------
    None
    """
    print(f'{os.path.basename(__file__)}: Export {GSDfilename} to .vtu')

    t = import_trajectory(GSDfilename=GSDfilename)

    for snapshot in t:
        pname = GSDfilename.replace('.gsd', '')

        if not os.path.exists(pname):
            os.makedirs(pname)

        filename = pname + '/' + pname + '_' + str(f'{snapshot.configuration.step:09d}')

        vtk(filename,
            np.array(snapshot.particles.position.T[0]),
            np.array(snapshot.particles.position.T[1]),
            np.array(snapshot.particles.position.T[2]),
            data={
                'Velocity x'              : np.array(snapshot.particles.velocity.T[0]),
                'Velocity y'              : np.array(snapshot.particles.velocity.T[1]),
                'Velocity z'              : np.array(snapshot.particles.velocity.T[2]),
                'TypeId'                  : np.array(snapshot.particles.typeid),
                'Slength'                 : np.array(snapshot.particles.slength),
                'Mass'                    : np.array(snapshot.particles.mass),
                'Density'                 : np.array(snapshot.particles.density),
                'Pressure'                : np.array(snapshot.particles.pressure),
                'Energy'                  : np.array(snapshot.particles.energy),
                'Ficticious Velocity x'   : np.array(snapshot.particles.auxiliary1.T[0]),
                'Ficticious Velocity y'   : np.array(snapshot.particles.auxiliary1.T[1]),
                'Ficticious Velocity z'   : np.array(snapshot.particles.auxiliary1.T[2]),
                'Backpressure corr x'     : np.array(snapshot.particles.auxiliary2.T[0]),
                'Backpressure corr y'     : np.array(snapshot.particles.auxiliary2.T[1]),
                'Backpressure corr z'     : np.array(snapshot.particles.auxiliary2.T[2]),
                'Transport Vel x'         : np.array(snapshot.particles.auxiliary3.T[0]),
                'Transport Vel y'         : np.array(snapshot.particles.auxiliary3.T[1]),
                'Transport Vel z'         : np.array(snapshot.particles.auxiliary3.T[2]),
                'Surface Force density x' : np.array(snapshot.particles.auxiliary4.T[0]),
                'Surface Force density y' : np.array(snapshot.particles.auxiliary4.T[1]),
                'Surface Force density z' : np.array(snapshot.particles.auxiliary4.T[2]),
            })