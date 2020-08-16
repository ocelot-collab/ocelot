# -*- coding: utf-8 -*-
"""
Created on Mon Oct 8 2018
@author: Mykola Veremchuk
Examples of converting synchrotron radiation (ocelot.rad.screen.Screen object) to radiation field
(ocelot.optics.wave.RadiationField object) which can be used with ocelot optics features
"""
__author__ = "Mykola Veremchuk"

import numpy as np
from ocelot.cpbd.elements import Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import Beam
from ocelot.rad.screen import Screen
from ocelot.rad.radiation_py import calculate_radiation
from ocelot.optics.wave import dfl_waistscan, screen2dfl, RadiationField
from ocelot.gui.dfl_plot import plot_dfl, plot_dfl_waistscan

# %%
# generating 2D synchrotron radiation (it will take about 1-3 minute)
# LOOK TUTORIAL ABOUT GENERATING SYNCHROTRON RADIATION IN demos/ipython_tutorials/9_synchrotron_radiation.ipynb
und = Undulator(Kx=0.43, nperiods=500, lperiod=0.007, eid="und")
lat = MagneticLattice(und)
beam = Beam()
beam.E = 2.5  # beam energy in [GeV]
beam.I = 0.1  # beam current in [A]
screen_2d = Screen()
screen_2d.z = 50.0  # distance from the begining of lattice to the screen
screen_2d.size_x = 0.0015  # half of screen size in [m] in horizontal plane
screen_2d.size_y = 0.0015  # half of screen size in [m] in vertical plane

screen_2d.ny = 200
screen_2d.nx = 200
screen_2d.start_energy = 7762  # [eV], starting photon energy
screen_2d.end_energy = 7762  # [eV], ending photon energy
screen_2d.num_energy = 1  # number of energy points[eV]

# calculate radiation
screen_2d = calculate_radiation(lat, screen_2d, beam)

# %%
# converting Screen to RadiationField (function generates new RadiationField() without changing Screen)
dfl_2d = screen2dfl(screen_2d,          # Screen object, electric field of which will be used to generate RadiationField
                    polarization='x')   # polarization for conversion to RadiationField ('x' or 'y')
plot_dfl(dfl_2d,
         domains='fs',
         fig_name='dfl_2d generated from screen_2d')

# scanning for waist position
wfs = dfl_waistscan(dfl_2d, np.linspace(-80, -20, 200))
plot_dfl_waistscan(wfs, fig_name='waist scan of dfl_2d')

# half analytical propagation to waist point
dfl_2d.prop_m(-48.25,
              m=0.05)
plot_dfl(dfl_2d,
         domains='fs',
         fig_name='dfl_2d at waist position')

# %%
# generating 3D synchrotron radiation (it will take up to 5 minute)
# LOOK TUTORIAL ABOUT GENERATING SYNCHROTRON RADIATION IN demos/ipython_tutorials/9_synchrotron_radiation.ipynb
screen_3d = Screen()
screen_3d.z = 50.0  # distance from the begining of lattice to the screen
screen_3d.size_x = 0.0015  # half of screen size in [m] in horizontal plane
screen_3d.size_y = 0.0015  # half of screen size in [m] in vertical plane

screen_3d.ny = 75
screen_3d.nx = 75
screen_3d.start_energy = 7730  # [eV], starting photon energy
screen_3d.end_energy = 7795  # [eV], ending photon energy
screen_3d.num_energy = 20  # number of energy points[eV]

# calculate radiation
screen_3d = calculate_radiation(lat, screen_3d, beam)

# %%
# converting Screen to RadiationField (function generates new RadiationField() without changing Screen)
dfl_3d = screen2dfl(screen_3d,          # Screen object, electric field of which will be used to generate RadiationField
                    polarization='x')   # polarization for conversion to RadiationField ('x' or 'y')
plot_dfl(dfl_3d,
         domains='fs',
         fig_name='dfl_3d generated from screen_3d in frequency-space domains',
         slice_xy=True)    # bool type variable, if True, slices will be plotted; if False, projections will be plotted
plot_dfl(dfl_3d,
         domains='ts',
         fig_name='dfl_3d generated from screen_3d in time-space domains',
         slice_xy=False)    # bool type variable, if True, slices will be plotted; if False, projections will be plotted

# scanning for waist position
wfs = dfl_waistscan(dfl_3d,
                    z_pos=np.linspace(-80, -20, 200))
plot_dfl_waistscan(wfs, fig_name='waist scan of dfl_3d')

# half analytical propagation to waist point
dfl_3d.prop_m(-48.25,
              m=0.05)
plot_dfl(dfl_3d,
         domains='fs',
         fig_name='dfl_3d at waist position in frequency-space domains',
         slice_xy=False)     # bool type variable, if True, slices will be plotted; if False, projections will be plotted
