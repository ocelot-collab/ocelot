#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:09:20 2021

@author: andrei
"""

__author__ = "Andrei Trebushinin"


from ocelot.optics.new_wave import *
from ocelot.gui.dfl_plot import *
from ocelot.rad.optics_elements import *
from ocelot.rad.optics_line import *
from ocelot.rad.transfer_function import *
from ocelot.rad.propagation import *

ocelog.setLevel(logging.DEBUG)

from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__)

E_pohoton = 200 #central photon energy [eV]
kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton), #[m] - central wavelength
        'rho':1.0e-4, 
        'shape':(501,501,1),             #(x,y,z) shape of field matrix (reversed) to dfl.fld
        'dgrid':(400e-5,400e-5,35e-6), #(x,y,z) [m] - size of field matrix
        'power_rms':(50e-6,50e-6,3e-6),#(x,y,z) [m] - rms size of the radiation distribution (gaussian)
        'power_center':(0,0,None),     #(x,y,z) [m] - position of the radiation distribution
        'power_angle':(0,0),           #(x,y) [rad] - angle of further radiation propagation
        'power_waistpos':(0,0),        #(Z_x,Z_y) [m] downstrean location of the waist of the beam
        'wavelength':None,             #central frequency of the radiation, if different from xlamds
        'zsep':None,                   #distance between slices in z as zsep*xlamds
        'freq_chirp':0,                #dw/dt=[1/fs**2] - requency chirp of the beam around power_center[2]
        'en_pulse':None,                #total energy or max power of the pulse, use only one
        'power':1e6,
        }

l_a, l_b = 100, 10
f = l_a*l_b / (l_a + l_b)

MirrorSurface = ImperfectMirrorSurface(hrms=3e-9, angle=10*np.pi/180, plane='x', eid='MirrorSurface')
Aperture = ApertureRect(lx=4000e-6, ly=4000e-6, eid='Aperture')
# Aperture = ApertureEllips(ax=4000e-6, ay=2000e-6, eid='Aperture')
Lens = ThinLens(fx=f, fy=f, eid='Lens')

FreeSpace_before_lens = FreeSpace(l=l_a - 10, mx=2, my=2, eid='free_space_before_lens')
FreeSpace_after_aperture = FreeSpace(l=10, eid='free_space_after_aperture')
FreeSpace_to_the_sample = FreeSpace(l=l_b, mx=0.1, my=0.1, eid='free_space_to_the_sample')

dfl = generate_gaussian_dfl(**kwargs);  #Gaussian beam defenition
# dfl.to_domain('sf')
plot_dfl(dfl, domains='s', fig_name='Gaussian beam', phase=True)
# plot_dfl(dfl, domains='k', fig_name='Gaussian beam', phase=True)

#%%
line = (FreeSpace_before_lens, Aperture, 
        FreeSpace_after_aperture, MirrorSurface, Lens, FreeSpace_to_the_sample)

dfl_prop = deepcopy(dfl)
lat = OpticsLine(line, start=FreeSpace_before_lens, stop=FreeSpace_before_lens)
dfl_prop = propagate(lat, dfl_prop)                
plot_dfl(dfl_prop, fig_name='{}'.format(FreeSpace_before_lens.eid), phase=True, domains='sf')

#%%
dfl_aperture = deepcopy(dfl_prop)

# line = (Aperture, FreeSpace_after_aperture)

lat = OpticsLine(line, start=Aperture, stop=Aperture)
dfl_aperture = propagate(lat, dfl_aperture)                
plot_dfl(dfl_aperture, fig_name='{}'.format(Aperture.eid), phase=True)

lat = OpticsLine(line, start=FreeSpace_after_aperture, stop=FreeSpace_after_aperture)
dfl_aperture = propagate(lat, dfl_aperture)                
plot_dfl(dfl_aperture, domains='s', fig_name='{}'.format(FreeSpace_after_aperture.eid), phase=True)
plot_dfl(dfl_aperture, domains='k', fig_name='{}'.format(FreeSpace_after_aperture.eid), phase=True)

#%%
dfl_lens = deepcopy(dfl_aperture)

# line = (MirrorSurface, Lens)

lat = OpticsLine(line, start=MirrorSurface, stop=Lens)
dfl_lens = propagate(lat, dfl_lens)                
plot_dfl(dfl_lens, fig_name='{}'.format(Lens.eid), phase=True)

#%%
# line = (FreeSpace_to_the_sample)

dfl_sample = deepcopy(dfl_lens)
lat = OpticsLine(line, start=FreeSpace_to_the_sample, stop=FreeSpace_to_the_sample)
dfl_sample = propagate(lat, dfl_sample)
                
plot_dfl(dfl_sample, fig_name='{}'.format(FreeSpace_to_the_sample.eid), phase=True)





















