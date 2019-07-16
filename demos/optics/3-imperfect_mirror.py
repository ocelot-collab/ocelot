# -*- coding: utf-8 -*-
"""
Created on Mon Oct 8 2018
@author: Mykola Veremchuk
"""

import matplotlib.pyplot as plt
import numpy as np
from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import plot_dfl, plot_1d_hprofile
ocelog.setLevel(logging.DEBUG)
# from copy import deepcopy

__author__ = "Svitozar Serkez, Mykola Veremchuk"



hprofile = generate_1d_profile(hrms=1e-9,              #
                               length=0.03,             #
                               points_number=1111,      #
                               wavevector_cutoff=0,     #
                               k=None,                  #
                               psd=None,                #
                               seed=666)

plot_1d_hprofile(hprofile, fig_name='mirror1 height profile and PSD')

dfl = generate_gaussian_dfl(1e-9, (1000, 1000, 1))
plot_dfl(dfl, phase=1, fig_name='radiation before mirror1')


dfl_reflect_surface(dfl, 
                    angle=np.pi * 2 / 180, 
                    hrms=None, 
                    height_profile=hprofile, 
                    axis='x')

plot_dfl(dfl, phase=1, fig_name='radiation after reflection from mirror1')


dfl.prop(z=10)
plot_dfl(dfl, phase=1, fig_name='radiation after reflection from mirror1 and propagation')


#%%

hprofile2 = generate_1d_profile(hrms=1e-9,              #
                               length=0.005,             #
                               points_number=150,      #
                               wavevector_cutoff=2000,     #
                               k=None,                  #
                               psd=None,                #
                               seed=123)

plot_1d_hprofile(hprofile2, fig_name='mirror2 height profile and PSD')


dfl2 = generate_gaussian_dfl(1e-9, (1000, 1000, 1))

dfl_reflect_surface(dfl2, 
                    angle=np.pi * 1 / 180, 
                    hrms=None, 
                    height_profile=hprofile2, 
                    axis='y')

plot_dfl(dfl2, phase=1, fig_name='radiation after mirror2')
dfl2.prop(z=10)
plot_dfl(dfl2, phase=1, fig_name='radiation after reflection from mirror2 and propagation')

#%%

dfl3 = generate_gaussian_dfl(1e-9, (1000, 1000, 1))
dfl_reflect_surface(dfl3, 
                    angle=np.pi * 2 / 180, 
                    hrms=1e-9, 
                    axis='x')
dfl3.prop(z=10)
plot_dfl(dfl3, fig_name='radiation after reflection internally generated mirror and propagation')


#TODO commented examples of write/read
