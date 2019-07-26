# -*- coding: utf-8 -*-
"""
Created on Mon Oct 8 2018
@author: Svitozar Serkez, Mykola Veremchuk
Examples of generating imperfect highly polished mirror surface, and its usage for reflection of RadiationField from
the mirror surface considering effects of mirror surface height errors
"""
__author__ = "Svitozar Serkez, Mykola Veremchuk"

import matplotlib.pyplot as plt
import numpy as np
from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import plot_dfl, plot_1d_hprofile
ocelog.setLevel(logging.DEBUG)
# from copy import deepcopy


# generating highly polished mirror by height errors RMS
hprofile = generate_1d_profile(hrms=1e-9,               # [m] height errors root mean square
                               length=0.03,             # [m] length of the mirror surface
                               points_number=1111,      # number of points (pixels) at the mirror surface
                               wavevector_cutoff=0,     # [1/m] point on k axis for cut off large wave lengths in the PSD (with default value 0 effects on nothing)
                               k=None,                  # [1/m] one dimension numpy.ndarray of wave vectors; if specified, psd will be calculated using this values as arguments
                               psd=None,                # [m^3] 1d array; power spectral density of surface (if not specified, will be generated)
                               seed=666)                # seed for np.random.seed() to allow reproducibility

# plotting 1d height profile
plot_1d_hprofile(hprofile, fig_name='mirror1 height profile and PSD')

# generating gaussian RadiationField
dfl = generate_gaussian_dfl(1e-9, (1000, 1000, 1))

# plotting generated RadiationField
plot_dfl(dfl, phase=1, fig_name='radiation before mirror1')

# reflecting generated RadiationField from the imperfect mirror
dfl_reflect_surface(dfl,                        # ocelot.optics.wave.RadiationField, which will be reflected from imperfect mirror surface (hprofile2)
                    angle=np.pi * 2 / 180,      # [radians] angle of incidence with respect to the surface
                    hrms=None,                  # [m] height errors root mean square
                    height_profile=hprofile,    # HeightProfile object of the reflecting surface (if not specified, will be generated using hrms)
                    axis='x')                   # direction along which reflection takes place

# plotting RadiationField after reflection
plot_dfl(dfl, phase=1, fig_name='radiation after reflection from mirror1')

# propagating RadiationField for 10 meters
dfl.prop(z=10)

# plotting RadiationField after propagation
plot_dfl(dfl, phase=1, fig_name='radiation after reflection from mirror1 and propagation')


#%%

# generating highly polished mirror by height errors RMS
hprofile2 = generate_1d_profile(hrms=1e-9,                  # [m] height errors root mean square
                               length=0.005,                # [m] length of the mirror surface
                               points_number=150,           # number of points (pixels) at the mirror surface
                               wavevector_cutoff=2000,      # [1/m] point on k axis for cut off large wave lengths in the PSD (with default value 0 effects on nothing)
                               k=None,                      # [1/m] one dimension numpy.ndarray of wave vectors; if specified, psd will be calculated using this values as arguments
                               psd=None,                    # [m^3] 1d array; power spectral density of surface (if not specified, will be generated)
                               seed=123)                    # seed for np.random.seed() to allow reproducibility

plot_1d_hprofile(hprofile2, fig_name='mirror2 height profile and PSD')


dfl2 = generate_gaussian_dfl(1e-9, (1000, 1000, 1))

dfl_reflect_surface(dfl2,                           # ocelot.optics.wave.RadiationField, which will be reflected from imperfect mirror surface (hprofile2)
                    angle=np.pi * 1 / 180,          # [radians] angle of incidence with respect to the surface
                    hrms=None,                      # [m] height errors root mean square
                    height_profile=hprofile2,       # HeightProfile object of the reflecting surface (if not specified, will be generated using hrms)
                    axis='y')                       # direction along which reflection takes place

plot_dfl(dfl2, phase=1, fig_name='radiation after mirror2')
dfl2.prop(z=10)
plot_dfl(dfl2, phase=1, fig_name='radiation after reflection from mirror2 and propagation')


#%%

dfl3 = generate_gaussian_dfl(1e-9, (1000, 1000, 1))
hprofile3 = dfl_reflect_surface(dfl3,                   # ocelot.optics.wave.RadiationField, which will be reflected from imperfect mirror surface (hprofile2)
                    angle=np.pi * 2 / 180,              # [radians] angle of incidence with respect to the surface
                    hrms=1e-9,                          # [m] height errors root mean square
                    axis='x',                           # direction along which reflection takes place
                    return_height_profile=True)         # boolean type variable; if it equals True the function will return internally generated height_profile

plot_1d_hprofile(hprofile3, fig_name='internally generated mirror3 height profile and PSD')
dfl3.prop(z=10)
plot_dfl(dfl3, fig_name='radiation after reflection internally generated mirror3 and propagation')

# saving HeightProfile to a file and then loading it back
# hprofile.save('imperfect_mirror.txt')
# hprofile_loaded = HeightProfile()
# hprofile_loaded.load('imperfect_mirror.txt')
# plot_1d_hprofile(hprofile_loaded, fig_name='mirror1 height profile and PSD loaded from file')
