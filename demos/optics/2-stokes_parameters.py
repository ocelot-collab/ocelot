'''
user interface for viewing Stokes parameters
'''
import matplotlib.pyplot as plt
import numpy as np
from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import plot_stokes_3d, plot_dfl, plot_stokes_angles, plot_stokes_values
ocelog.setLevel(logging.DEBUG)
from copy import deepcopy
    
#%% Generation of Field and Stokes

#l1 = 8; l2 = 13 #drft from source
##l1 = 8; l2 = -8
#
#grid_tr = 15e-4 #transverse grid size [m]
#shape_xy = 51*2+1 #number of mesh points
#shape_z = 51*2+1
#rms_xy = 19e-6 #ixe of source waist
#E_ph = 300 #photon energy
#dE_ph = 0.1 #photon energy difference
##dE_ph = 0
#
#lamds =  h_eV_s * speed_of_light / E_ph
#dlamds = h_eV_s * speed_of_light / (E_ph + dE_ph) - lamds
#
#dflR = generate_gaussian_dfl(xlamds=lamds, 
#                    wavelength = lamds, 
#                    dgrid=(grid_tr, grid_tr, 50e-6), 
#                    shape=(shape_xy, shape_xy, shape_z), 
#                    power_rms=(rms_xy, rms_xy, 5e-6), 
##                    power_center = (300e-6, -150e-6, 0.25e-4), 
#                    power = (1e6), 
#                    phase_chirp_lin = 0)
#
#dflL = generate_gaussian_dfl(xlamds=lamds, 
#                    wavelength = lamds + dlamds, 
#                    dgrid=(grid_tr, grid_tr, 50e-6), 
#                    shape=(shape_xy, shape_xy, shape_z), 
#                    power_rms=(rms_xy,  rms_xy,  5e-6), 
##                    power_center = (300e-6, -150e-6, 0.25e-4), 
#                    power = (1e6), 
#                    phase_chirp_lin = 0)
#
#dflR.prop_m(l1, m=1)
#dflL.prop_m(l2, m=1)
#plot_dfl(dflR, cmap='Greys')
#
#S = calc_stokes_dfl(dflR, dflL, basis = 'rl')



L = 7
l = 7
grid_xy = 15e-4
grid_z = 50e-6
shape_xy = 251*1
shape_z = 51
rms_xy = 19e-6
rms_z = 5e-6
E_ph = 300
dE_ph = 0.03
power_center = (300e-6, -150e-6, 0.25e-4)
lamds =  h_eV_s * speed_of_light / E_ph
dlamds = h_eV_s * speed_of_light / (E_ph + dE_ph) - lamds

dflR = generate_gaussian_dfl(xlamds=lamds, 
                    wavelength = lamds, 
                    dgrid=(grid_xy, grid_xy, grid_z), 
                    shape=(shape_xy, shape_xy, shape_z), 
                    power_rms=(rms_xy, rms_xy, rms_z), 
                    power_center = power_center, 
                    power = (1e6), 
                    phase_chirp_lin = 0)
dflL = generate_gaussian_dfl(xlamds=lamds, 
                    wavelength = lamds + dlamds, 
                    dgrid=(grid_xy, grid_xy, grid_z), 
                    shape=(shape_xy, shape_xy, shape_z), 
                    power_rms=(rms_xy,  rms_xy,  rms_z), 
                    power_center = power_center, 
                    power = (1e6), 
                    phase_chirp_lin = 0)

dflR.prop_m( 5, m=1)
dflL.prop_m(-10, m=1)

plot_dfl(dflR, cmap='Greys')

S = calc_stokes_dfl(dflR, dflL, basis = 'rl')

#%% fix
#plot_stokes_values(S[25, np.newaxis, 125, np.newaxis, :], direction='x', fig=plt.figure('Stokes_x_vals'))
#plot_stokes_angles(S[25, np.newaxis, 125, np.newaxis, :], direction='x', fig=plt.figure('Stokes_x_ang'))
#%% fix
#plot_stokes_values(S[:, 125, np.newaxis, 125, np.newaxis], direction='z', fig=plt.figure('Stokes_z_vals'))
#plot_stokes_angles(S[:, 125, np.newaxis, 125, np.newaxis], direction='z', fig=plt.figure('Stokes_z_ang'))
#%%
#normalization = 's0' #normalize to s0 (degree of polarization is depicted, without hint of radiation intensity)
normalization = 's0_max' #normalize to the maximum value of s0 (shall be rewritten in terms of power)
plot_stokes_3d(S, cmap_lin='brightwheel', x_plane='max_slice', y_plane='max_slice', z_plane='max_slice', fig_name='slice_properties', cbars=True, normalization=normalization)
plot_stokes_3d(S, cmap_lin='brightwheel', x_plane='proj', y_plane='proj', z_plane='proj', fig_name='projected_properties', cbars=True, normalization=normalization)