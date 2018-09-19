'''
user interface for viewing Stokes parameters
'''
import matplotlib.pyplot as plt
import numpy as np
from ocelot.optics.wave import *
from ocelot.gui.genesis_plot import plot_stokes_3d, plot_dfl
ocelog.setLevel(logging.DEBUG)
    
#%% Generation of Field and Stokes
L = 7
l = 7
grid_tr = 15e-4
shape_xy = 251*1
shape_z = 51
rms_xy = 19e-6
E_ph = 300
dE_ph = 0.03
lamds =  h_eV_s * speed_of_light / E_ph
dlamds = h_eV_s * speed_of_light / (E_ph + dE_ph) - lamds

dflR = generate_dfl(xlamds=lamds, 
                    wavelength = lamds, 
                    dgrid=(grid_tr, grid_tr, 50e-6), 
                    shape=(shape_xy, shape_xy, shape_z), 
                    power_rms=(rms_xy, rms_xy, 5e-6), 
                    power_center = (300e-6, -150e-6, 0.25e-4), 
                    power = (1e6), 
                    phase_chirp_lin = 0)
dflL = generate_dfl(xlamds=lamds, 
                    wavelength = lamds + dlamds, 
                    dgrid=(grid_tr, grid_tr, 50e-6), 
                    shape=(shape_xy, shape_xy, shape_z), 
                    power_rms=(rms_xy,  rms_xy,  5e-6), 
                    power_center = (300e-6, -150e-6, 0.25e-4), 
                    power = (1e6), 
                    phase_chirp_lin = 0)

dflR.prop_m( 5, m=1)
dflL.prop_m(-10, m=1)
plot_dfl(dflR, cmap='Greys')

S = calc_stokes_dfl(dflR, dflL, basis = 'rl')

plot_stokes_3d(S, x_plane='max_slice', y_plane='max_slice', z_plane='max_slice', fig_name='slice_properties', cbars=True)
plot_stokes_3d(S, x_plane='proj', y_plane='proj', z_plane='proj', fig_name='projected_properties', cbars=False)