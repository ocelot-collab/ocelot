#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 15:15:01 2023

@author: trebushi
"""

from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import *
from ocelot import ocelog
from ocelot.common.globals import *  # import of constants like "h_eV_s" and

import math
import logging
import matplotlib
import matplotlib.gridspec as gridspec

from copy import deepcopy

#%%
# create RadiationField object with Gaussian beam 
# can be adjusted to be plane wave just increase Gaussian beam waist

THz = 3
xlamds = speed_of_light/THz/1e12

b = 0.3
a = 0.055

sigma_g = 333.5
dfl = generate_gaussian_dfl(xlamds=xlamds, shape=(351, 351, 1), dgrid=(4*a, 4*a,  10e-6), power_rms=(sigma_g * 1e-2/2, sigma_g*1e-2/2,  10e-6),
                          power_center=(0, 0, None), power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                          zsep=None, freq_chirp=0, en_pulse=None, power=1e6)


dfl.fld = dfl.fld/np.sqrt(dfl.E())
    
plot_dfl(dfl, fig_name='dfl at source domain', domains = "fs")
# plot_dfl(dfl, fig_name='dfl at source k domain', domains = "fk")
#%%
# propagate through two cells of an iris line
dfl_iris_exit = deepcopy(dfl)

Nf = a**2/dfl.xlamds/b
print('Nf =', Nf)

N = 2
n_iter = 100
dfl_iris_exit, E_x, E_y, rad_left_array, loss_per_cell_array = dfl_prop_iris(dfl, N=N, a=a, center=([25e-3, 0e-3], [-20e-3, 0e-3]), 
                                                                     b=b, n_iter=n_iter, absorption_outer_pipe=0,
                                                                     acount_first_cell_loss=1)

I_x_array = np.real(E_x)**2 + np.imag(E_x)**2

#%%
# plot
r = dfl_iris_exit.scale_y()*1e3
Z = np.linspace(0, N*b, n_iter)

# Set the font size for axes labels and tick labels
fontsize = 18
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize

# Set up the figure and axes using gridspec
fig = plt.figure(figsize=(12, 6))
gs = gridspec.GridSpec(1, 3)  # Divide the figure into a 1x3 grid

# Allocate space for the main plot and subplot
ax1 = fig.add_subplot(gs[0, :2])  # Main plot occupies the first two columns
ax2 = fig.add_subplot(gs[0, 2])  # Subplot occupies the third column

# Plot the main plot
ax1.pcolormesh(Z, r, np.log(I_x_array), cmap='Greys', vmax=np.max(np.log(I_x_array))*0.5, vmin=-30)
# ax1.pcolormesh(Z, r, I_x_array1/np.max(I_x_array1, axis=0), cmap='Greys')

ax1.set_xlabel('z, m', fontsize=fontsize)
ax1.set_ylabel('x, mm', fontsize=fontsize)
ax1.set_ylim(np.min(r), np.max(r))

# Plot the subplot
ax2.plot(I_x_array[:, -1]/np.max(I_x_array), r)

ax2.set_xlabel('arb. units', fontsize=fontsize)
ax2.set_ylabel('', fontsize=fontsize)
ax2.set_ylim(np.min(r), np.max(r))

# Show the plot
plt.tight_layout()
plt.show()
