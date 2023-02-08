# -*- coding: utf-8 -*-
"""
Created on Fri May 28 2021
@author: sserkez
"""

import numpy as np
from ocelot.utils.xfel_utils import *
from ocelot.gui.genesis_plot import *
from copy import deepcopy

fntsz = 6
params = {'image.cmap': def_cmap, 'backend': 'ps', 'axes.labelsize': 3 * fntsz, 'font.size': 3 * fntsz, 'legend.fontsize': 4 * fntsz, 'xtick.labelsize': 4 * fntsz,  'ytick.labelsize': 4 * fntsz, 'text.usetex': False}
rcParams.update(params)

#%% two pulses

# E_photon = 500 #central photon energy
# E_photon_range = 40 #photon energy range

# npoints = 700 #points in time domain (more points = wider frequency scale)
# td_scale = np.linspace(0,9e-6,npoints)
# td_env = np.zeros_like(td_scale)
# idx1 = np.logical_and(td_scale>0.5e-6,td_scale<2.5e-6)
# idx2 = np.logical_and(td_scale>5.5e-6,td_scale<7.5e-6)
# idx = np.logical_or(idx1,idx2)
# td_env[idx] = 1

# fd_scale = np.linspace(E_photon-E_photon_range/2, E_photon+E_photon_range/2, npoints) #frequency domain scale (here is it overruled by fit_scale = 'td', so won't change the scale of the output)
# fd_env = np.exp(-(fd_scale-499)**2 / 1.5**2)

# td_phase2 = np.linspace(0,800,npoints)
# td_phase1 = np.linspace(0,-800,npoints)
# td_phase = td_phase2
# td_phase[0:250] = td_phase1[0:250]

# n_events=10
#%% Sine

# E_photon = 500 #central photon energy
# E_photon_range = 40 #photon energy range
# npoints = 700 #points in time domain (more points = wider frequency scale)

# td_scale = np.linspace(0,16e-6,npoints)
# td_env = np.ones_like(td_scale)

# fd_scale = np.linspace(E_photon-E_photon_range/2, E_photon+E_photon_range/2, npoints) #frequency domain scale (here is it overruled by fit_scale = 'td', so won't change the scale of the output)
# fd_env = np.exp(-(fd_scale-499)**2 / 1.5**2)

# td_phase2 = np.linspace(0,800,npoints)


# td_phase = np.sin(td_phase2 / 63*1.5) * 63/1.5
# n_events=10
#%% Chirped pulse

E_photon = 500 #central photon energy
E_photon_range = 40 #photon energy range

npoints = 300 #points in time domain (more points = wider frequency scale)
chirp_val_lin = 100 #abs lin chirp value, arb.units
# chirp_val_lin = 0 #abs chirp value, arb.units

chirp_val_quad = -100 #abs quad chirp value, arb.units
# chirp_val_quad = 0 #abs quad chirp value, arb.units

n_events=10 #for statistics, final wigner is averaged

td_scale = np.linspace(0,6e-6,npoints) #time domain scale in um
td_env = np.ones_like(td_scale) #time domain power envelope, flat-top
td_env[td_scale<1e-6] = 0
td_env[td_scale>5e-6] = 0

fd_scale = np.linspace(E_photon-E_photon_range/2, E_photon+E_photon_range/2, npoints) #frequency domain scale (here is it overruled by fit_scale = 'td', so won't change the scale of the output)
fd_env = np.exp(-(fd_scale-E_photon)**2 / 1.5**2) #instantaneous bandwidth, constant over the pulse

td_phase = np.sign(chirp_val_lin) * np.linspace(-np.abs(chirp_val_lin)**(1/2), np.abs(chirp_val_lin)**(1/2), npoints)**2 #Phase chirp for linear frequency chirp squared
td_phase += np.sign(chirp_val_quad) * np.linspace(0, np.abs(chirp_val_quad)**(1/4), npoints)**(2+2) #Phase chirp for quadratic frequency chirp squared

#%%
(td_scale, td, fd_scale, fd) = imitate_1d_sase_like(td_scale, td_env, fd_scale, fd_env, 
                                                    #td_phase = td_phase, #does not work in fit_scale = 'td'
                                                    fd_phase = None, phen0 = None, en_pulse = None, 
                                                    fit_scale = 'td', 
                                                    n_events = n_events,
                                                    seed = 1 #This value determines random number generation for stochastic shot noise generation. Change to obtain another ensemble. 
                                                    )

#td is time domain field, it is a matrix of shape (npoints, nevents), np.abs(td) and np.angle(td) give amplitude and phase
#td_scale is time domain scale in um
#fd* is the same, but in the frequency domain

td = td * np.exp(1j * td_phase[:, np.newaxis]) #adding frequency chirp
## Now td is the complex field in time domain with introduced chirp. Now we ignore frequency domain, as we didn't add the chirp there

#%% the code below is just for plotting

plt.close(100)
plt.figure(100) #plotting power and phase of 1st event
plt.clf()
plt.plot(td_scale * 1e6, np.abs(td[:,0])**2)
plt.plot(td_scale * 1e6, np.angle(td[:,0]))
plt.legend(['intensity [arb.units]', 'phase [rad]'])
plt.xlabel('s [um]')
plt.show()


#%% Wigner averaging

for idx, tdi in enumerate(td.T):
    wig = WignerDistribution()
    wig.field = tdi
    wig.s = td_scale
    wig.xlamds = (speed_of_light * h_eV_s) / np.mean(fd_scale)
    wig.z = 0
    wig.eval()
    if idx == 0:
        W = wig.wig
    else:
        W += wig.wig

W /= n_events #normalization
#%% Wigner plotting
wig_av = deepcopy(wig)
wig_av.wig = W

plot_wigner(wig_av, figsize=3, x_units='fs', plot_text=0, plot_proj=1, plot_moments=0, fig_name='wigner_averaged', autoscale=0)
