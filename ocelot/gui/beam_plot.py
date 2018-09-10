

import sys
import os
import csv
import time
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
from ocelot.adaptors.genesis import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.common.math_op import *  # import of mathematical functions
from ocelot.utils.xfel_utils import *
from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.optics.wave import *

# from pylab import rc, rcParams #tmp
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

def_cmap = 'viridis'
# def_cmap = 'Greys'

fntsz = 4
params = {'image.cmap': def_cmap, 'backend': 'ps', 'axes.labelsize': 3 * fntsz, 'font.size': 3 * fntsz, 'legend.fontsize': 4 * fntsz, 'xtick.labelsize': 4 * fntsz,  'ytick.labelsize': 4 * fntsz, 'text.usetex': False}
rcParams.update(params)

def plot_fel_power_z(fel, z=None, fig=None, magn_coeff=1, fs=False):
    '''
    plots estimated FEL power at position z
    fel is FelParamterArray object
    '''
    if fig.__class__ == int or fig is None:
        fig = plt.figure(fig)
    ax = fig.gca()
    
    if z is None:
        z = fel.z_sat_min
    
    if fs:
        ax.plot(fel.s * speed_of_light *1e15, fel.P(z))
        ax.set_xlabel('t [fs]')
    else:
        ax.plot(fel.s * 1e6, fel.P(z))
        ax.set_xlabel('s [um]')
    # z_ind = find_nearest_idx(out.z, fel.z_sat_min * magn_coeff )
    # plt.plot(out.s * 1e6, out.p_int[:,z_ind])
    ax.set_title('Pulse energy @ %.2fm = %.2e J' %(z * magn_coeff, fel.E(z)))
    ax.set_ylabel('P [W]')
    
    #fig.savefig(exp_dir + 'power_sat.png', format = 'png')
    fig.show()

def plot_fel_spectrogram(fel, z=None, fig=None, magn_coeff=1):
    
    if fig.__class__ == int or fig is None:
        fig = plt.figure(fig)
    ax = fig.gca()
    
    if z is None:
        z = fel.z_sat_min
    Psat = fel.P(z)
    Psat[np.isnan(Psat)]=0
    idx = fel.idx

    phen0 = fel.phen0
    dphen = phen0 * fel.rho3
    dp = dphen[idx] / 10
    s_arr = fel.s * 1e6
    phen_arr = np.arange(np.amin(phen0 - 3 * dphen), np.amax(phen0 + 3 * dphen), dp)
    spec = np.zeros((s_arr.size, phen_arr.size))
    for i in range(s_arr.size):
        if dphen[i] != 0:
            spec[i] = np.exp(-(phen_arr - phen0[i])**2 / 2 / dphen[i]**2) / np.sqrt(2 * np.pi * dphen[i]**2)
    spec = spec * Psat[:, np.newaxis]

    
    ax.pcolormesh(s_arr, phen_arr, spec.T)
    ax.set_ylabel('E [eV]')
    ax.set_xlabel('s [um]')
    ax.autoscale(tight=1)
    fig.show()
    
    
    
    
    
    
    