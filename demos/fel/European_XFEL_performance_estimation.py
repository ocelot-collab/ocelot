# -*- coding: utf-8 -*-
"""
Created on Fri Jul 6 15:00:00 2018

@author: sserkez
"""

import ocelot
from ocelot.rad.fel import *
from ocelot import ocelog
from ocelot.common.globals import *
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import matplotlib.colors as colors
from ocelot.common.math_op import find_nearest_idx
from ocelot.common.ocelog import ocelog
from ocelot.rad.undulator_params import *

ocelog.setLevel(logging.ERROR)

#%% Calculation
cmap = 'terrain'
class tmp():
        pass

# E_beam = 7.8 # GeV # CW operation energy
# E_beam = 17.5 # GeV # maximum
# E_beam = 14 # GeV # standard
# E_photon = 10000 # photon energy

dE_beam = 3e-3 # GeV # energy spread
emit_beam = 0.4 # normalized emittance [mm mrad]
beta_beam = 30 # average beta function [m]
I_beam = 5000 # peak current [A]
t_pulse = 1e-15 # flattop e-beam "duration" [s]
L_magn = 35*5 # magnetic undulator length [m] (SASE1/2)
und_duty_factor = 5/(5+1.2) # undulator duty factor
xlamd = 0.04 # undulator period [m] (SASE1/2)

nharm = 3

method = 'mxie' # Estimation method: {M.Xie} or {Saldin Shneidmiller Yurkov}: 'ssy_opt'

E_beam_arr = np.geomspace(7.5, 17.5, 25)
E_photon_arr = np.geomspace(3000, 40000, 55) # (SASE1/2)

## uncomment for _SASE3_ parameters
# L_magn = 35*5 # m (SASE3)
# xlamd = 0.068 # m (SASE3)
# E_photon_arr = np.geomspace(250, 3000, 15) # (SASE3)

x_arr = E_photon_arr / 1e3
y_arr = E_beam_arr
x_label = '$E_{photon}$ [keV]'
y_label = '$E_{beam}$ [GeV]'
xlim = [np.amin(x_arr), np.amax(x_arr)]
ylim = [np.amin(y_arr), np.amax(y_arr)]

arr = np.zeros((len(y_arr), len(x_arr)))
K_arr = deepcopy(arr)
P_Lmagn_arr = deepcopy(arr)
P_sat = deepcopy(arr)
sat_power_arr = deepcopy(arr)
beta_arr = deepcopy(arr)
coherence = deepcopy(arr)
medium_accuracy_arr = deepcopy(arr)
phen_arr = deepcopy(arr)
beam_power_arr = deepcopy(arr)
rho_arr = deepcopy(arr)
lg3d_arr = deepcopy(arr)
z_sat_magn_arr = deepcopy(arr)
aw_arr = deepcopy(arr)
beta_arr = deepcopy(arr)
N_phot_arr = deepcopy(arr)


# for xi in [0]:
#     for yi in [0]:
for xi, E_photon in enumerate(E_photon_arr):
    for yi, E_beam in enumerate(E_beam_arr):
        
        xlamds = eV2lambda(E_photon)
        
        gamma = E_beam / m_e_GeV
        delgam = dE_beam / m_e_GeV
        aw = np.sqrt(2 * gamma**2 * xlamds / xlamd - 1)
        K_arr[yi, xi] = aw * np.sqrt(2)
        
        tmp.aw0 = aw        
        tmp.iwityp = 0 # 0 only!!! 1 = helical, 0 = planar
        
        tmp.curpeak = I_beam
        tmp.gamma0 = gamma
        tmp.delgam = delgam
        tmp.xlamd = xlamd
        
        tmp.emitx = emit_beam * 1e-6
        tmp.emity = emit_beam * 1e-6
        tmp.betax = beta_beam
        tmp.betay = beta_beam
        tmp.Lg_mult = 1 # empirical lain length multiplication parameter
        tmp.P_mult = 0.5 # empirical saturation power multiplication parameter
        
        tmp.hn = nharm # harmonic number for harmonic lasing
        tmp.qf = 1 # quantun fluctuation effects
        
        fel = calculateFelParameters(tmp, method=method)
        
        # if fel.delta_q >= 1:
        #     print(fel.delta_q)
        #     print(fel.phenh)
        # fel.beta_opt(apply=1, method=fel.method) # calculate optimal beta function (time consuming if method is "mxie", "ssy_opt" automatically assumes optimal beta)
        
        phen_arr[yi, xi] = fel.phenh
        beam_power_arr[yi, xi] = fel.Pb
        rho_arr[yi, xi] = fel.rho3
        lg3d_arr[yi, xi] = fel.lg3
        
        en = 2 * np.pi * (fel.emitx+fel.emity) / 2 / fel.lambda0 / fel.gamma0
        coherence[yi, xi] = 1.1 * en**(1/4) / (1+0.15*en**(9/4))
        # coherence[yi, xi] = (np.log(fel.Nc/en)/4/en)**2
        
        z_sat_magn_arr[yi, xi] = fel.z_sat_magn
        P_Lmagn_arr[yi, xi] = fel.P(L_magn)
        P_sat[yi, xi] = fel.P()
        sat_power_arr[yi, xi] = fel.P(fel.z_sat_magn)
        beta_arr[yi, xi] = np.mean([fel.betax, fel.betay])
        medium_accuracy_arr[yi,xi] = fel.inaccurate

x_arr *= nharm
xlim[0] *= nharm
xlim[1] *= nharm

P_arr = P_Lmagn_arr # power at saturation limited by undulator magnetic length
# P_arr = P_sat #power at saturation if unlimited undulator is available

P_arr[np.isnan(P_arr)] = 0
E_arr = P_arr * t_pulse * 0.5
N_phot =  E_arr / q_e / phen_arr


# printFelParameters(fel)   
# print('z = {:.2f}, P={:.2f}'.format(fel.z_sat_magn / 0.8, fel.P() / 1e9))
# plt.figure(5000)
# plt.scatter(fel.z_sat_magn / 0.8, fel.P() / 1e9)
# plt.xlim(xmin=0)
# plt.show()

#%% Plotting

## Available range of K parameter
K_avail = np.logical_and(K_arr > 1.65, K_arr < 3.9) #current SASE12
# K_avail = np.logical_and(K_arr > 1.3, K_arr < 3.9) #extended SASE12
# K_avail = np.logical_and(K_arr > 4, K_arr < 9) #current SASE3
# K_avail = np.logical_and(K_arr > 1.0, K_arr < 2.5) #SCU18
# K_avail = np.logical_and(K_arr > 1.0, K_arr < 2.7) #SCU20
# K_avail = np.logical_and(K_arr > 1.0, K_arr < 10) #SCU40
# K_avail = np.ones_like(K_avail)


plt.figure('Z_sat')
plt.clf()
plt.title('Saturation length (physical) [m]')
plt.pcolormesh(x_arr, y_arr, (z_sat_magn_arr / und_duty_factor) * K_avail, cmap=cmap, vmin=0)#, vmax=L_magn / und_duty_factor * 2)
plt.colorbar()

if np.any(z_sat_magn_arr / L_magn < 1) and np.any(z_sat_magn_arr / L_magn > 1):
    cs = plt.contour(x_arr, y_arr, z_sat_magn_arr / L_magn, [1], colors='r')
    plt.clabel(cs, inline=1, fontsize=10, fmt=' undulator {} m'.format(L_magn/und_duty_factor))
    
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel(x_label)
plt.ylabel(y_label)


plt.figure('Energy')
plt.clf()
plt.title('Power at saturation [uJ/fs]')
plt.pcolormesh(x_arr, y_arr, E_arr*1e6 * K_avail, cmap=cmap)#, norm=colors.LogNorm())#vmin=1e8, vmax=1e12
plt.colorbar()
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel(x_label)
plt.ylabel(y_label)


# plt.figure('beta')
# plt.clf()
# plt.title('Beta function [m]')
# plt.pcolormesh(x_arr, y_arr, beta_arr, cmap=cmap, vmin=0, vmax=30)
# plt.colorbar()
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(x_label)
# plt.ylabel(y_label)      
        

# plt.figure('N_phot')
# plt.clf()
# plt.title('N_phot [number]')
# plt.pcolormesh(x_arr, y_arr, N_phot * t_pulse/1e-15 * K_avail, cmap=cmap, norm=colors.LogNorm())#vmin=1e8, vmax=1e12
# plt.colorbar()
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(x_label)
# plt.ylabel(y_label)


# plt.figure('Power')
# plt.clf()
# plt.title('Power at saturation [GW]')
# plt.pcolormesh(x_arr, y_arr, P_arr/1e9 * K_avail, cmap=cmap)#, norm=colors.LogNorm())#vmin=1e8, vmax=1e12
# plt.colorbar()
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(x_label)
# plt.ylabel(y_label)



# plt.figure('K_avail')
# plt.clf()
# plt.title('K available')
# plt.pcolormesh(x_arr, y_arr, K_avail, cmap=cmap)#, vmin=0, vmax=1)
# # plt.colorbar()
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(x_label)
# plt.ylabel(y_label)


# plt.figure('accurate')
# plt.clf()
# plt.title('High estimation accuracy')
# plt.pcolormesh(x_arr, y_arr, np.logical_not(medium_accuracy_arr), cmap=cmap)#, vmin=0, vmax=1)
# # plt.colorbar()
# plt.xlim(xlim)
# plt.ylim(ylim)
# plt.xlabel(x_label)
# plt.ylabel(y_label)



plt.show()