#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 17:33:00 2021

@author: andrei
"""
import numpy as np
import matplotlib.pyplot as plt
from ocelot.optics.wave import * # import OCELOT propagation module
from ocelot.gui.dfl_plot import *  # import plotting routine
from ocelot.common.globals import *  # import of constants like "h_eV_s" and

ocelog.setLevel(logging.INFO)
_logger = logging.getLogger(__name__)

n_s = 200 # number of periods in an undulator
l_w = 0.018 # [m] undulator period 
L_w = l_w * n_s # undulator length

E_ph = 2167 # resonance energy
w = E_ph / hr_eV_s 
xlamds = 2 * np.pi * speed_of_light / w

sigma_r = np.sqrt(2*xlamds*L_w)/4/np.pi #natural radiation size in the waist
sigma_rp = np.sqrt(xlamds/2/L_w) #natural radiation divergence at the waist

#### #1
ebeam_sigma_x = 38e-06 # x electron beam size 
ebeam_sigma_y = 4.68e-06 # y electron beam size
ebeam_sigma_xp = 25e-06 # x electron beam divergence
ebeam_sigma_yp = 20e-06 # y electron beam divergence

ebeam_sigma_z = 1 # Electron beam duration. Not physical for the current version of
                  # SERVAL as here is considered the case of fully monochromatic radiation
                  # long beam duration corresponds to an infinitesimally small spectral line

N_b = 150 #number of statistical realizations
Nz, Ny, Nx = N_b, 601, 601 # the shape of the dfl.fld

e_beam_param = 'N_x = {}, '.format(round((ebeam_sigma_x)**2/xlamds/L_w, 3)) + 'N_y = {}, '.format(round((ebeam_sigma_y)**2/xlamds/L_w, 3)) + \
               'D_x = {}, '.format(round((ebeam_sigma_xp)**2 * L_w/xlamds, 3)) + 'D_y = {}, '.format(round((ebeam_sigma_yp)**2 * L_w/xlamds, 3)) + \
               'N_b = {} '.format(N_b)
print(e_beam_param)

Lz, Ly, Lx = ebeam_sigma_z, 300e-6, 400e-6 #size of realspace grid [m]
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz

# create radiation field
dfl_SERVAL = RadiationField((N_b, Ny, Nx))  
dfl_SERVAL.dx, dfl_SERVAL.dy, dfl_SERVAL.dz = dx, dy, dz
dfl_SERVAL.xlamds = xlamds # SVEA carrieer frequency

dfl_SERVAL.fld = np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) + 1j * np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) # Gaussian noise

filePath = 'SERVAL_field'
dfl_SERVAL.filePath = filePath+'.dfl'

fieldname_SERVAL = '0-at_source_SERVAL'
dfl_SERVAL = undulator_field_dfl_SERVAL(dfl_SERVAL, L_w=L_w, 
                                        sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, 
                                        sig_xp=ebeam_sigma_xp, sig_yp=ebeam_sigma_yp,
                                        k_support = 'intensity', s_support='intensity', showfig=False)

plot_dfl(dfl_SERVAL, domains='sf', phase=True, fig_name = fieldname_SERVAL)
plot_dfl(dfl_SERVAL, domains='kf', phase=True, fig_name = fieldname_SERVAL)
#%%
dfl_prop_SERVAL = deepcopy(dfl_SERVAL)

dfl_prop_SERVAL.prop_m(25, m=[15, 20])
dfl_prop_SERVAL.to_domain(domains='sf') 

fieldname_SERVAL = '1-far_zone_25_m_SERVAL'
plot_dfl(dfl_prop_SERVAL, domains='sf', phase=True, fig_name = fieldname_SERVAL)
plot_dfl(dfl_prop_SERVAL, domains='kf', phase=True, fig_name = fieldname_SERVAL)
#%%
dfl2_SERVAL = deepcopy(dfl_prop_SERVAL)
dfl2_SERVAL.to_domain(domains='sf')
ap_x, ap_y = 1e-3, 1e-3
dfl_ap_rect(dfl2_SERVAL, ap_x=ap_x, ap_y=ap_y)
dfl3_SERVAL = deepcopy(dfl2_SERVAL)

f = 25*25/(25+25)
dfl3_SERVAL.curve_wavefront(r=f)
dfl3_SERVAL.prop_m(12.5, m=0.25)

fieldname_SERVAL = '3-far_zone_12_5_m_after_aperture_SERVAL'
plot_dfl(dfl3_SERVAL, domains='sf', phase=True, fig_name = fieldname_SERVAL)
plot_dfl(dfl3_SERVAL, domains='kf', phase=True, fig_name = fieldname_SERVAL)
#%%
dfl4_SERVAL = deepcopy(dfl3_SERVAL)
m_x, m_y = 0.3, 0.3
dfl4_SERVAL.prop_m(12.5, m=[m_x, m_y])

fieldname_SERVAL = '4-in_focus_SERVAL'
plot_dfl(dfl4_SERVAL, domains='sf', phase=True, fig_name = fieldname_SERVAL)
plot_dfl(dfl4_SERVAL, domains='kf', phase=True, fig_name = fieldname_SERVAL)



