#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 13:21:22 2021

@author: andrei
"""

import numpy as np
import matplotlib.pyplot as plt
from analytics_to_OCELOT import *
from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
import scipy.special as sc
from scipy import signal
from scipy import misc
import copy

sim_dir = r'd:\DESYcloud\projects\2020_Partially_coherent\sim_dir'

ocelog.setLevel(logging.INFO)

n_s = 500
l_w = 0.018 # [m] undulator period 
L_w = l_w * n_s

E_ph = 1240 # eV
w = E_ph / hr_eV_s 
xlamds = 2 * np.pi * speed_of_light / w

sigma_r = np.sqrt(2*xlamds*L_w)/4/np.pi #natural radiation size in the waist
sigma_rp = np.sqrt(xlamds/2/L_w) #natural radiation divergence at the waist


#### #1
# ebeam_sigma_x = 1.5e-05
# ebeam_sigma_y = 5e-06
# ebeam_sigma_xp = 5e-07
# ebeam_sigma_yp = 7e-06
#### #2
# ebeam_sigma_x = 0.0001
# ebeam_sigma_y = 2e-05
# ebeam_sigma_xp = 2.5e-06
# ebeam_sigma_yp = 2.5e-05
#### #3

# k_support_SERVAL = 'intensity'
k_support_SERVAL = 'amplitude'
s_support_SERVAL='conv'
# s_support_SERVAL='beam'

ebeam_sigma_x = 0.1e-06
ebeam_sigma_y = 0.1e-06
ebeam_sigma_xp = 0.1e-06
ebeam_sigma_yp = 0.1e-06


# ebeam_sigma_x = 50e-06
# ebeam_sigma_y = 10e-06
# ebeam_sigma_xp = 0.1e-06
# ebeam_sigma_yp = 0.1e-06

ebeam_sigma_x = 0.1e-06
ebeam_sigma_y = 0.1e-06
ebeam_sigma_xp = 10e-06
ebeam_sigma_yp = 20e-06

# ebeam_sigma_x = 50e-06
# ebeam_sigma_y = 10e-06
# ebeam_sigma_xp = 10e-06
# ebeam_sigma_yp = 20e-06


ebeam_sigma_z = 2000e-6
ebeam_sigma_gamma = 1e-4 #TODO: relative electron energy spread

N_b = 100 #number of statistical realizations
N_e = 100 #number of macro electrons 
Nz, Ny, Nx = N_b, 101, 101 # the shape of the dfl.fld

str_simulation_param = 'ebeam_sigma_x = {}\n'.format(ebeam_sigma_x) + \
                       'ebeam_sigma_y = {}\n'.format(ebeam_sigma_y) + \
                       'ebeam_sigma_xp = {}\n'.format(ebeam_sigma_xp) + \
                       'ebeam_sigma_yp = {}\n'.format(ebeam_sigma_yp) + \
                       'N_b = {}\n'.format(N_b) + \
                       'N_e = {}\n'.format(N_e) + \
                       'grid mesh x = {}\n'.format(Nx) + 'grid mesh y = {}\n'.format(Ny) 

script_name = os.path.basename(__file__)
simulation_name = "{:.2E}".format(ebeam_sigma_x) + '_um_' + \
                  "{:.2E}".format(ebeam_sigma_y) + '_um_' + \
                  "{:.2E}".format(ebeam_sigma_xp) + '_urad_' + \
                  "{:.2E}".format(ebeam_sigma_yp) + '_urad_' + str(script_name.split('.')[0])

e_beam_param = r'$N_x$ = {}, '.format(round((ebeam_sigma_x)**2/xlamds/L_w, 3)) + r'$N_y$ = {}, '.format(round((ebeam_sigma_y)**2/xlamds/L_w, 3)) + \
               r'$D_x$ = {}, '.format(round((ebeam_sigma_xp)**2 * L_w/xlamds, 3)) + r'$D_y$ = {}, '.format(round((ebeam_sigma_yp)**2 * L_w/xlamds, 3)) + \
               r'$N_b$ = {} '.format(N_b) + r'$N_e = {}$'.format(N_e)
print(e_beam_param)

#%% 
### make a directory on your machine        
###saving simulation parameters in a .txt file
filePath = sim_dir + simulation_name + '/'
os.makedirs(filePath, exist_ok=True)
f = open(filePath + 'prm.txt', "w")
f.write(str_simulation_param)
f.close()

script_dir = os.getcwd() + '/' + script_name
new_script_dir = filePath + script_name
### seed for comparing fields
seed = 1234
###
#%% Monte Calro
Lz, Ly, Lx = 1000e-6, 800e-6, 800e-6 #size of realspace grid [m]
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz

### creating RadiationField object
dfl_MP = RadiationField((Nz, Ny, Nx))
dfl_MP.dx, dfl_MP.dy, dfl_MP.dz = dx, dy, dz
dfl_MP.xlamds = xlamds
dfl_MP.filePath = filePath
dfl_MP.to_domain('sf')

fieldname_MC = ''
# approximation = "far_field"
approximation = "near_field"

dfl_MP = undulator_field_dfl_MP(dfl_MP, z=5, L_w=L_w, E_ph=E_ph, N_e=N_e, N_b=N_b,
                                            sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, sig_xp=ebeam_sigma_xp, sig_yp=ebeam_sigma_yp, 
                                            approximation=approximation, mode='incoh', seed=seed)
dfl_MP.prop(-5)
plot_dfl(dfl_MP, column_3d=1, phase=1, domains='sf', fig_name='MP')
plot_dfl(dfl_MP, column_3d=1, phase=1, domains='kf', fig_name='MP')

#%% SERVAL


dfl_SERVAL = RadiationField((Nz, Ny, Nx))  
dfl_SERVAL.dx, dfl_SERVAL.dy, dfl_SERVAL.dz = dx, dy, dz
dfl_SERVAL.xlamds = xlamds # SVEA carrieer frequency
dfl_SERVAL.domain_z = 'f'

dfl_SERVAL.fld = np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) + 1j * np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) # Gaussian noise

dfl_SERVAL.filePath = filePath+'.dfl'

ocelog.warning('starting')

fieldname_SERVAL = ''
dfl_SERVAL = undulator_field_dfl_SERVAL(dfl_SERVAL, L_w=L_w, 
                                        sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, 
                                        sig_xp=ebeam_sigma_xp, sig_yp=ebeam_sigma_yp,
                                        k_support = k_support_SERVAL,
                                        s_support = s_support_SERVAL,
                                        showfig=0)

plot_dfl(dfl_SERVAL, column_3d=1, phase=1, domains='sf', fig_name='SERVAL')
plot_dfl(dfl_SERVAL, column_3d=1, phase=1, domains='kf', fig_name='SERVAL')


#%%

dfl_SERVAL.to_domain(domains='sf') 
dfl_MP.to_domain(domains='sf') 

dfl_VL_corr_s = dfl_xy_corr(dfl_SERVAL)
dfl_MP_corr_s = dfl_xy_corr(dfl_MP)

plot_two_dfls(dfl_SERVAL, dfl_MP, domains='s', label_first='SERVAL', label_second='Multiple Particles', title='Field_S', fig_name='Field_S')

# plot_dfl(dfl_VL_corr_s, phase=1, domains='sf', fig_name='SERVAL_correlation')
# plot_dfl(dfl_MP_corr_s, phase=1, domains='sf', fig_name='MP_coherence')

plot_two_dfls(dfl_VL_corr_s, dfl_MP_corr_s, domains='s', label_first='SERVAL', label_second='Multiple Particles', title='Correlation_S', fig_name='Correlation_S')

dfl_SERVAL.to_domain(domains='kf') 
dfl_MP.to_domain(domains='kf') 

plot_two_dfls(dfl_SERVAL, dfl_MP, domains='k', label_first='SERVAL', label_second='Multiple Particles', title='Field_K', fig_name='Field_K')

dfl_VL_corr_k = dfl_xy_corr(dfl_SERVAL)
dfl_MP_corr_k = dfl_xy_corr(dfl_MP)

# plot_dfl(dfl_VL_corr_k, phase=1, domains='kf', fig_name='SERVAL_correlation')
# plot_dfl(dfl_MP_corr_k, phase=1, domains='kf', fig_name='MP_correlation')

plot_two_dfls(dfl_VL_corr_k, dfl_MP_corr_k, domains='k', label_first='SERVAL', label_second='Multiple Particles', title='Correlation_K', fig_name='Correlation_K')



#%% correlation calc

# dfl_MP.to_domain(domains='sf') 
# dfl_corr = dfl_xy_corr(dfl_MP)


# plot_dfl(dfl_MP, column_3d=1, phase=1, domains='sf', fig_name='s_field')
# plot_dfl(dfl_corr, phase=1, domains='sf', fig_name='s_coherence')
# plot_dfl(dfl_corr, phase=1, domains='kf',fig_name='s_coherence_shape?')

# dfl_MP.to_domain(domains='kf')
# dfl_corr = dfl_xy_corr(dfl_MP)

# plot_dfl(dfl_MP, column_3d=1, phase=1, domains='kf', fig_name='f_field')
# plot_dfl(dfl_corr, phase=1, domains='kf',fig_name='f_coherence')
# plot_dfl(dfl_corr, phase=1, domains='sf', fig_name='f_coherence_shape?')





#%%
# approximation = 'far_field'
# dfl_MP_ff = undulator_field_dfl_MP(dfl_MP, z=5, L_w=L_w, E_ph=E_ph, N_e=N_e, N_b=N_b,
#                                             sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, sig_xp=ebeam_sigma_xp, sig_yp=ebeam_sigma_yp, 
#                                             approximation=approximation, mode='incoh', seed=seed)
# dfl_MP_ff.to_domain('k')
#%%
# # dfl_MP.to_domain(domains='sf') 
# dfl_MP.to_domain(domains='kf')





# dflc = dfl_interp(dfl_MP, 1, (0.4,0.4), return_result=1)
# # dflc = dfl_interp(dfl_MP, (0.5,0.3), 1, return_result=1)


# # dflc.fld = dflc.fld[:,:,dflc.Nx() // 2, np.newaxis]

# MCF = dflc.mut_coh_func(norm=0)


# Nxh = dflc.Nx() // 2
# Nyh = dflc.Ny() // 2

# dfl_coh = deepcopy(dflc)
# dfl_coh.fld = MCF[np.newaxis,Nyh, Nxh]

# plot_dfl(dfl_coh, fig_name='MCF_2d')


# plt.figure('MCF')
# plt.close()
# plt.figure('MCF')



# plt.plot(dflc.scale_y(), np.abs(MCF[Nyh,Nxh,:  ,Nxh]), label='y abs', linestyle='-', color='b')
# plt.plot(dflc.scale_y(), np.angle(MCF[Nyh,Nxh,:  ,Nxh]), label='y angle', linestyle='--', color='b')
# plt.plot(dflc.scale_x(), np.abs(MCF[Nyh,Nxh,Nyh,:  ]), label='x abs', linestyle='-', color='r')
# plt.plot(dflc.scale_x(), np.angle(MCF[Nyh,Nxh,Nyh,:  ]), label='x angle', linestyle='--', color='r')
# plt.legend()
# plt.show()



# dfl_corr = RadiationField()

# #%%


    
    


# dfl_MP.to_domain(domains='sf') 
# # dfl_MP.to_domain(domains='kf')

# dfl_tmp = dfl_interp(dfl_MP, 1, 1, return_result=1)

# # dfl1 = dfl_tmp.fld[:, :, :, np.newaxis, np.newaxis].conjugate()
# # dfl2 = dfl_tmp.fld[:, np.newaxis, np.newaxis, :, :]

# Nxh = dfl_tmp.Nx() // 2
# Nyh = dfl_tmp.Ny() // 2

# dfl1 = dfl_tmp.fld[:, Nyh, Nxh, np.newaxis, np.newaxis].conjugate()
# dfl2 = dfl_tmp.fld[:, :, :]

# print(dfl1.shape)
# print(dfl2.shape)

# J = np.mean(dfl1 * dfl2,
#                 axis=0)

# print(J.shape)


# dfl_coh = deepcopy(dfl_tmp)
# dfl_coh.fld = J[np.newaxis,:,:]

# plot_dfl(dfl_coh, fig_name='J_orig')

# I = dfl_tmp.int_xy()
# J = J / (I[Nyh, np.newaxis :] * I[:, Nxh, np.newaxis])

# dfl_coh.fld = J[np.newaxis,:,:]
# plot_dfl(dfl_coh, fig_name='J_norm')
#%%

# filePath = sim_dir + simulation_name + '/'
# write_field_file(dfl_MP, filePath=filePath, fileName=fieldname_MC)

# #%%
# # Define mesh size
# Lz, Ly, Lx = 10000e-6, 1200e-6, 1000e-6 #size of realspace grid [m]
# dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz

# # create radiation field
# dfl_SERVAL = RadiationField((50, Ny, Nx))  
# dfl_SERVAL.dx, dfl_SERVAL.dy, dfl_SERVAL.dz = dx, dy, dz
# dfl_SERVAL.xlamds = xlamds # SVEA carrieer frequency

# dfl_SERVAL.fld = np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) + 1j * np.random.randn(dfl_SERVAL.Nz(), dfl_SERVAL.Ny(), dfl_SERVAL.Nx()) # Gaussian noise

# dfl_SERVAL.filePath = filePath+'.dfl'
# # dfl_omega_0 = 2*np.pi * speed_of_light / dfl.xlamds # set undulator resonance to SVEA carrier frequency (to the middle of the photon energy domain mesh)
# # radiation_omega_resonance = dfl_omega_0 * 1.01 # shift undulator resonance to the right

# fieldname_SERVAL = '1-far_zone_50_m_SERVAL'
# dfl_SERVAL = SERVAL_undulator_field_generator(dfl_SERVAL, L_w=L_w, 
#                                         sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, 
#                                         sig_xp=ebeam_sigma_xp, sig_yp=ebeam_sigma_yp,
#                                         showfig=True)

# #%%
# dfl_prop_SERVAL = deepcopy(dfl_SERVAL)

# dfl_prop_SERVAL.prop_m(50, m=[15, 10])
# dfl_prop_SERVAL.to_domain('st')

# dfl_prop_SERVAL.to_domain(domains='sf') 

# plot_dfl(dfl_prop_SERVAL, domains='sf', phase=True, fig_name = fieldname_SERVAL)
# # plot_dfl(dfl_prop_SERVAL, domains='kf', phase=True, fig_name = fieldname_SERVAL)

# filePath = sim_dir + simulation_name + '/'
# write_field_file(dfl_prop_SERVAL, filePath=filePath, fileName=fieldname_SERVAL)

# #%%

# fieldname_MC = '1-far_zone_50_m_MC'
# fieldname_SERVAL = '1-far_zone_50_m_SERVAL'

# filepath_MC = filePath + fieldname_MC
# filepath_SERVAL = filePath + fieldname_SERVAL

# dfl_MP = read_field_file(filepath_MC)
# dfl_SERVAL = read_field_file(filepath_SERVAL)

# simulation_name = simulation_name.replace('.', '_')

# #%%
# plot_two_dfls(dfl_MP, dfl_SERVAL, domains='s', figname='1-far_zone_50_m' + simulation_name, 
#               slice_xy=False, phase=False, label_first='MC', 
#               label_second='SERVAL', title=e_beam_param, filePath=filePath, savefig=True)


# #%%
# dfl2_MC = deepcopy(dfl_MP)
# dfl2_SERVAL = deepcopy(dfl_SERVAL)

# dfl2_MC.to_domain(domains='sf')
# dfl2_SERVAL.to_domain(domains='sf')

# ap_x = 2.5e-3
# ap_y = 2.5e-3
# dfl_ap_rect(dfl2_MC, ap_x=ap_x, ap_y=ap_y)
# dfl_ap_rect(dfl2_SERVAL, ap_x=ap_x, ap_y=ap_y)

# # interpL = 0.25
# # interpN = 4
# # dfl_interp(dfl2_MC, interpN=(interpN, interpN), interpL=(interpL, interpL), method='quintic')
# # dfl_interp(dfl2_SERVAL, interpN=(interpN, interpN), interpL=(interpL, interpL), method='quintic')

# dfl2_MC.prop(10)
# dfl2_SERVAL.prop(10)
# #%%
# fieldname = '2-far_zone_60_m_after_aperture'
# # plot_dfl(dfl2_SRW, domains='sf', phase=True, fig_name = fieldname)
# # plot_dfl(dfl2_SRW, domains='kf', phase=True, fig_name = fieldname)

# # plot_dfl(dfl2_MC, domains='s', phase=True, fig_name = filePath + fieldname)


# plot_two_dfls(dfl2_MC, dfl2_SERVAL, domains='s', figname=fieldname + simulation_name, 
#               slice_xy=True, phase=False, label_first='MC', 
#               label_second='SERVAL', title=e_beam_param, filePath=filePath, savefig=True)

# #%%
# dfl3_MC = deepcopy(dfl2_MC)
# dfl3_SERVAL = deepcopy(dfl2_SERVAL)

# f = 60*10/(60+10)
# dfl3_MC.curve_wavefront(r=f)
# dfl3_SERVAL.curve_wavefront(r=f)

# m_x = 0.03
# m_y = 0.08
# dfl3_MC.prop_m(10, m=[m_x, m_y])
# dfl3_SERVAL.prop_m(10, m=[m_x, m_y])

# #%%
# fieldname = '3-70_m_focal_plane'
# plot_two_dfls(dfl3_MC, dfl3_SERVAL, domains='s', figname=fieldname  + simulation_name, 
#               slice_xy=True, phase=False, label_first='MC', 
#               label_second='SERVAL', title=e_beam_param, filePath=filePath, savefig=True)


# # save of your python script in the simulation directory
# write_script(script_dir, new_script_dir)



















