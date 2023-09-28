#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 10:49:46 2023

@author: trebushi
"""

from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import *
from ocelot import ocelog
from ocelot.common.globals import *  # import of constants like "h_eV_s" and

import math
import logging
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import special

from copy import deepcopy

def plot_iris_line_prop_result(E_x, E_y, dfl_iris_exit, a, b, N, n_iter, rad_left_array, loss_per_cell_array, 
                               direction = 'x', savefig=False, filePath='', loss_per_cell_max = 0.1):

    fontsize=18
    avr_cells = 500
    direction = direction
    
    Z = np.linspace(0, N*b, n_iter)
    fig, axs = plt.subplots(3, 2, figsize=(15,8), gridspec_kw={'width_ratios': [4, 1]})
    fig.subplots_adjust(top=0.95, bottom=0.2, left=0.09, right=0.95)  # adjust this value to fit your needs
    
    cells = np.arange(0, N-1)
    
    if n_iter/N == 1:
        rad_left_array_ext = rad_left_array
        rad_left_array_ext = np.append(rad_left_array_ext, rad_left_array[-1])
    else:
        rad_left_array_ext = np.repeat(rad_left_array, (n_iter/N))#.tolist()
    
    if direction=='x':
        I_x_array1 = np.real(E_x)**2 + np.imag(E_x)**2
        r = dfl_iris_exit.scale_x()*1e3
        I_2D = I_x_array1[:,:]/np.max(I_x_array1[:,:]) * (1/rad_left_array_ext)
    elif direction=='y':
        I_y_array1 = np.real(E_y)**2 + np.imag(E_y)**2
        r = dfl_iris_exit.scale_y()*1e3
        I_2D = I_y_array1[:,:]/np.max(I_y_array1[:,:]) * (1/rad_left_array_ext)
    
        
    axs[0, 0].pcolormesh(Z, r, np.log(I_2D), cmap='Greys')
    axs[1, 0].pcolormesh(Z, r, I_2D, cmap='Greys', vmax=np.max(I_2D)* 0.5)
    
    axs[0, 1].semilogx(I_2D[:,-1]/np.max(I_2D[:,-1]), r)
    axs[0, 1].set_xlim(0)
    axs[0, 1].grid()
    
    axs[1, 1].plot(I_2D[:,-1]/np.max(I_2D[:,-1]), r)
    
    # axs[0, 0].set_xticks([])
    axs[0, 0].set_xticklabels([])
    axs[0, 0].set_xlabel('') 
    axs[0, 0].set_ylabel('x, [mm]', fontsize=fontsize) 
    
    # axs[1, 0].set_xticks([])
    axs[1, 0].set_xticklabels([])
    axs[1, 0].set_xlabel('') 
    axs[1, 0].set_ylabel('x, [mm]', fontsize=fontsize) 
    
    axs[0, 1].set_xlabel('log(I), arb.units', fontsize=fontsize) 
    axs[1, 1].set_xlabel('I, arb.units', fontsize=fontsize) 
    
    axs[1, 1].set_xlim(0)
    axs[1, 1].grid()
    
    axs[2, 0].plot(cells, loss_per_cell_array, '--', color='black', label='Losses per cell, %')
    mean_loos = np.mean(loss_per_cell_array[-avr_cells:-1])
    axs[2, 0].axhline(y=mean_loos, color='gray', alpha=0.5, linewidth=1, linestyle='--')
    axs[2, 0].scatter(cells[0], mean_loos, color='black', s=50)   
    
    label_position = (cells[0], mean_loos)
    label = '{}%'.format(round(mean_loos, 3))
    bbox_props = dict(boxstyle="round", facecolor="white", alpha=0.95)  # Set alpha value to 0.5 for semitransparency
    axs[2, 0].annotate(label, xy=label_position, textcoords="offset points", xytext=(30,-18), ha='center', rotation=0, backgroundcolor='white', bbox=bbox_props)
    
    # Hide the existing x-axis of axs[2, 0]
    # axs[2, 0].spines['bottom'].set_visible(False)
    # axs[2, 0].xaxis.set_ticks([])
    
    # axs[2, 0].legend()
    axs[2, 0].set_xlim(0, N)
    axs[2, 0].set_ylim(0, loos_per_cell_max)
    axs[2, 0].set_ylabel('Losses per cell, %', fontsize=fontsize) 
    axs[2, 0].set_xlabel('iris number', fontsize=fontsize) 
    axs[2, 0].grid()
    
    ax2 = axs[2, 0].twinx()
    ax2.plot(cells, rad_left_array, color='red', label='Power left, %')   
    ax2.scatter(cells[-1], rad_left_array[-1], color='red', s=50)  
    ax2.yaxis.set_tick_params(labelcolor='red', color='red')
    
    label_position = (cells[-2], rad_left_array[-1])
    label = '{}%'.format(round(rad_left_array[-1], 1))
    
    bbox_props = dict(boxstyle="round", facecolor="white", alpha=0.9)  # Set alpha value to 0.5 for semitransparency
    
    ax2.annotate(label, xy=label_position, textcoords="offset points", xytext=(-20,10), ha='center', color='r', backgroundcolor='white', bbox=bbox_props)
    
    # ax2.plot(cells, total_loss_array, label='Total losses')
    # ax2.plot(np.arange(2, N)*z, L_p[1:], label='theory')
    # ax2.legend()
    ax2.set_ylabel('Power left, %', fontsize=fontsize, c='r')#, rotation=-90)
    
    # Get the handles and labels from both axes
    handles, labels = axs[2, 0].get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    
    # Combine handles and labels
    handles.extend(handles2)
    labels.extend(labels2)
    
    # Create a single legend using the combined handles and labels
    ax2.legend(handles, labels, bbox_to_anchor=(1.42,0.7))#loc='upper left')
    
    # ax2.set_ylim(0)
    ax2.set_xlim(0, np.max(cells))
    # ax2.set_yscale('log')
    
    # ax2.grid()
    fig.delaxes(axs[2, 1])
    
    # Create a new subplot for the secondary x-axis
    ax3 = fig.add_subplot(4, 2, 7)#, sharey=axs[2, 0])  # 4 rows, 2 columns, 7th subplot
    
    # Adjust the position of the new subplot
    box = axs[2, 0].get_position()
    ax3.set_position([box.x0, box.y0 - box.height * 0.55, box.width, box.height * 0.2])
    # ax3.plot(Z[1:], loss_per_cell_array, '--', color='black', label='Losses per cell, %')
    ax3.set_xlim([np.min(Z), np.max(Z)])
    # Hide the spines, ticks, and labels of the new subplot
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.yaxis.set_ticks([])
    ax3.yaxis.set_label_coords(-0.1, 0.5)
    
    # Set any other properties of ax3 here
    ax3.set_xlabel('z, m', fontsize=fontsize)  # Add a label for the new x-axis
    
    labelss = ['0', ' ', '20', ' ', '40', '', '60', '', '80', '', '100']
    ax2.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ax2.set_yticklabels(labelss)
    ax2.set_xlabel('iris number', fontsize=fontsize) 
    
    plt.tight_layout()
       
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        # _logger.debug(ind_str + 'saving *{:}.{:}'.format(suffix, savefig))
        fig.savefig(filePath + suffix + '.' + str(savefig), format=savefig, dpi=300)
    # _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))

    plt.show()

#%%

THz = 3
xlamds = speed_of_light/THz/1e12

N = 2000
n_iter = N
b = 0.3
a = 0.055

dfl = generate_gaussian_dfl(xlamds=xlamds, shape=(201, 201, 1), dgrid=(3*a, 3*a, 10), power_rms=(205.5e-2, 205.5e-2, 1),
                          power_center=(0, 0, None), power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                          zsep=None, freq_chirp=0, en_pulse=None, power=1e6)
Nf = a**2/dfl.xlamds/b
print('Nf =', Nf)
    
# plot_dfl(dfl,fig_name='dfl at source domain', domains = "fs")
# plot_dfl(dfl, fig_name='dfl at source k domain', domains = "fk")

dfl_iris_exit = deepcopy(dfl)

#%%
dfl_iris_exit, E_x1, E_y1, rad_left_array1, loss_per_cell_array1 = dfl_prop_iris(dfl, N=N, a=a, b=b, 
                                                                                 n_iter=n_iter, acount_first_cell_loss=1)
#%%
plot_iris_line_prop_result(E_x1, E_y1, dfl_iris_exit, a, b, N, n_iter, rad_left_array1, loss_per_cell_array1, 
                               direction = 'x')

#%%

N = 1000
n_iter = N
b = 0.3
a = 0.055

dfl_waveguide2 = deepcopy(dfl_iris_exit)

rad_left_array2 = 0 
loss_per_cell_array2 = 0

dfl_waveguide2, E_x2, E_y2, rad_left_array2, loss_per_cell_array2 = dfl_prop_iris(dfl_iris_exit, N=N, a=a, b=b, 
                                                                              n_iter=n_iter, acount_first_cell_loss=1)

#%%
plot_iris_line_prop_result(E_x2, E_y2, dfl_waveguide2, a, b, N, n_iter, rad_left_array2, loss_per_cell_array2, 
                               direction = 'x')

#%%

dfl_waveguide2.to_domain('sf')
x = dfl_waveguide2.scale_x()

M = (8 * np.pi * a**2 / dfl_waveguide2.xlamds / b)**(-1/2)
beta0 = 0.824 

v_nj = special.jn_zeros(0, 1)[-1]
k_nj = v_nj*(1 - (1 + 1j)*beta0*M)/a

E = special.jv(0, k_nj * x) #* np.exp(-1j*k_z * z)
E_gauss = np.exp(-(k_nj)**2 * x**2 / 4)
V_max = np.max(abs(E)**2/np.sum(abs(E)**2))

plt.figure()

plt.plot(x/a, abs(E)**2/np.sum(abs(E)**2) / V_max, label='Fundamental mode')
plt.plot(x/a, abs(E_gauss)**2/np.sum(abs(E_gauss)**2) / V_max, color='black', linestyle='--', label='Gaussian')
plt.plot(x/a, dfl_waveguide2.intensity()[0, dfl_waveguide2.Ny()//2, :]/np.sum(dfl_waveguide2.intensity()[0, dfl_waveguide2.Ny()//2, :])/V_max, 
         color='red', label='Numerical')

plt.xlim(0,  1)
plt.ylim(0)
plt.xlabel('r/a', fontsize=18)
plt.ylabel('arb. units', fontsize=18)
plt.tight_layout()
plt.legend()
plt.show()






