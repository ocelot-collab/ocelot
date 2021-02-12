#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 14:25:09 2020

@author: tae
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
import copy
from scipy import misc, signal, special


from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.gui.genesis_plot import *

from ocelot.common.globals import *
from ocelot.common.math_op import find_nearest_idx, fwhm, std_moment, bin_scale, bin_array, mut_coh_func
from ocelot.common.py_func import filename_from_path

from ocelot.common.ocelog import *
_logger = ocelog
# _logger = logging.getLogger(__name__)


#### TODO ####
#1 generate the field in 'f' domain #done
#2 whats whrong eith shift fucntion????
#3 .'far_field'.replace('_', '')
#4 line 508 garammar
#5 522 E = 0 !bad!
#6 524-537 remove coh incoh conditions -- just check if phi
#7  
def write_script(script_dir, new_script_dir):
    # save of your python script in the simulation directory
    
    f = open(script_dir, 'r')
    script = []
    for line in f:
        script.append(str(line))
    f.close()
    
    f_out = open(new_script_dir, 'w')
    f_out.writelines(script)
    f_out.close()

def undulator_field_far(theta_x, theta_y, l_x, l_y, eta_x, eta_y, z, L_w, E_ph, phi=False): 
    '''
    Undulator field generator for the far field approximation
    
    ### field expressions look at the following article
    Geloni, Gianluca, et al. "Fourier treatment of near-field synchrotron radiation theory." Optics communications 276.1 (2007): 167-179.
    ###
    
    Parameters
    ----------
    theta_x : float
        grid along x having a shape [Nx, Ny]
    theta_y : float
        grid along y having a shape [Nx, Ny]
    l_x : float
        random x offset of an electron
    l_y : float
        random y offset of the electron
    eta_x : float
        random x deflection angle of the electron
    eta_y : float
        random y deflection angle of the electron
    L_w : float, optional
        Undulator length in meteres
    E_ph : int, optional
        Photon energy in eV
    phi : double. Default is None 
        longitudinal phase of the electron. If not None -- long beam approximation: fields sum up with a random phase; 
        if None -- short beam approximation: fields sum up in phase.

    Returns
    -------
    an array with a shape (Nx, Ny)
        the field emitted by the electron
    '''
    w = E_ph / hr_eV_s 

    _logger.debug(ind_str + 'calculating undulator field in the far field approximation')
    
    if z != 0: 
        z_0 = z
    else: 
        raise AttributeError('"phi" is bool type ')
        _logger.error('"z" must not be zero')

    
    start_time = time.time()
    E = -np.sinc(L_w * w * ((theta_x - (l_x/z_0) - eta_x)**2 + (
                            theta_y - (l_y/z_0) - eta_y)**2) / 4 /speed_of_light / np.pi) * np.exp(
                            1j * w  * z_0 * ((theta_x - (l_x/z_0))**2 + (theta_y - (l_y/z_0))**2) / 2 / speed_of_light)
    if phi != None:
        _logger.debug(ind_str + 'adding a random phase to the field')
        E = E * np.exp(1j * phi)
    
    _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
    return E        

def undulator_field_near(x, y, l_x, l_y, eta_x, eta_y, z, L_w, E_ph, phi=None):
    '''
    Field generator for near field approximation
    
    ### field expressions look at the following article
    Geloni, Gianluca, et al. "Fourier treatment of near-field synchrotron radiation theory." Optics communications 276.1 (2007): 167-179.
    ###
    
    Parameters
    ----------
    theta_x : float
        grid along x
    theta_y : float
        grid along y
    l_x : float
        x offset of an electron
    l_y : float
        y offset of the electron
    eta_x : float
        x deflection angle of the electron
    eta_y : float
        y deflection angle of the electron
    phi : double. Default is None 
        longitudinal phase of the electron. The default is None -- short beam approximation: fields sum up in phase. 
        If not None -- long beam approximation: fields sum up with a random phase. 

    Returns
    -------
    array with a shape (Nx, Ny)
        the field emitted by the electron

    '''
    w = E_ph / hr_eV_s 

    if z != 0: 
        z_0 = z
    else: 
        raise AttributeError('"phi" is bool type ')
        _logger.error('"z" must not be zero')
   
    _logger.debug('calculating undulator field in the near field approximation')
    start_time = time.time() 
    
    E = (sc.expi(1j * w * ((x - l_x - z_0*eta_x)**2 + (y - l_y - z_0*eta_y)**2) / (2 * z_0 * speed_of_light - L_w*speed_of_light)) - sc.expi(
        1j * w * ((x - l_x - z_0*eta_x)**2 + (y - l_y - z_0*eta_y)**2) / (2 * z_0 * speed_of_light + L_w*speed_of_light))) * np.exp(
            1j * w * (((x - l_x)**2 + (y - l_y)**2) - (x - l_x - z_0*eta_x)**2 - (y  - l_y - z_0*eta_y)**2) / 2 / speed_of_light / z_0)  
    
    # the field has a singularity at r = 0 and z_0 = L_w
    ### this block is for excluding nan value
    _logger.debug(ind_str + 'looking for None value because of features at r = 0 and z_0 = L_w, may appear for filament beams')
    none_indx = np.where(np.isnan(E))
    if not none_indx:    
        _logger.WARNING(ind_str + "none_indx = {}".format(none_indx))
        avr = E[(none_indx[0] - 1)[0]:(none_indx[0] + 2)[0],(none_indx[1] - 1)[0]:(none_indx[1] + 2)[0]]
        avr = np.average(avr[~np.isnan(avr)])
        E[none_indx] = avr
        _logger.debug(ind_str + 'nan value found')
    else:
        _logger.debug(ind_str + 'no nan value found')
    ###
    
    if phi != None:
        _logger.debug(ind_str + 'adding a random phase to the field')
        _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
        return E * np.exp(1j * phi)
    else:
        _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
        return E

def undulator_field_dfl_SP(dfl, z, L_w, E_ph, N_e=1, sig_x=0, sig_y=0, sig_xp=0, sig_yp=0, approximation='far_field', mode='incoh', seed=None):
    
    """  
    SRW-like Undulator radaition generator - single particle radiation is propagated through the optical system
    
    Parameters
    ----------
    dfl : RadiationField object
        3d or 2d coherent radiation distribution, *.fld variable is the same as Genesis dfl structure
        
        Transverse dfl.fld.shape defines the field dimension -- (Nx, Ny).
        
        Returned dfl.fld represents a simulation where longitudibnal scale is N_e single electron fields, 
        each field has its own amplitude and the phase; and propagets seperatly. 
        
        After propagation intesities should be summed up to find the resulting field 
    
    z : float
        Distance from source at which the field is calculated
    N_e : int, optional
        Number of macro particles in an electron bunch. The default is 1.
    L_w : float, optional
        Undulator length in meteres
    E_ph : int, optional
        Photon energy in eV
    sig_x : float, optional
        r.m.s. of the x electron offset assuming Gaussian distribution. The default is 0.
    sig_y : float, optional
        r.m.s. of the y electron offset assuming Gaussian distribution. The default is 0.
    sig_xp : float, optional
        r.m.s. of the x electron deflection assuming Gaussian distribution. The default is 0.
    sig_yp : float, optional
        r.m.s. of the y electron deflection assuming Gaussian distribution. The default is 0.
    approximation : str, optional
        Used approximation to calculate the field ('near or farfield'). The default is 'far_field'.
    mode : str, optional
       A mode to sum up the field. For the SRW generator physical meaning has only incoherent addition of fields. The default is 'incoh'.
    seed : int, optional
        A seed for the random generator. The default is None.

    Returns
    -------
    dfl : RadiationField object with dfl.shape = (N_e, dfl.Ny(), dfl.Nx())
    """
        
    _logger.info('Calculating undulator field with SRW-like approach')
    _logger.info(ind_str + 'used approximation is ' + approximation.replace('_', ' '))

    start_time = time.time()

    if approximation == 'near_field':
        dfl.to_domain(domains='sf', indexing='ij') 
        x = dfl.scale_kx()
        y = dfl.scale_ky()
        x, y = np.meshgrid(x, y)    
    elif approximation == 'far_field':
        dfl.to_domain(domains='sf', indexing='ij') 
        kx = dfl.scale_x()/z
        ky = dfl.scale_y()/z
        theta_x, theta_y  = np.meshgrid(kx, ky)
    else: 
        _logger.error(ind_str + '"approximation" must be whether "near_field" of "far_field"')
        raise AttributeError('"approximation" must be whether "near_field" of "far_field"')

    if seed != None:
        random.seed(seed)
        _logger.debug(ind_str + 'seed is {}'.format(seed))
    
    l_x   = np.random.normal(0, sig_x, (N_e))
    eta_x = np.random.normal(0, sig_xp, (N_e))
    l_y   = np.random.normal(0, sig_y, (N_e))
    eta_y = np.random.normal(0, sig_yp, (N_e))
    phi = np.random.uniform(low=0, high=2*np.pi, size=N_e)

    E_list = []
    i=1
    for l_xi, l_yi, eta_xi, eta_yi, phi_i in zip(l_x, l_y, eta_x, eta_y, phi):
        if approximation == 'far_field':
            if mode == 'incoh':
                E = undulator_field_far(theta_x, theta_y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph, phi=phi_i)
        elif approximation == 'near_field':
            if mode == 'incoh':
                E = undulator_field_near(x, y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph, phi=phi_i)
        else: 
            _logger.error(ind_str + 'attribute "method" must be "far_field" or "near_field"')
            raise AttributeError('attribute "method" must be "far_field" or "near_field"') 
        E_list.append(E)     
        i = i + 1
        print('fields simulated = {0} out of {1}'.format(i, N_e))
    
    E_list = np.array(E_list)
    dfl.fld = E_list
    
    _logger.debug(ind_str + 'undulator field calculated in {:.2f} seconds'.format(time.time() - start_time))
    _logger.info(ind_str + 'done')

    return dfl

def undulator_field_dfl_MP(dfl, z, L_w, E_ph, N_b=1, N_e=1, sig_x=0, sig_y=0, sig_xp=0, sig_yp=0, approximation='far_field', mode='incoh', seed=None):
    
    """  
    Multiple particle Undulator radaition generator - fields from numerous macroparticles is added and propagated as one
    
    Parameters
    ----------
    dfl : RadiationField object
        3d or 2d coherent radiation distribution, *.fld variable is the same as Genesis dfl structure
        
        Transverse dfl.fld.shape defines the field dimension -- dfl.shape = (N_e, dfl.Ny(), dfl.Nx())
        
        Returned dfl.fld represents a simulation where eahc longitudinal slice is a statistical realization of a field from N_e macro electrons, 
        each realization has its own amplitude and phase; and propagets seperatly. 
        
        After propagation of the realizations the intensities should be summed up to find the resulting field.
    
    z : float
        Distance at which calculate the field
    N_b : int, optional
        Number of realizations to generate. The default is 1.
    N_e : int, optional
        Number of macro electrons in an electron bunch. The default is 1.
    L_w : float, optional
        Undulator length in meteres
    E_ph : int, optional
        Photon energy in eV
    sig_x : float, optional
        r.m.s. of the x electron offset for the gaussian distribution. The default is 0.
    sig_y : float, optional
        r.m.s. of the y electron offset for the gaussian distribution. The default is 0.
    sig_xp : float, optional
        r.m.s. of the x electron deflection for the gaussian distribution. The default is 0.
    sig_yp : float, optional
        r.m.s. of the y electron deflection for the gaussian distribution. The default is 0.
    approximation : str, optional
        Used approximation to culculate the field ('near or farfield'). The default is 'far_field'.
    mode : str, optional
       A mode to sum up the field.
    seed : int, optional
        A seed for the random generator. The default is None.

    Returns
    -------
    dfl : RadiationField object
    """
    _logger.info(ind_str + 'Calculating undulator field with Monte Carlo method...')
    _logger.info(ind_str + 'used approximation is ' + approximation.replace('_', ' '))

    start_time = time.time()

    if approximation == 'near_field':
        dfl.to_domain(domains='sf', indexing='ij') 
        x = dfl.scale_kx()
        y = dfl.scale_ky()
        x, y = np.meshgrid(x, y)    
    elif approximation == 'far_field':
        dfl.to_domain(domains='sf', indexing='ij') 
        kx = dfl.scale_x()/z
        ky = dfl.scale_y()/z
        theta_x, theta_y  = np.meshgrid(kx, ky)
    else: 
        _logger.error(ind_str + '"approximation" must be whether "near_field" of "far_field"')
        raise AttributeError('"approximation" must be whether "near_field" of "far_field"')
  
    if seed != None:
        random.seed(seed)
        _logger.info('seed is {}'.format(seed))
    
    E_list = []   
    for bunch in range(N_b):
        l_x   = np.random.normal(0, sig_x, (N_e))
        eta_x = np.random.normal(0, sig_xp, (N_e))
        l_y   = np.random.normal(0, sig_y, (N_e))
        eta_y = np.random.normal(0, sig_yp, (N_e))
        phi = np.random.uniform(low=0, high=2*np.pi, size=N_e)
        
        E = np.zeros((dfl.shape()[1], dfl.shape()[2]))
        i=1
        for l_xi, l_yi, eta_xi, eta_yi, phi_i in zip(l_x, l_y, eta_x, eta_y, phi):
            if approximation == 'far_field':
                if mode == 'incoh':
                    E = E + undulator_field_far(theta_x, theta_y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph, phi=phi_i)   
                elif mode == 'coh':
                    E = E + undulator_field_far(theta_x, theta_y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph)
            elif approximation == 'near_field':
                if mode == 'incoh':
                    E = E + undulator_field_near(x, y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph, phi=phi_i)
                elif mode == 'coh':
                    E = E + undulator_field_near(x, y, l_xi, l_yi, eta_xi, eta_yi, z, L_w=L_w, E_ph=E_ph)
            else: 
                _logger.error(ind_str + 'attribute "method" must be "far_field" or "near_field"')
                raise AttributeError('attribute "method" must be "far_field" or "near_field"')
            if N_b == 1:
                i = i + 1
                _logger.debug(ind_str + 'electrons simulated = {0} out of {1}'.format(i, N_e))
        E = E/N_e
        E_list.append(E)
        # _logger.debug(ind_str + 'realizations number = {0} out of {1}'.format(bunch, N_b))
        _logger.info(ind_str + 'realizations number = {0} out of {1}'.format(bunch, N_b))
    
    E_list = np.array(E_list)
    dfl.fld = E_list

    _logger.info(ind_str + 'undulator field calculated in {:.2f} seconds'.format(time.time() - start_time))

    return dfl

# def SERVAL_undulator_field_generator(dfl, L_w, sig_x=0, sig_y=0, sig_xp=0, sig_yp=0, showfig=False, seed=None):

#     w_0 = 2*np.pi * speed_of_light / dfl.xlamds
    
#     if showfig:
#         plot_dfl(dfl, line_off_xy = False, fig_name = '1-noise')
    
#     dfl.to_domain('sf')
    
#     x, y = np.meshgrid(dfl.scale_x(), dfl.scale_y())#, indexing='ij')
    
#     mask_xy_ebeam = np.exp(- x**2 / 4 / sig_x**2 - y**2 / 4 / sig_y**2) # 4 because amplitude, not intensity
    
#     mask_xy_ebeam /= np.sum(mask_xy_ebeam)
    
#     # mask_xy_radiation = np.sqrt((1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0]))**2)
#     mask_xy_radiation = 1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0])

#     mask_xy = signal.fftconvolve(mask_xy_radiation, mask_xy_ebeam, mode='same')
#     mask_xy = mask_xy_ebeam
    
#     ocelog.info('Multiplying by near field mask')
#     dfl.fld *= mask_xy
#     # dfl.fld *= np.sqrt(mask_xy)
#     ocelog.info(ind_str +'done')
    
#     if showfig:
#         plot_dfl(dfl, domains='s', line_off_xy = False, fig_name = '2-X_e-beam-size')
#         plot_dfl(dfl, domains='k', line_off_xy = False, fig_name = '2-X_e-beam-size')
        
#     dfl.to_domain('kf')
     
#     k_x, k_y = np.meshgrid(dfl.scale_x(), dfl.scale_y())
    
    
#     mask_kxky_ebeam = np.exp(-k_y**2 / 4 / sig_yp**2 - k_x**2 / 4 / sig_xp**2 ) # 4 because amplitude, not intensity
    
#     mask_kxky_ebeam /= np.sum(mask_kxky_ebeam)
    
        
#     # mask_kxky_radiation = np.sqrt((np.sinc(w_0 * L_w * (k_x**2 + k_y**2) / 4 / speed_of_light / np.pi))**2)# Geloni2018 Eq.3, domega/omega = 2dgamma/gamma, divided by pi due to np.sinc definition
#     # mask_kxky_radiation = (np.sinc(w_0 * L_w * (k_x**2 + k_y**2) / 4 / speed_of_light / np.pi))# Geloni2018 Eq.3, domega/omega = 2dgamma/gamma, divided by pi due to np.sinc definition

#     mask_kxky_radiation = np.fft.ifftshift(mask_xy_radiation)#, axes=(0,1))
#     mask_kxky_radiation = np.fft.ifft2(mask_kxky_radiation)#, axes=(0,1))
#     # mask_kxky_radiation = np.fft.ifftshift(mask_kxky_radiation)#, axes=(0,1))
#     mask_kxky_radiation /= np.sqrt(dfl.Nx() * dfl.Ny())
    
#     if True:
#         plt.figure()
#         plt.pcolormesh(np.real(mask_kxky_radiation)**2 + np.imag(mask_kxky_radiation)**2)
#         plt.show()

#     if True:
#         plt.figure()
#         plt.pcolormesh(np.real(mask_kxky_ebeam)**2 + np.imag(mask_kxky_ebeam)**2)
#         plt.show()

#     mask_kxky_radiation = np.fft.fftshift(mask_kxky_radiation)    
#     mask_kxky = signal.fftconvolve(mask_kxky_ebeam, mask_kxky_radiation, mode='same')

#     mask_kxky /= np.sum(mask_kxky)
    
#     if True:
#         plt.figure()
#         plt.pcolormesh(np.real(mask_kxky)**2 + np.imag(mask_kxky)**2)
#         plt.show()
    
#     dfl.fld *= mask_kxky[np.newaxis, :, :]
#     plot_dfl(dfl, domains='k', fig_name = '4-X_e-beam-divergence_after')
#     plot_dfl(dfl, domains='s', fig_name = '4-X_e-beam-divergence_after')

#     # dfl.fld *= np.sqrt(mask_kxky[np.newaxis, :, :])

#     if showfig:
#         passs
#         # plot_dfl(dfl, domains='s', fig_name = '4-X_e-beam-divergence')
#         # plot_dfl(dfl, domains='k', fig_name = '4-X_e-beam-divergence')

#     return dfl 

def undulator_field_dfl_SERVAL(dfl, L_w, sig_x=0, sig_y=0, sig_xp=0, sig_yp=0, k_support = 'intensity', s_support='conv', showfig=False, seed=None):
    
    _logger.info('Generating undulator field with Serval algorithm')
    w_0 = 2*np.pi * speed_of_light / dfl.xlamds
    
    if showfig:
        plot_dfl(dfl, line_off_xy = False, fig_name = '1-X_noise')
    
    dfl.to_domain('sf')
    
    x, y = np.meshgrid(dfl.scale_x(), dfl.scale_y())#, indexing='ij')
    
    mask_xy_ebeam = np.exp(- x**2 / 4 / sig_x**2 - y**2 / 4 / sig_y**2) # 4 because amplitude, not intensity
    # mask_xy_ebeam = np.exp(- x**2 / 2 / sig_x**2 - y**2 / 2 / sig_y**2) # 4 because amplitude, not intensity

    mask_xy_ebeam /= np.sum(mask_xy_ebeam)
    
    # mask_xy_radiation = np.sqrt((1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0]))**2)
    # mask_xy_radiation = 1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0])

    # mask_xy_radiation = (1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0]))**2
    if s_support == 'conv':
        _logger.info(ind_str +'s_support == "conv"')
        mask_xy_radiation = 1j*(np.pi - 2*special.sici(w_0*(x**2 + y**2)/speed_of_light/L_w)[0])
        mask_xy = signal.fftconvolve(mask_xy_radiation, mask_xy_ebeam, mode='same')
    else:
        _logger.info(ind_str +'s_support == "beam"')
        mask_xy = mask_xy_ebeam
    
    _logger.info(ind_str +'Multiplying by real space mask')
    dfl.fld *= mask_xy
    # dfl.fld *= np.sqrt(mask_xy)
    _logger.info(2*ind_str +'done')

    if showfig:
        plot_dfl(dfl, domains='s', line_off_xy = False, fig_name = '2-X_e-beam-size')
        plot_dfl(dfl, domains='k', line_off_xy = False, fig_name = '2-X_e-beam-size')
                
    dfl.to_domain('kf')

    k_x, k_y = np.meshgrid(dfl.scale_x(), dfl.scale_y())
    mask_kxky_ebeam = np.exp(-k_y**2 / 4 / sig_yp**2 - k_x**2 / 4 / sig_xp**2 ) # 4 because amplitude, not intensity
    # mask_kxky_ebeam = np.exp(-k_y**2 / 2 / sig_yp**2 - k_x**2 / 2 / sig_xp**2 ) # 2 because intensity
    mask_kxky_ebeam /= np.sum(mask_kxky_ebeam)
    
    # mask_kxky_radiation = np.sqrt((np.sinc(w_0 * L_w * (k_x**2 + k_y**2) / 4 / speed_of_light / np.pi))**2)# Geloni2018 Eq.3, domega/omega = 2dgamma/gamma, divided by pi due to np.sinc definition
    # mask_kxky_radiation = (np.sinc(w_0 * L_w * (k_x**2 + k_y**2) / 4 / speed_of_light / np.pi))# Geloni2018 Eq.3, domega/omega = 2dgamma/gamma, divided by pi due to np.sinc definition
        
    mask_kxky_radiation = np.sinc(w_0 * L_w * (k_x**2 + k_y**2) / 4 / speed_of_light / np.pi)# Geloni2018 Eq.3, domega/omega = 2dgamma/gamma, divided by pi due to np.sinc definition

    if k_support == 'intensity':
        _logger.info(ind_str +'k_support == "intensity"')
        mask_kxky = signal.fftconvolve(mask_kxky_ebeam**2, mask_kxky_radiation**2, mode='same')
        mask_kxky = np.sqrt(mask_kxky[np.newaxis, :, :])
        mask_kxky /= np.sum(mask_kxky)
    elif k_support == 'amplitude':
        _logger.info(ind_str +'k_support == "amplitude"')
        mask_kxky = signal.fftconvolve(mask_kxky_ebeam, mask_kxky_radiation, mode='same')
        mask_kxky /= np.sum(mask_kxky)
    else:
        raise ValueError('k_support should be either "intensity" or "amplitude"')
    
    # dfl.fld *= mask_kxky[np.newaxis, :, :]
    _logger.info(ind_str +'Multiplying by inverse space mask')
    dfl.fld *= mask_kxky
    _logger.info(2*ind_str +'done')

    if showfig:
        plot_dfl(dfl, domains='s', fig_name = '3-X_radaition_size')
        plot_dfl(dfl, domains='k', fig_name = '3-X_radiation_divergence')

    return dfl 









