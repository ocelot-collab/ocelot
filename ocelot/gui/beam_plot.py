

import sys
import os
import csv
import time
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import logging


from ocelot.adaptors.genesis import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.common.math_op import *  # import of mathematical functions
from ocelot.utils.xfel_utils import *
from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.optics.wave import *
from ocelot.gui.settings_plot import *
from ocelot.cpbd.beam import BeamFormFactor

_logger = logging.getLogger(__name__)

# from pylab import rc, rcParams #tmp
# from matplotlib import rc, rcParams
# from mpl_toolkits.axes_grid1 import make_axes_locatable

# def_cmap = 'viridis'
# # def_cmap = 'Greys'

# fntsz = 4
# params = {'image.cmap': def_cmap, 'backend': 'ps', 'axes.labelsize': 3 * fntsz, 'font.size': 3 * fntsz, 'legend.fontsize': 4 * fntsz, 'xtick.labelsize': 4 * fntsz,  'ytick.labelsize': 4 * fntsz, 'text.usetex': False}
# rcParams.update(params)

@if_plottable
def plot_beam(beam, figsize=3, showfig=True, savefig=False, fig=None, plot_xy=None, debug=0):
    _logger.info('plotting beam')

    if showfig == False and savefig == False:
        return

    # if beam.__class__ != GenesisBeam:
    #     raise ValueError('wrong beam object: should be GenesisBeam')

    fontsize = 15

    if plot_xy == None:
        if np.mean(beam.x) == 0 and np.mean(beam.y) == 0 and np.mean(beam.xp) == 0 and np.mean(beam.yp) == 0:
            plot_xy = 0
        else:
            plot_xy = 1

    if fig == None:
        fig = plt.figure()
    fig.clf()

    idx = beam.idx_max()

    g0 = np.mean(beam.g).astype(int)  # mean
    g_dev = beam.g - g0  # deviation from mean

    fig.set_size_inches((4 * figsize, (3 + plot_xy) * figsize), forward=True)
    ax_I = fig.add_subplot(2 + plot_xy, 2, 1)
    plt.grid(True)
    ax_I.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.I, 'r', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.I[idx], 'bs')
    ax_I.set_ylim(ymin=0)

    if hasattr(beam, 'eloss'):
        if (beam.eloss != 0).any():
            ax_loss = ax_I.twinx()
            p2, = plt.plot(1.e6 * np.array(beam.s), 1.e-3 * np.array(beam.eloss), 'g', lw=3)
            ax_loss.legend([p1, p2], [r'$I [A]$', r'Wake $[KV/m]$'], fontsize=fontsize, loc='best')
            ax_loss.set_ylim(auto=True)
        else:
            ax_I.legend([r'$I [A]$'], fontsize=fontsize, loc='best')
    else:
        ax_I.legend([r'$I [A]$'], fontsize=fontsize, loc='best')

    ax_I.text(0.02, 0.98, r'Q = {:.2f} pC'.format(beam.charge() * 1e12), horizontalalignment='left',
              verticalalignment='top', transform=ax_I.transAxes)
    # ax.set_xlim([np.amin(beam.s),np.amax(beam.x)])
    ax = fig.add_subplot(2 + plot_xy, 2, 2, sharex=ax_I)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
    # p1,= plt.plot(1.e6 * np.array(beam.s),1.e-3 * np.array(beam.eloss),'r',lw=3)
    p1, = plt.plot(1.e6 * np.array(beam.s), g_dev, 'r', lw=3)
    plt.plot(1.e6 * beam.s[idx], g_dev[idx], 'bs')
    ax = ax.twinx()
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.dg, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.dg[idx], 'bs')
    ax.legend([p1, p2], [r'$\gamma$ + ' + str(g0), r'$\delta \gamma$'], loc='best')

    ax = fig.add_subplot(2 + plot_xy, 2, 3, sharex=ax)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.emit_xn * 1e6, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.emit_yn * 1e6, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.emit_xn[idx] * 1e6, 'bs')
    ax.set_ylim(ymin=0)

    ax.legend([p1, p2], [r'$\varepsilon_x [\mu m]$', r'$\varepsilon_y [\mu m]$'], fontsize=fontsize, loc='best')
    # ax3.legend([p3,p4],[r'$\varepsilon_x$',r'$\varepsilon_y$'])

    ax = fig.add_subplot(2 + plot_xy, 2, 4, sharex=ax)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.beta_x, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.beta_y, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.beta_x[idx], 'bs')
    ax.set_ylim(ymin=0)
    ax.legend([p1, p2], [r'$\beta_x [m]$', r'$\beta_y [m]$'], fontsize=fontsize, loc='best')

    if plot_xy:
        ax = fig.add_subplot(3, 2, 5, sharex=ax)
        plt.grid(True)
        ax.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
        p1, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.x), 'r', lw=3)
        p2, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.y), 'g', lw=3)

        ax.legend([p1, p2], [r'$x [\mu m]$', r'$y [\mu m]$'], fontsize=fontsize, loc='best')

        # beam_beta = sqrt(1 - (1 / beam.g0**2))
        # beam_p = beam.g0 * beam_beta
        # # p=beam.g0*m_e_eV/speed_of_light
        # pz = sqrt(beam_p**2 - beam.px**2 - beam.py**2)
        # xp = beam.px / pz
        # yp = beam.py / pz

        ax = fig.add_subplot(3, 2, 6, sharex=ax)
        plt.grid(True)
        ax.set_xlabel(r'$s [\mu m]$', fontsize=fontsize)
        p1, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.xp), 'r', lw=3)
        p2, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.yp), 'g', lw=3)

        ax.legend([p1, p2], [r'$x_p [\mu rad]$', r'$y_p [\mu rad]$'], fontsize=fontsize, loc='best')

    ax.set_xlim([1.e6 * np.amin(beam.s), 1e6 * np.amax(beam.s)])

    fig.subplots_adjust(hspace=0.2, wspace=0.3)

    # if savefig != False:
    #     if hasattr(beam,'filePath'):
    #         if savefig == True:
    #             filetype = 'png'
    #         else:
    #
    #         path = beam.filePath
    #         name = beam.fileName()
    #     else:
    #         if if savefig == True:
    #
    #         if debug > 1:
    #             print('      saving ' + beam.fileName() + '.' + savefig)
    #         plt.savefig(beam.filePath + '.' + savefig, format=savefig)

    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving ' + beam.filePath + '.' + savefig)
        plt.savefig(beam.filePath + '.' + savefig, format=savefig)

    if showfig:
        rcParams["savefig.directory"] = os.path.dirname(beam.filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)

    _logger.debug(ind_str + 'done')


def plot_beam_form_factor(form_factor: BeamFormFactor, x_axis_units='THz', x_axis_limits=(None, None), y_axis_limits=(None, None),
                          figsize=(2 * 3, 3 * 3), normalize=True, fig_name='Beam Form Factor', log_scale=False, showtext=True,
                          savefig=False, showfig=True, **kwargs):
    '''
    plots a BeamFormFactor object
    '''
    if not showfig and not savefig:
        return

    start_time = time.time()
    _logger.info('plotting Beam Form Factor')
    fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches(figsize, forward=True)

    ax_current_profile = fig.add_subplot(2, 1, 1)
    s, current = form_factor.beam_array.s, form_factor.beam_array.I
    ax_current_profile.plot(s * 1e6, current * 1e-3, **kwargs)
    ax_current_profile.set_title('Current profile', fontsize=15)
    ax_current_profile.set_xlabel('s [$\mu$m]')
    ax_current_profile.set_ylabel('peak current [kA]')
    if showtext:
        ax_current_profile.text(0.97, 0.97, 'Q={:.2f}pC'.format(form_factor.beam_array.charge()*1e12), 
                                horizontalalignment='right', verticalalignment='top', transform=ax_current_profile.transAxes, fontsize=12, color='black')

    form_factor.calc()
    frequency, ffactor =  form_factor.frequency, form_factor.modulus
    if normalize:
        ffactor /= ffactor[0]
        
    xlable = '$\omega$ [Hz]'
    if x_axis_units in ['MHz', 'mhz']:
        frequency = frequency * 1e-6
        xlable = '$\omega$ [MHz]'
    elif x_axis_units in ['GHz', 'ghz']:
        frequency = frequency * 1e-9
        xlable = '$\omega$ [GHz]'
    elif x_axis_units in ['THz', 'thz', 'f']:
        frequency = frequency * 1e-12
        xlable = '$\omega$ [THz]'
    elif x_axis_units in ['keV', 'kev']:
        frequency = frequency * 2 * np.pi * hr_eV_s * 1e-3
        xlable = '$\epsilon$ [keV]'
    elif x_axis_units in ['eV', 'ev']:
        frequency = frequency * 2 * np.pi * hr_eV_s
        xlable = '$\epsilon$ [eV]'
    elif x_axis_units in ['meV', 'mev']:
        frequency = frequency * 2 * np.pi * hr_eV_s * 1e3
        xlable = '$\epsilon$ [meV]'
    else:
        _logger.warning(ind_str + 'x_axis_units must one of \'MHz\', \'GHz\', \'THz\', \'keV\', \'eV\', \'meV\'')

    ax_form_factor = fig.add_subplot(2, 1, 2)
    if log_scale:
        ax_form_factor.semilogy(frequency, ffactor**2, **kwargs)
    else:
        ax_form_factor.plot(frequency, ffactor**2, **kwargs)
    ax_form_factor.set_title('Form factor profile', fontsize=15)
    ax_form_factor.set_ylabel('|form factor|^2 [a. u.]')
    ax_form_factor.set_xlabel(xlable)
    ax_form_factor.set_xlim(*x_axis_limits)
    ax_form_factor.set_ylim(*y_axis_limits)
    
    fig.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.draw()
    _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
    if savefig:
        _logger.debug(ind_str + 'saving figure ' + fig_name)
        if not isinstance(savefig, str):
            fig.savefig(fig_name + '.png')
        else:
            fig.savefig(fig_name + '.' + savefig)

    if showfig:
        _logger.debug('showing figure ' + fig_name)
        plt.show()
    else:
        plt.close('all')


@if_plottable
@save_show
def plot_estimator_power_z(fel, z=None, fig=None, und_duty_factor=1, fs=False):
    '''
    plots estimated FEL power at position z
    fel is FelParamterArray object
    und_duty_factor <= 1
    '''
    if fig == None:
        fig = plt.figure()
    fig.clf()
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
    ax.set_title('Pulse energy @ %.2fm = %.2e J' %(z / und_duty_factor, fel.E(z)))
    ax.set_ylabel('P [W]')
    return fig
    
    #fig.savefig(exp_dir + 'power_sat.png', format = 'png')
    return fig

@if_plottable
@save_show
def plot_estimator_spectrogram(fel, z=None, fig=None, cmap='Reds', **kwargs):
    
    if fig == None:
        fig = plt.figure()
    fig.clf()
    
    ax = fig.gca()
    
    # if z is None:
        # z = fel.z_sat_min
    # Psat = fel.P(z)
    # Psat[np.isnan(Psat)]=0
    # idx = fel.idx

    # phen0 = fel.phen0
    # dphen = phen0 * fel.rho3
    # dp = dphen[idx] / 10
    # s_arr = fel.s * 1e6
    # phen_arr = np.arange(np.amin(phen0 - 3 * dphen), np.amax(phen0 + 3 * dphen), dp)
    # spec = np.zeros((s_arr.size, phen_arr.size))
    # for i in range(s_arr.size):
        # if dphen[i] != 0:
            # spec[i] = np.exp(-(phen_arr - phen0[i])**2 / 2 / dphen[i]**2) / np.sqrt(2 * np.pi * dphen[i]**2)
    # spec = spec * Psat[:, np.newaxis]
    
    s_arr, phen_arr, spectrogram = fel.spectrogram(z = z)
    ax.pcolormesh(s_arr * 1e6, phen_arr, spectrogram, cmap=cmap)
    ax.set_xlabel('s [um]')
    ax.set_ylabel('E [eV]')
    ax.autoscale(tight=1)
    return fig

@if_plottable
@save_show
def plot_estimator_spectrum(fel, z=None, fig=None, **kwargs):
    
    if fig == None:
        fig = plt.figure()
    fig.clf()
    
    ax = fig.gca()
    
    phen_arr, spectrum = fel.spectrum(z = z)
    ax.plot(phen_arr, spectrum/np.amax(spectrum))
    ax.set_xlabel('E [eV]')
    ax.set_ylabel('spec. density')
    ax.autoscale(tight=1)
    return fig


@if_plottable
@save_show
def plot_estimator_power_evo(fel, fig=None, und_duty_factor=1, **kwargs):
    
    if fig == None:
        fig = plt.figure()
    fig.clf()
    
    ax = fig.gca()
    
    P=[]
    E=[]
    z = np.linspace(0,fel.z_sat_min,50)
    for zi in z:
        P.append(np.amax(fel.P(zi)))
        E.append(fel.E(zi))
    
    plt.title('Energy@{:.2f}m = {:.2e} J'.format(fel.z_sat_min / und_duty_factor, E[-1]))
    ax.semilogy(z  / und_duty_factor, np.array(P))
    ax.set_ylabel('P [W]')
    ax.set_xlabel('z [m]')
    return fig


@if_plottable
@save_show
def plot_estimator_energy_evo(fel, fig=None, und_duty_factor=1, **kwargs):
    
    if fig == None:
        fig = plt.figure()
    fig.clf()
    
    ax = fig.gca()
    
    E=[]
    z = np.linspace(0,fel.z_sat_min,50)
    for zi in z:
        E.append(fel.E(zi))
        
    plt.title('Energy@{:.2f}m = {:.2e} J'.format(fel.z_sat_min / und_duty_factor, E[-1]))
    ax.semilogy(z / und_duty_factor, np.array(E))
    ax.set_ylabel('E [J]')
    ax.set_xlabel('z [m]')
    return fig
    
    