"""
user interface for viewing genesis simulation results
"""

"""
MEMO

plt.gcf() to get current figure
plt.gca() to get current axis

ax.set_xlabel('')
ax.get_xlim()
ax.set_xlim([0, 1])
ax.set_ylim(ymin=0)

"""

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

from ocelot.gui.colormaps2d.colormap2d import *
from ocelot.gui.settings_plot import *
from ocelot.gui.dfl_plot import *  # tmp
from ocelot.gui.beam_plot import plot_beam  # tmp

# from pylab import rc, rcParams #tmp

# in order to run decorators properly
import functools

_logger = logging.getLogger(__name__)

__author__ = "Svitozar Serkez"


@if_plottable
def plot_gen_out_all_paral(exp_dir, stage=1, savefig='png', debug=1):
    print('start')
    from ocelot.utils.xfel_utils import background
    i = 0
    dir = exp_dir + 'run_' + str(i) + '/'

    while (os.path.exists(dir)):
        print(i)
        file = dir + 'run.' + str(i) + '.s' + str(stage) + '.gout'
        if (file):
            print('good', i)
            background(
                """plot_gen_out_all(""""+file+"""", choice=(1,1,1,1,0,0,0,0,0,0,0),debug=""" + str(debug) + """)""")

        i += 1
        dir = exp_dir + 'run_' + str(i) + '/'
        print(dir)

    return


# plot_gen_stat(proj_dir=exp_dir, run_inp=[], stage_inp=[], param_inp=[], s_param_inp=['p_int','energy','r_size_weighted'], z_param_inp=[], dfl_param_inp=[], s_inp=['max'], z_inp=[0,'end'], savefig=1, saveval=1, showfig=0, debug=0)

@if_plottable
def plot_gen_out_all(handle=None, savefig='png', showfig=False, choice='all', vartype_dfl=complex128, debug=1, *args,
                     **kwargs):
    """
    plots all possible output from the genesis output
    handle is either:
        genesis output object
        path to genesis output file
        path to folders with genesis output files

    choice=(1,1,1,1,6.05,1,0,0,0,0,0)
            0 1 2 3 4    5 6 7 8 9 10
        0 - electron evolution
        1 - radiation evolution
        2 - profile at z=0m
        3 - profile at the end
        4 - profile every m meters
        5 - dfl at the end, space    -time      domain
        6 -                 inv.space-time      domain
        7 -                 space    -frequency domain
        8 -                 inv.space-frequency domain
        9 - dpa as edist at the end, smeared
        10 - dpa as edist at the end, not smeared
        11 - wigner distribution at end,
        12 - ebeam bucket at max power

    #picks as an input "GenesisOutput" object, file path of directory as strings.
    #plots e-beam evolution, radiation evolution, initial and final simulation window
    #If folder path is provided, all *.gout and *.out files are plotted
    """

    _logger.info('plotting all genesis output')
    _logger.debug('choice = ' + str(choice))
    plotting_time = time.time()

    # plt.ioff()

    if savefig is True:
        savefig = 'png'

    if choice == 'all':
        choice = (1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, -1, 1)
    elif choice == 'gen':
        choice = (1, 1, 1, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0)

    if len(choice) > 13:
        choice = choice[:13]
    elif len(choice) < 13:
        choice += tuple((np.zeros(13 - len(choice)).astype(int)))

    if os.path.isdir(str(handle)):
        handles = []
        for root, dirs, files in os.walk(handle):
            for name in files:
                if name.endswith('.gout') or name.endswith('.out'):
                    handles.append(os.path.join(root, name))
        _logger.info('\n  plotting all files in {}'.format(str(handle)))
    else:
        handles = [handle]

    for handle in handles:

        if os.path.isfile(str(handle)):
            # if os.path.getsize(str(handle)) > 0:
            _logger.info('plotting {}'.format(str(handle)))
            try:
                handle = read_out_file(handle, read_level=2, debug=debug)
            except (IOError, ValueError):
                continue

        if isinstance(handle, GenesisOutput):
            if choice[0]:
                f0 = plot_gen_out_e(handle, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)
            if choice[1]:
                f1 = plot_gen_out_ph(handle, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)
            if choice[2]:
                f2 = plot_gen_out_z(handle, z=0, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)
            if choice[3]:
                f3 = plot_gen_out_z(handle, z=inf, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)
            if choice[11] != 0:
                if choice[11] == -1:
                    try:
                        W = wigner_out(handle, pad=2)
                        plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    except:
                        _logger.warning('could not plot wigner')
                else:
                    if choice[11] == 1:
                        _logger.warning(
                            'choice[11] in plot_gen_out_all defines interval of Wigner plotting. To plot at the end set to "-1"')
                    try:
                        for z in np.arange(0, np.amax(handle.z), choice[11]):
                            W = wigner_out(handle, z=z, pad=2)
                            plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    except:
                        _logger.warning('could not plot wigner')
            if choice[4] != 0 and choice[4] != []:
                for z in np.arange(choice[4], np.amax(handle.z), choice[4]):
                    plot_gen_out_z(handle, z=z, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)
            if choice[12]:
                try:
                    dpa = read_dpa_file_out(handle)
                    plot_dpa_bucket_out(handle, dpa, scatter=0, slice_pos='max_P', repeat=3, showfig=showfig,
                                        savefig=savefig, cmap=def_cmap)
                except:
                    _logger.warning('could not plot particle buckets')

        if os.path.isfile(handle.filePath + '.dfl') and any(choice[5:8]):
            dfl = read_dfl_file_out(handle, debug=debug)
            if dfl.Nz() == 0:
                _logger.warning('empty dfl, skipping')
            else:
                if choice[5]:
                    f5 = plot_dfl(dfl, showfig=showfig, savefig=savefig, debug=debug)
                if choice[6]:
                    f6 = plot_dfl(dfl, domains='tk', auto_zoom=0, showfig=showfig, savefig=savefig, debug=debug)
                if choice[7]:
                    f7 = plot_dfl(dfl, domains='sf', auto_zoom=0, showfig=showfig, savefig=savefig, debug=debug)
                if choice[8]:
                    f8 = plot_dfl(dfl, domains='kf', auto_zoom=0, showfig=showfig, savefig=savefig, debug=debug)

        if os.path.isfile(handle.filePath + '.dpa') and (choice[9] or choice[10]) and handle('itdp') == True:
            dpa = read_dpa_file_out(handle, debug=debug)
            if np.size(dpa.ph) == 0:
                _logger.warning('empty dpa, skipping')
            else:
                if choice[9]:
                    edist = dpa2edist(handle, dpa, num_part=5e4, smear=1, debug=debug)
                    f9 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig, bins=100,
                                    debug=debug)
                if choice[10]:
                    edist = dpa2edist(handle, dpa, num_part=5e4, smear=0, debug=debug)
                    f10 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig,
                                     bins=(100, 100, 100, 100), debug=debug)

    if savefig is not False:
        _logger.info(ind_str + 'plots recorded to *.' + str(savefig) + ' files')

    if showfig:
        _logger.info(ind_str + 'showing plots, close all to proceed')
        plt.show()
    # else:
    #     plt.close('all')

    _logger.info(ind_str + 'total plotting time {:.2f} seconds'.format(time.time() - plotting_time))


@if_plottable
def plot_gen_out_z(g, z=np.inf, params=['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'rad_spec'],
                   figsize=3.5, x_units='um', y_units='ev', legend=False, fig_name=None, savefig=False, showfig=True,
                   debug=1, *args, **kwargs):
    """
    radiation parameters at distance z
    g/out = GenesisOutput() object
    z distance along undulator [m]
    params = parameters of interest:
        'rad_power+el_current' - radiation power and electron beam current
        'el_energy+el_espread+el_bunching' - electron beam energy +/- spread and bunching
        'rad_phase' - phase of radiation
        'rad_spec' - on-axis spectrum
    out_z_params overrides params
    figsize - np.size of figure (unit-less)
    x_units - units of time domain ('um' of 'fs')
    y_units - units of frequency domain ('nm' of 'ev')
    legend - plot legend - tbd
    fig_name - override figure name
    savefig - save figure
    showfig - show figure
    showtext - print text with additional info
    """

    import matplotlib.ticker as ticker

    if showfig is False and savefig is False:
        return

    t_domain = ['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'el_energy+el_espread', 'rad_phase']
    f_domain = ['rad_spec']
    # t_domain_i = list(set(t_domain).intersection(params))
    # f_domain_i = list(set(f_domain).intersection(params))
    # t_domain_n = len(t_domain_i)
    # f_domain_n = len(f_domain_i)
    # add sorting of f_domain to the end params += [params.pop(i)]

    if 'out_z_params' in kwargs:
        params = kwargs['out_z_params']

    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')

    if z == np.inf:
        # print ('Showing profile parameters at the end of undulator')
        z = np.amax(g.z)

    elif z > np.amax(g.z):
        # print ('Z parameter too large, setting to the undulator end')
        z = np.amax(g.z)

    elif z < np.amin(g.z):
        # print ('Z parameter too small, setting to the undulator entrance')
        z = np.amin(g.z)

    zi = np.where(g.z >= z)[0][0]
    z = g.z[zi]

    if os.path.isfile(str(g)):
        g = read_out_file(g, read_level=2)
    # add check for output object

    # if fig_name is None:
    #     if g.fileName() == '':
    #         fig = plt.figure(params_str)
    #         if debug > 0:
    #             print('    plotting ' + params_str)
    #     else:
    #         fig = plt.figure(g.fileName() + '_' + params_str)
    #         if debug > 0:
    #             print('    plotting ' + g.fileName() + '_' + params_str)
    # else:
    #     fig = plt.figure(fig_name)
    #     if debug > 0:
    #         print('    plotting ' + fig_name)
    if fig_name is None:
        if g.fileName() == '':
            fig = plt.figure('Bunch profile at z={:.3f} [m]'.format(z))
        else:
            fig = plt.figure('Bunch profile at z={:.3f} [m]'.format(z) + g.fileName())
    else:
        fig = plt.figure(fig_name)
    _logger.info('plotting bunch profile at z={:.3f} [m]'.format(z))

    if np.size(figsize) == 1:
        figsize = (3 * figsize, (len(params) + 0.5) * figsize)
    fig.set_size_inches(figsize, forward=True)
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    plt.clf()
    fig.subplots_adjust(hspace=0)
    ax = []

    if g('itdp') == False:
        _logger.error('plotting bunch profile is not applicable for steady-state')
        return

    params_t = list(set(params).difference(f_domain))
    params_t = [v for v in params if v in params_t]
    params_f = list(set(params).difference(t_domain))
    params_f = [v for v in params if v in params_f]
    _logger.debug(ind_str + 'params_t: {:}'.format(params_t))
    _logger.debug(ind_str + 'params_f: {:}'.format(params_f))

    for index, param in enumerate(params_t):
        if len(ax) == 0:
            ax.append(fig.add_subplot(len(params), 1, index + 1))
        else:
            ax.append(fig.add_subplot(len(params), 1, index + 1, sharex=ax[0]))

        if param == 'rad_power+el_current':
            subfig_z_power_curr(ax[-1], g, zi=zi, x_units=x_units, legend=legend, *args, **kwargs)
        elif param == 'el_energy+el_espread+el_bunching':
            subfig_z_energy_espread_bunching(ax[-1], g, zi=zi, x_units=x_units, legend=legend, *args, **kwargs)
        elif param == 'el_energy+el_espread':
            subfig_z_energy_espread(ax[-1], g, zi=zi, x_units=x_units, legend=legend, *args, **kwargs)
        elif param == 'rad_phase':
            subfig_z_phase(ax[-1], g, zi=zi, x_units=x_units, legend=legend, *args, **kwargs)
        # elif param == 'el_energy':
        #     subfig_evo_el_energy(ax[-1], g, legend)
        else:
            pass
            # print('! wrong parameter ' + param)

    for axi in ax:
        if axi is not ax[-1]:
            axi.set_xlabel('')
            for label in axi.get_xticklabels():
                label.set_visible(False)
    axt = len(ax)

    for index, param in enumerate(params_f):
        if len(ax) - axt == 0:
            ax.append(fig.add_subplot(len(params), 1, index + axt + 1))
        else:
            ax.append(fig.add_subplot(len(params), 1, index + axt + 1, sharex=ax[-1]))
        if param == 'rad_spec':
            subfig_z_spec(ax[-1], g, zi=zi, y_units=y_units, estimate_ph_sp_dens=True, legend=legend, *args, **kwargs)
        else:
            _logger.warning(ind_str + 'wrong parameter ' + param)
    axf = len(ax) - axt
    # ax[0].set_xlim(g.z[0], g.z[-1])
    # ax[-1].set_xlabel('z [m]')
    if axt != 0 and axf != 0:
        fig.subplots_adjust(top=0.95, bottom=0.2, right=0.8, left=0.15)
    else:
        fig.subplots_adjust(top=0.95, bottom=0.1, right=0.8, left=0.15)

    for axi in ax[axt:]:
        if axt != 0:
            pos1 = axi.get_position()  # get the original position
            pos2 = [pos1.x0 + 0, pos1.y0 - 0.1, pos1.width / 1.0, pos1.height / 1.0]
            axi.set_position(pos2)
            if axi is not ax[-1]:
                axi.set_xlabel('')
                for label in axi.get_xticklabels():
                    label.set_visible(False)

    plt.draw()
    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        fig.savefig(g.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)

    if showfig:
        dir_lst = g.filePath.split(os.path.sep)
        dir = os.path.sep.join(dir_lst[0:-1]) + os.path.sep
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


@if_plottable
def subfig_z_power_curr(ax_curr, g, zi=None, x_units='um', legend=False, *args, **kwargs):
    ax_curr.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_curr.set_xlabel(r's [$\mu$m]')
        x = g.t * speed_of_light * 1.0e-15 * 1e6
    elif x_units == 'fs':
        ax_curr.set_xlabel(r't [fs]')
        x = g.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi is None:
        zi = -1

    ax_curr.plot(x, g.I / 1e3, 'k--')
    ax_curr.set_ylabel(r'I [kA]')
    ax_curr.set_ylim(ymin=0)
    if kwargs.get('showtext', True):
        ax_curr.text(0.02, 0.98, "Q= {:.2f} pC".format(g.beam_charge * 1e12), fontsize=12, horizontalalignment='left',
                     verticalalignment='top', transform=ax_curr.transAxes,
                     color='black')  # horizontalalignment='center', verticalalignment='center',
    ax_curr.grid(kwargs.get('grid', True))

    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(x, g.p_int[:, zi], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(g.p_int[:,zi])!=np.amin(g.p_int[:,zi]):
    #     ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    if kwargs.get('showtext', True):
        if 'n_photons' in dir(g):
            ax_curr.text(0.98, 0.98, "E= {:.2e} J\nN$_{{phot}}$= {:.2e}".format(g.pulse_energy[zi], g.n_photons[zi]),
                         fontsize=12, horizontalalignment='right', verticalalignment='top', transform=ax_curr.transAxes,
                         color='green')  # horizontalalignment='center', verticalalignment='center',
        else:
            ax_curr.text(0.98, 0.98, "E= {:.2e} J".format(g.pulse_energy[zi]), fontsize=12, horizontalalignment='right',
                         verticalalignment='top', transform=ax_curr.transAxes,
                         color='green')  # horizontalalignment='center', verticalalignment='center',

    ax_curr.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())

    ax_power.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_energy_espread_bunching(ax_energy, g, zi=None, x_units='um', legend=False, *args, **kwargs):
    ax_energy.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_energy.set_xlabel(r's [$\mu$m]')
        x = g.t * speed_of_light * 1.0e-15 * 1e6
    elif x_units == 'fs':
        ax_energy.set_xlabel(r't [fs]')
        x = g.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    ax_energy.plot(x, g.el_energy[:, zi] * m_e_GeV, 'b-', x, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV,
                   'r--', x, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(kwargs.get('grid', True))
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(x, g.bunching[:, zi], 'grey', linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.set_ylim(ymin=0)
    ax_bunching.grid(False)

    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')

    ax_bunching.tick_params(axis='y', which='both', colors='grey')
    ax_bunching.yaxis.label.set_color('grey')

    ax_energy.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_energy_espread(ax_energy, g, zi=None, x_units='um', legend=False, *args, **kwargs):
    ax_energy.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_energy.set_xlabel(r's [$\mu$m]')
        x = g.t * speed_of_light * 1.0e-15 * 1e6
    elif x_units == 'fs':
        ax_energy.set_xlabel(r't [fs]')
        x = g.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    ax_energy.plot(x, g.el_energy[:, zi] * m_e_GeV, 'b-', x, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV,
                   'r--', x, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(kwargs.get('grid', True))
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')

    ax_energy.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_phase(ax_phase, g, zi=None, x_units='um', legend=False, *args, **kwargs):
    ax_phase.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_phase.set_xlabel(r's [$\mu$m]')
        x = g.t * speed_of_light * 1.0e-15 * 1e6
    elif x_units == 'fs':
        ax_phase.set_xlabel(r't [fs]')
        x = g.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi is None:
        zi = -1

    if "rewrap" in kwargs:
        _logger.warning(ind_str + '"rewrap" argument is obsolete')

    if 'subfig_z_phase_Ephase' in kwargs:
        Ephase = kwargs['subfig_z_phase_Ephase']
        if not hasattr(g, 'spec'):
            g.calc_spec()
        g.phase_fix(phen=Ephase)  # creates phi_mid_disp attribute wrt. photon energy Ephase

    if hasattr(g, 'phi_mid_disp'):
        phase_disp = g.phi_mid_disp[:, zi]

    else:
        phase_disp = g.phi_mid[:, zi]
    #     phase = unwrap(g.phi_mid[:, zi])
    #     phase_cor = np.arange(g.nSlices) * (maxspectrum_wavelength - g('xlamds')) / g('xlamds') * g('zsep') * 2 * pi
    #     phase_fixed = phase + phase_cor
    #     phase_fixed -= power[maxspower_index, zi]
    #     n = 1
    #     phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
    # else:
    #     phase_fixed = g.phi_mid[:, zi]
    ax_phase.plot(x, phase_disp, 'k-', linewidth=0.5)
    if kwargs.get('showtext', True):
        if hasattr(g, 'phi_mid_disp'):
            _txt = r'(on axis, rewrapped)'
        else:
            _txt = r'(on axis)'
            ax_phase.text(0.98, 0.98, _txt, fontsize=10, horizontalalignment='right', verticalalignment='top',
                          transform=ax_phase.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_phase.set_ylabel(r'$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])
    ax_phase.grid(kwargs.get('grid', True))

    ax_phase.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_phase.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_spec(ax_spectrum, g, zi=None, y_units='ev', estimate_ph_sp_dens=True, legend=False, *args, **kwargs):
    number_ticks = 6
    # n_pad = 1

    mode = kwargs.get('mode', 'mid')
    mode = kwargs.get('spec_mode', mode)

    if zi is None:
        zi = -1

    if hasattr(g, 'spec'):
        if g.spec_mode != mode:
            g.calc_spec(mode=mode)
    else:
        g.calc_spec(mode=mode)
    # TMP
    # g.calc_spec(mode='int')

    if 'spec' not in dir(g):
        g.calc_spec()

    if y_units == 'nm':
        x = g.freq_lamd
    elif y_units in ['ev', 'eV']:
        x = g.freq_ev

    if estimate_ph_sp_dens:
        y_units = 'ev'
        spec = g.spec_phot_density[:, zi]
        # spec = calc_ph_sp_dens(g.spec[:, zi], g.freq_ev, g.n_photons[zi])
    else:
        spec = g.spec[:, zi]

    # power = np.pad(g.p_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')
    # phase = np.pad(g.phi_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')  # not supported by the numpy 1.6.2

    ax_spectrum.plot(x, spec, 'r-')
    if kwargs.get('showtext', True):
        if mode == 'mid':
            ax_spectrum.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right',
                             verticalalignment='top',
                             transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
        else:
            ax_spectrum.text(0.98, 0.98, r'(integrated assuming on-axis phases)', fontsize=10,
                             horizontalalignment='right', verticalalignment='top',
                             transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',

    ax_spectrum.set_ylim(ymin=0)
    ax_spectrum.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_spectrum.get_yaxis().get_major_formatter().set_scientific(True)
    ax_spectrum.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    ax_spectrum.grid(kwargs.get('grid', True))
    if np.amin(x) != np.amax(x):
        ax_spectrum.set_xlim([np.amin(x), np.amax(x)])

    maxspectrum_index = np.argmax(spec)
    # maxspower_index = np.argmax(power[:, zi])
    maxspectrum_value = x[maxspectrum_index]

    spec_width = None
    if np.sum(spec) != 0:
        pos, width, arr = fwhm3(spec)
        if width is not None:
            if arr[0] == arr[-1]:
                dx = abs(x[pos] - x[pos - 1])
            else:
                dx = abs((x[arr[0]] - x[arr[-1]]) / (arr[0] - arr[-1]))
            spec_width = dx * width / x[
                pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

    if spec_width is not None and maxspectrum_value is not None:
        if y_units == 'nm':
            if kwargs.get('showtext', True):
                ax_spectrum.text(0.02, 0.98,
                                 r"$\lambda^{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (
                                     maxspectrum_value * 1e-9, spec_width), fontsize=12, horizontalalignment='left',
                                 verticalalignment='top', transform=ax_spectrum.transAxes,
                                 color='red')  # horizontalalignment='center', verticalalignment='center',
            if estimate_ph_sp_dens:
                ax_spectrum.set_ylabel(r'[$N_{phot}$/nm](estim)')
            else:
                ax_spectrum.set_ylabel(r'P($\lambda$) [a.u.]')
            ax_spectrum.set_xlabel(r'$\lambda$ [nm]')
        elif y_units in ['ev', 'eV']:
            if kwargs.get('showtext', True):
                ax_spectrum.text(0.02, 0.98, r"$E_{ph}^{max}$= %.2f eV " "\n" "$(\Delta E/E)_{fwhm}$= %.2e" % (
                    maxspectrum_value, spec_width), fontsize=12, horizontalalignment='left', verticalalignment='top',
                                 transform=ax_spectrum.transAxes,
                                 color='red')  # horizontalalignment='center', verticalalignment='center',
            if estimate_ph_sp_dens:
                ax_spectrum.set_ylabel(r'[$N_{phot}$/eV] (estim)')
            else:
                ax_spectrum.set_ylabel(r'P($E_{ph}$) [a.u.]')
            ax_spectrum.set_xlabel(r'$E_{photon}$ [eV]')
        # ax_spectrum.text(0.02, 0.98, r"$\lambda_{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (maxspectrum_value, spec_width), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_spectrum.transAxes, color='red')  # horizontalalignment='center', verticalalignment='center',

    ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spectrum.tick_params(axis='y', which='both', colors='r')
    ax_spectrum.yaxis.label.set_color('r')
    ax_spectrum.yaxis.get_offset_text().set_color(ax_spectrum.yaxis.label.get_color())


@if_plottable
def plot_gen_out_z_old(g, figsize=(10, 14), x_units='um', y_units='ev', legend=True, fig_name=None, z=inf,
                       savefig=False, showfig=1, debug=1, **kwargs):
    print('soon this function will be replaced by plot_gen_out_z_new (currently being tested)')
    number_ticks = 6

    if showfig is False and savefig is False:
        return

    if g('itdp') is False:
        print('    plotting bunch profile at ' + str(z) + ' [m]')
        print('!     not applicable for steady-state')
        return

    import matplotlib.ticker as ticker

    if z == inf:
        # print ('Showing profile parameters at the end of undulator')
        z = np.amax(g.z)

    elif z > np.amax(g.z):
        # print ('Z parameter too large, setting to the undulator end')
        z = np.amax(g.z)

    elif z < np.amin(g.z):
        # print ('Z parameter too small, setting to the undulator entrance')
        z = np.amin(g.z)

    zi = np.where(g.z >= z)[0][0]
    z = g.z[zi]

    if debug > 0:
        print('    plotting bunch profile at ' + str(z) + ' [m]')

    font_size = 1
    if fig_name is None:
        if g.fileName() == '':
            fig = plt.figure('Bunch profile at ' + str(z) + 'm')
        else:
            fig = plt.figure('Bunch profile at ' + str(z) + 'm ' + g.fileName())
    else:
        fig = plt.figure(fig_name)
    fig.set_size_inches(figsize, forward=True)
    # plt.rc('axes', grid=True)

    # left, width = 0.1, 0.85

    plt.clf()

    ax_curr = fig.add_subplot(4, 1, 1)
    ax_curr.clear()
    ax_energy = fig.add_subplot(4, 1, 2, sharex=ax_curr)
    ax_energy.clear()
    ax_phase = fig.add_subplot(4, 1, 3, sharex=ax_curr)
    ax_phase.clear()
    ax_spectrum = fig.add_subplot(4, 1, 4)
    ax_spectrum.clear()

    for ax in ax_curr, ax_phase, ax_spectrum, ax_energy:
        if ax != ax_spectrum and ax != ax_phase:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    fig.subplots_adjust(hspace=0)

    s = g.t * speed_of_light * 1.0e-15 * 1e6

    ax_curr.plot(s, g.I / 1e3, 'k--')
    ax_curr.set_ylabel(r'I [kA]')
    ax_curr.set_ylim(ymin=0)
    ax_curr.text(0.02, 0.98, "Q= %.2f pC" % (g.beam_charge * 1e12), fontsize=12, horizontalalignment='left',
                 verticalalignment='top',
                 transform=ax_curr.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_curr.text(0.98, 0.98, "E= %.2e J" % (g.pulse_energy[zi]), fontsize=12, horizontalalignment='right',
                 verticalalignment='top', transform=ax_curr.transAxes,
                 color='green')  # horizontalalignment='center', verticalalignment='center',
    ax_curr.grid(kwargs.get('grid', True))

    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(s, g.p_int[:, zi], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(g.p_int[:,zi])!=np.amin(g.p_int[:,zi]):
    #     ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]

    # ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:, zi] * m_e_GeV, 'b-', s, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV,
                   'r--', s, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(kwargs.get('grid', True))
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(s, g.bunching[:, zi], 'grey', linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.set_ylim(ymin=0)
    ax_bunching.grid(False)

    n_pad = 1
    power = np.pad(g.p_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)],
                   mode='constant')
    phase = np.pad(g.phi_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)],
                   mode='constant')  # not supported by the numpy 1.6.2

    # spectrum = abs(np.fft.fft(np.sqrt(np.array(power)) * np.exp(1.j * np.array(phase)), axis=0))**2 / sqrt(g.nSlices) / (2 * g.leng / g('ncar'))**2 / 1e10
    # e_0 = 1239.8 / g('xlamds') / 1e9
    # g.freq_ev1 = h_eV_s * np.fft.fftfreq(len(spectrum), d=g('zsep') * g('xlamds') * g('ishsty') / speed_of_light) + e_0
    # lamdscale = 1239.8 / g.freq_ev1

    # lamdscale_array = np.swapaxes(np.tile(lamdscale, (g.nZ, 1)), 0, 1)

    # for std calculation
    # spectrum_lamdpos=np.sum(spectrum*lamdscale_array/np.sum(spectrum,axis=0),axis=0)
    # spectrum_lamdwidth=sqrt(np.sum(spectrum*(lamdscale_array-spectrum_lamdpos)**2/np.sum(spectrum,axis=0),axis=0))

    # ax_spectrum.plot(np.fft.fftshift(lamdscale), np.fft.fftshift(spectrum[:,zi]), 'r-')
    ax_spectrum.plot(g.freq_lamd, g.spec[:, zi], 'r-')
    ax_spectrum.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top',
                     transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel(r'P($\lambda$) [a.u.]')
    ax_spectrum.set_xlabel(r'$\lambda$ [nm]')
    ax_spectrum.set_ylim(ymin=0)
    ax_spectrum.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_spectrum.get_yaxis().get_major_formatter().set_scientific(True)
    ax_spectrum.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    ax_spectrum.grid(kwargs.get('grid', True))
    if np.amin(g.freq_lamd) != np.amax(g.freq_lamd):
        ax_spectrum.set_xlim([np.amin(g.freq_lamd), np.amax(g.freq_lamd)])
    ax_phase.set_xlabel(r's [$\mu$m]')

    maxspectrum_index = np.argmax(g.spec[:, zi])
    maxspower_index = np.argmax(power[:, zi])
    maxspectrum_wavelength = g.freq_lamd[maxspectrum_index] * 1e-9

    spectrum_lamdwidth_fwhm = None
    if np.sum(g.spec[:, zi]) != 0:
        pos, width, arr = fwhm3(g.spec[:, zi])
        if width != None:
            if arr[0] == arr[-1]:
                dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos - 1])
            else:
                dlambda = abs((g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]))
            spectrum_lamdwidth_fwhm = dlambda * width / g.freq_lamd[
                pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

    if spectrum_lamdwidth_fwhm is not None and maxspectrum_wavelength is not None:
        ax_spectrum.text(0.02, 0.98, r"$\lambda_{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (
            maxspectrum_wavelength, spectrum_lamdwidth_fwhm), fontsize=12, horizontalalignment='left',
                         verticalalignment='top', transform=ax_spectrum.transAxes,
                         color='red')  # horizontalalignment='center', verticalalignment='center',

    phase = unwrap(g.phi_mid[:, zi])

    phase_cor = np.arange(g.nSlices) * (maxspectrum_wavelength - g('xlamds')) / g('xlamds') * g('zsep') * 2 * pi
    phase_fixed = phase + phase_cor
    phase_fixed -= power[maxspower_index, zi]
    n = 1
    phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
    ax_phase.plot(s, phase_fixed, 'k-', linewidth=0.5)
    ax_phase.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top',
                  transform=ax_phase.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_phase.set_ylabel(r'$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])
    ax_phase.grid(kwargs.get('grid', True))
    # plt.yticks(plt.yticks()[0][0:-1])

    # ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_phase.xaxis.major.locator.set_params(nbins=number_ticks)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_curr.yaxis.major.locator.set_params(nbins=number_ticks)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(s[0], s[-1])

    fig.subplots_adjust(top=0.95, bottom=0.2, right=0.85, left=0.15)

    # fig.set_size_inches((8,8),forward=True)

    pos1 = ax_spectrum.get_position()  # get the original position
    pos2 = [pos1.x0 + 0, pos1.y0 - 0.1, pos1.width / 1.0, pos1.height / 0.9]
    ax_spectrum.set_position(pos2)

    ax_spectrum.tick_params(axis='y', which='both', colors='r')
    ax_spectrum.yaxis.label.set_color('r')
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')

    ax_bunching.tick_params(axis='y', which='both', colors='grey')
    ax_bunching.yaxis.label.set_color('grey')

    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())
    ax_spectrum.yaxis.get_offset_text().set_color(ax_spectrum.yaxis.label.get_color())

    plt.draw()
    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        fig.savefig(g.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')


@if_plottable
def plot_gen_out_e(g, legend=False, figsize=4, fig_name='Electrons', savefig=False, showfig=True, debug=1, *args,
                   **kwargs):
    fig = plot_gen_out_evo(g, params=['und_quad', 'el_size', 'el_energy', 'el_bunching'], figsize=figsize,
                           legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)


@if_plottable
def plot_gen_out_ph(g, legend=False, figsize=4, fig_name='Radiation', savefig=False, showfig=True, debug=1, *args,
                    **kwargs):
    if g('itdp'):
        fig = plot_gen_out_evo(g, params=['rad_pow_en_log', 'rad_pow_en_lin', 'rad_spec_log', 'rad_size'],
                               figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig,
                               debug=debug)
    else:
        fig = plot_gen_out_evo(g, params=['rad_pow_log', 'rad_size'], figsize=figsize, legend=legend, fig_name=fig_name,
                               savefig=savefig, showfig=showfig, debug=debug)


@if_plottable
def plot_gen_out_slip(g, legend=False, figsize=4, fig_name='Slippage', savefig=False, showfig=True, debug=1, *args,
                      **kwargs):
    if g('itdp'):
        fig = plot_gen_out_evo(g, params=['rad_spec_evo_n', 'rad_pow_evo_n'], figsize=figsize, legend=legend,
                               fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug, *args, **kwargs)


@if_plottable
def plot_gen_out_evo(g, params=['und_quad', 'el_size', 'el_pos', 'el_energy', 'el_bunching', 'rad_pow_en_log',
                                'rad_pow_en_lin', 'rad_spec_log', 'rad_size', 'rad_spec_evo_n', 'rad_pow_evo_n'],
                     figsize=4, legend=False, fig_name='', savefig=False, showfig=True, debug=1, *args, **kwargs):
    """
    plots evolution of given parameters from genesis output with undulator length
    """
    import matplotlib.ticker as ticker

    if showfig is False and savefig is False:
        return

    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')

    if os.path.isfile(str(g)):
        g = read_out_file(g, read_level=2)
    # add check for output object
    if fig_name == '':
        if g.fileName() == '':
            fig = plt.figure(params_str)
            _logger.info('plotting ' + params_str)
        else:
            fig = plt.figure(g.fileName() + '_' + params_str)
            _logger.info('plotting ' + g.fileName() + '_' + params_str)
    else:
        fig = plt.figure(fig_name)
        _logger.info('plotting ' + fig_name)

    if np.size(figsize) == 1:
        figsize = (3 * figsize, (len(params) + 0.5) * figsize)

    fig.set_size_inches(figsize, forward=True)
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    plt.clf()
    fig.subplots_adjust(hspace=0)

    ax = []
    is_tdp = g('itdp')
    for index, param in enumerate(params):
        if len(ax) == 0:
            ax.append(fig.add_subplot(len(params), 1, index + 1))
        else:
            ax.append(fig.add_subplot(len(params), 1, index + 1, sharex=ax[0]))
        # ax[-1]
        if param == 'und_quad':
            subfig_evo_und_quad(ax[-1], g, legend, **kwargs)
        elif param == 'und':
            subfig_evo_und(ax[-1], g, legend, **kwargs)
        elif param == 'el_size':
            subfig_evo_el_size(ax[-1], g, legend, **kwargs)
        elif param == 'el_pos':
            subfig_evo_el_pos(ax[-1], g, legend, **kwargs)
        elif param == 'el_energy':
            subfig_evo_el_energy(ax[-1], g, legend, **kwargs)
        elif param == 'el_bunching':
            subfig_evo_el_bunching(ax[-1], g, legend, harm=1, **kwargs)
        elif param.startswith('el_bunching_h'):
            bunch_harm = int(param.replace('el_bunching_h', ''))
            subfig_evo_el_bunching(ax[-1], g, legend, harm=bunch_harm, **kwargs)
        elif param == 'rad_pow_en_log':
            if not is_tdp:
                subfig_evo_rad_pow(ax[-1], g, legend, **kwargs)
            else:
                subfig_evo_rad_pow_en(ax[-1], g, legend, **kwargs)
        elif param == 'rad_pow_en_lin':
            if not is_tdp:
                subfig_evo_rad_pow(ax[-1], g, legend, log=0, **kwargs)
            else:
                subfig_evo_rad_pow_en(ax[-1], g, legend, log=0, **kwargs)
        elif param == 'rad_pow_log':
            subfig_evo_rad_pow(ax[-1], g, legend, **kwargs)
        elif param == 'rad_pow_lin':
            subfig_evo_rad_pow(ax[-1], g, legend, log=0, **kwargs)
        elif param == 'rad_size':
            subfig_rad_size(ax[-1], g, legend, **kwargs)
        elif param == 'rad_spec_log':
            if is_tdp:
                subfig_evo_rad_spec(ax[-1], g, legend, **kwargs)
        elif param == 'rad_spec_nolog':
            if is_tdp:
                subfig_evo_rad_spec(ax[-1], g, legend, log=0, **kwargs)
        elif param == 'rad_spec_evo_n':
            if is_tdp:
                subfig_evo_rad_spec_sz(ax[-1], g, legend, norm=1, **kwargs)
        elif param == 'rad_pow_evo_n':
            if is_tdp:
                subfig_evo_rad_pow_sz(ax[-1], g, legend, norm=1, **kwargs)
        elif param == 'rad_spec_evo':
            if is_tdp:
                subfig_evo_rad_spec_sz(ax[-1], g, legend, norm=0, **kwargs)
        elif param == 'rad_pow_evo':
            if is_tdp:
                subfig_evo_rad_pow_sz(ax[-1], g, legend, norm=0, **kwargs)
        else:
            _logger.warning('wrong parameter ' + param)

    ax[0].set_xlim(g.z[0], g.z[-1])
    ax[-1].set_xlabel('z [m]')
    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.8, left=0.15)

    for axi in ax[0:-1]:
        for label in axi.get_xticklabels():
            label.set_visible(False)

    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        if fig_name == 'Electrons':
            savepath = g.filePath + '_elec.' + str(savefig)
        elif fig_name == 'Radiation':
            savepath = g.filePath + '_rad.' + str(savefig)
        elif fig_name == '':
            savepath = g.filePath + '_' + params_str + '.' + str(savefig)
        else:
            savepath = g.filePath + '_' + fig_name + '.' + str(savefig)
        _logger.debug('saving figure to {}'.format(savepath))
        fig.savefig(savepath, format=savefig)

    plt.draw()
    if showfig is True:
        dir_lst = g.filePath.split(os.path.sep)
        dir = os.path.sep.join(dir_lst[0:-1]) + os.path.sep
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


@if_plottable
def subfig_evo_und_quad(ax_und, g, legend, **kwargs):
    number_ticks = 6

    ax_und.plot(g.z, g.aw, 'b-', linewidth=1.5)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(kwargs.get('grid', True))

    ax_quad = ax_und.twinx()
    ax_quad.plot(g.z, g.qfld, 'r-', linewidth=1.5)
    ax_quad.set_ylabel('Quad')
    ax_quad.grid(False)

    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_quad.yaxis.major.locator.set_params(nbins=number_ticks)

    if np.amax(g.aw) != 0:
        aw_tmp = np.array(g.aw)[np.array(g.aw) != 0]
        if np.amax(aw_tmp) != np.amin(aw_tmp):
            diff = np.amax(aw_tmp) - np.amin(aw_tmp)
            ax_und.set_ylim([np.amin(aw_tmp) - diff / 10, np.amax(aw_tmp) + diff / 10])
    else:
        ax_und.set_ylim([0, 1])
    ax_und.tick_params(axis='y', which='both', colors='b')
    ax_und.yaxis.label.set_color('b')
    ax_quad.tick_params(axis='y', which='both', colors='r')
    ax_quad.yaxis.label.set_color('r')


@if_plottable
def subfig_evo_und(ax_und, g, legend, **kwargs):
    number_ticks = 6

    ax_und.plot(g.z, g.aw, 'b-', linewidth=1.5)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(kwargs.get('grid', True))

    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)

    if np.amax(g.aw) != 0:
        aw_tmp = np.array(g.aw)[np.array(g.aw) != 0]
        if np.amax(aw_tmp) != np.amin(aw_tmp):
            diff = np.amax(aw_tmp) - np.amin(aw_tmp)
            ax_und.set_ylim([np.amin(aw_tmp) - diff / 10, np.amax(aw_tmp) + diff / 10])
    else:
        ax_und.set_ylim([0, 1])
    ax_und.tick_params(axis='y', which='both', colors='b')
    ax_und.yaxis.label.set_color('b')


@if_plottable
def subfig_evo_el_size(ax_size_tsize, g, legend, which='both', **kwargs):
    number_ticks = 6

    if which == 'both' or which == 'averaged':
        ax_size_tsize.plot(g.z, np.average(g.xrms, axis=0, weights=g.I) * 1e6, 'g-', g.z,
                           np.average(g.yrms, axis=0, weights=g.I) * 1e6, 'b-')
    if which == 'both' or which == 'peak_curr':
        idx_pk = np.where(g.I == np.amax(g.I))[0][0]
        ax_size_tsize.plot(g.z, g.xrms[idx_pk, :] * 1e6, 'g--', g.z, g.yrms[idx_pk, :] * 1e6, 'b--')
    ax_size_tsize.set_ylabel(r'$\sigma_{x,y}$ [$\mu$m]')

    ax_size_tsize.set_ylim(ymin=0)
    ax_size_tsize.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_size_tsize.grid(kwargs.get('grid', True))


@if_plottable
def subfig_evo_el_pos(ax_size_tpos, g, legend, which='both', **kwargs):
    if hasattr(g, 'x') and hasattr(g, 'y'):
        if which == 'both' or which == 'averaged':
            ax_size_tpos.plot(g.z, np.average(g.x, axis=0, weights=g.I) * 1e6, 'g-', g.z,
                              np.average(g.y, axis=0, weights=g.I) * 1e6, 'b-')
        if which == 'both' or which == 'peak_curr':
            idx_pk = np.where(g.I == np.amax(g.I))[0][0]
            ax_size_tpos.plot(g.z, g.x[idx_pk, :] * 1e6, 'g--', g.z, g.y[idx_pk, :] * 1e6, 'b--')
        ax_size_tpos.set_ylabel(r'$x,y$ [$\mu$m]')
        ax_size_tpos.grid(kwargs.get('grid', True))


@if_plottable
def subfig_evo_el_energy(ax_energy, g, legend, **kwargs):
    number_ticks = 6

    el_energy = g.el_energy * m_e_MeV
    el_energy_av = int(np.mean(el_energy))
    ax_energy.plot(g.z, np.average(el_energy - el_energy_av, axis=0), 'b-', linewidth=1.5)
    ax_energy.set_ylabel('E + ' + str(el_energy_av) + '[MeV]')
    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.grid(kwargs.get('grid', True))

    ax_spread = ax_energy.twinx()
    ax_spread.plot(g.z, np.average(g.el_e_spread * m_e_GeV * 1000, weights=g.I, axis=0), 'm--', g.z,
                   np.amax(g.el_e_spread * m_e_GeV * 1000, axis=0), 'r--', linewidth=1.5)
    ax_spread.set_ylabel(r'$\sigma_E$ [MeV]')
    ax_spread.grid(False)
    ax_spread.set_ylim(ymin=0)

    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spread.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')
    ax_spread.tick_params(axis='y', which='both', colors='r')
    ax_spread.yaxis.label.set_color('r')


@if_plottable
def subfig_evo_el_bunching(ax_bunching, g, legend, harm=1, **kwargs):
    number_ticks = 6

    if harm == 1:
        ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I, axis=0), 'k-', g.z, np.amax(g.bunching, axis=0),
                         'grey', linewidth=1.5)
    else:
        ax_bunching.plot(g.z, np.average(getattr(g, 'h{}_bunching'.format(harm)), weights=g.I, axis=0), 'k-', g.z,
                         np.amax(getattr(g, 'h{}_bunching'.format(harm)), axis=0), 'grey', linewidth=1.5)
    # ax_bunching.plot(g.z, np.amax(g.bunching, axis=0), 'grey',linewidth=1.5) #only max
    ax_bunching.set_ylabel(r'Bunching h{}'.format(harm))
    ax_bunching.set_ylim(ymin=0)
    # ax_bunching.set_ylim([0,0.8])
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.grid(kwargs.get('grid', True))


@if_plottable
def subfig_evo_rad_pow_en(ax_rad_pow, g, legend, log=1, **kwargs):
    ax_rad_pow.plot(g.z, np.amax(g.p_int, axis=0), 'g-', linewidth=1.5)
    ax_rad_pow.set_ylabel(r'P [W]')
    ax_rad_pow.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_pow.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(g.p_int) > 0 and log:
        ax_rad_pow.set_yscale('log')
    if not log:
        ax_rad_pow.set_ylim(ymin=0)
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_en = ax_rad_pow.twinx()
    ax_rad_en.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_en.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(g.p_int) > 0 and log:
        ax_rad_en.plot(g.z, g.pulse_energy, 'k--', linewidth=1.5)
        ax_rad_en.set_ylabel(r'E [J]')
        ax_rad_en.set_yscale('log')
    if not log:
        if np.amax(g.pulse_energy) < 1e-4:
            ax_rad_en.plot(g.z, g.pulse_energy * 1e6, 'k--', linewidth=1.5)
            ax_rad_en.set_ylabel(r'E [$\mu$J]')
        else:
            ax_rad_en.plot(g.z, g.pulse_energy * 1e3, 'k--', linewidth=1.5)
            ax_rad_en.set_ylabel(r'E [mJ]')
        ax_rad_en.set_ylim(ymin=0)
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_pow.grid(kwargs.get('grid', True))  # , which='minor')
    # ax_rad_pow.grid(False, which="minor")
    ax_rad_pow.tick_params(axis='y', which='both', colors='g')
    ax_rad_pow.yaxis.label.set_color('g')
    ax_rad_en.tick_params(axis='y', which='both', colors='k')
    ax_rad_en.yaxis.label.set_color('k')
    ax_rad_en.grid(False)
    # ax_rad_en.grid(False, which='minor')
    ax_rad_pow.yaxis.get_offset_text().set_color(ax_rad_pow.yaxis.label.get_color())
    ax_rad_en.yaxis.get_offset_text().set_color(ax_rad_en.yaxis.label.get_color())

    if kwargs.get('showtext', True):
        ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W ' '\n' r'$E_{end}$= %.2e J' % (np.amax(g.p_int[:, -1]),
                                                                                       np.mean(g.p_int[:, -1],
                                                                                               axis=0) * g(
                                                                                           'xlamds') * g(
                                                                                           'zsep') * g.nSlices / speed_of_light),
                        fontsize=12, horizontalalignment='right', verticalalignment='bottom',
                        transform=ax_rad_pow.transAxes)


@if_plottable
def subfig_evo_rad_pow(ax_rad_pow, g, legend, log=1, **kwargs):
    ax_rad_pow.plot(g.z, np.amax(g.p_int, axis=0), 'g-', linewidth=1.5)
    ax_rad_pow.set_ylabel('P [W]')
    ax_rad_pow.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_pow.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(g.p_int) > 0 and log:
        ax_rad_pow.set_yscale('log')
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_pow.grid(False)  # , which='minor')
    ax_rad_pow.tick_params(axis='y', which='both', colors='g')
    ax_rad_pow.yaxis.label.set_color('g')
    ax_rad_pow.yaxis.get_offset_text().set_color(ax_rad_pow.yaxis.label.get_color())
    if kwargs.get('showtext', True):
        ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W' % (np.amax(g.p_int[:, -1])), fontsize=12,
                        horizontalalignment='right', verticalalignment='bottom', transform=ax_rad_pow.transAxes)


@if_plottable
def subfig_evo_rad_spec(ax_spectrum, g, legend, log=1, **kwargs):
    if not hasattr(g, 'spec'):
        g.calc_spec()

    ax_spectrum.plot(g.z, np.amax(g.spec, axis=0), 'r-', linewidth=1.5)
    ax_spectrum.text(0.5, 0.98, r"(on axis)", fontsize=10, horizontalalignment='center', verticalalignment='top',
                     transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel(r'P$(\lambda)_{max}$ [a.u.]')
    plt.yticks(plt.yticks()[0][0:-1])

    if np.amax(np.amax(g.spec, axis=0)) > 0 and log:
        ax_spectrum.set_yscale('log')
    ax_spectrum.grid(kwargs.get('grid', True))

    spectrum_lamdwidth_fwhm = np.zeros_like(g.z)
    spectrum_lamdwidth_std = np.zeros_like(g.z)

    for zz in range(g.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(g.spec[:, zz]) != 0:
            pos, width, arr = fwhm3(g.spec[:, zz])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos - 1])
                else:
                    dlambda = abs((g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]))
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / g.freq_lamd[pos]
                # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

            spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz],
                                                                                           0, 1)

        # try:
        #     peak = fwhm3(g.spec[:, zz])
        #     spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[0] - g.freq_lamd[1]) * peak[1] / g.freq_lamd[peak[0]]  # the FWHM of spectral line (error when paekpos is at the edge of lamdscale)
        # except:
        #     spectrum_lamdwidth_fwhm[zz] = 0
        #
        # try:
        #     spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)
        # except:
        #     spectrum_lamdwidth_std[zz] = 0

    ax_spec_bandw = ax_spectrum.twinx()
    ax_spec_bandw.plot(g.z, spectrum_lamdwidth_fwhm * 100, 'm--', label="fwhm")
    ax_spec_bandw.plot(g.z, 2 * spectrum_lamdwidth_std * 100, 'm:', label="std")
    ax_spec_bandw.grid(False)
    plt.yticks(plt.yticks()[0][0:-1])
    ax_spec_bandw.set_ylim(ymin=0)

    if legend:
        ax_spec_bandw.legend()
    ax_spec_bandw.set_ylabel(r'$\Delta\lambda/\lambda, \%$' + '\n' + r'(-- fwhm, $\cdots2\sigma$)')
    if spectrum_lamdwidth_fwhm[-1] != None and kwargs.get('showtext', True):
        ax_spec_bandw.text(0.98, 0.98, r"$(\Delta\lambda/\lambda)_{end}^{fwhm}$= %.2e" % (spectrum_lamdwidth_fwhm[-1]),
                           fontsize=12, horizontalalignment='right', verticalalignment='top',
                           transform=ax_spec_bandw.transAxes)

    ax_spectrum.yaxis.label.set_color('r')
    ax_spectrum.tick_params(axis='y', which='both', colors=ax_spectrum.yaxis.label.get_color())
    ax_spectrum.yaxis.get_offset_text().set_color(ax_spectrum.yaxis.label.get_color())

    ax_spec_bandw.yaxis.label.set_color('m')
    ax_spec_bandw.tick_params(axis='y', which='both', colors=ax_spec_bandw.yaxis.label.get_color())
    ax_spec_bandw.yaxis.get_offset_text().set_color(ax_spec_bandw.yaxis.label.get_color())

    # yticks = ax_spec_bandw.yaxis.get_major_ticks()
    # # print(yticks)
    # print (yticks[1].label.get_text())
    # yticks[1].label.set_text('')
    # print (yticks[1].label.get_text())


@if_plottable
def subfig_rad_size(ax_size_t, g, legend, **kwargs):
    if g.nSlices == 1:
        ax_size_t.plot(g.z, g.r_size.T * 2 * 1e6, 'b-', linewidth=1.5)
        ax_size_t.plot([np.amin(g.z), np.amax(g.z)], [g.leng * 1e6, g.leng * 1e6], 'b-', linewidth=1.0)
        ax_size_t.set_ylabel('transverse $[\mu m]$')
    else:
        if hasattr(g, r'rad_t_size_weighted'):
            ax_size_t.plot(g.z, g.rad_t_size_weighted, 'b-', linewidth=1.5)
        else:
            if np.amax(g.p_int) > 0:
                weight = g.p_int + np.amin(g.p_int[g.p_int != 0]) / 1e6
            else:
                weight = np.ones_like(g.p_int)

            ax_size_t.plot(g.z, np.average(g.r_size * 2 * 1e6, weights=weight, axis=0), 'b-', linewidth=1.5)

    ax_size_t.set_ylim(ymin=0)
    # ax_size_t.set_ylabel(r'$\sim$size$_{transv}$ [$\mu$m]'+'\n'+r'($2\sigma$)')
    ax_size_t.set_ylabel(r'$\sim$size$_{transv}$ [$\mu$m]' + '\n' + u'(\u2014 $2\sigma$)')
    ax_size_t.grid(kwargs.get('grid', True))
    plt.yticks(plt.yticks()[0][0:-1])

    if g.nSlices > 1:
        ax_size_s = ax_size_t.twinx()
        size_long_fwhm = np.zeros_like(g.z)
        size_long_std = np.zeros_like(g.z)
        s = g.t * speed_of_light * 1.0e-15 * 1e6
        delta_s = (s[1] - s[0])
        for zz in range(g.nZ):
            # size_long_fwhm[zz] = fwhm(g.s,g.p_int[:, zz])
            if np.sum(g.p_int[:, zz]) != 0:
                # try:
                _, width, _ = fwhm3(g.p_int[:, zz])
                if width != None:
                    size_long_fwhm[zz] = abs(delta_s) * width
                else:
                    size_long_fwhm[zz] = None
                # except:
                #     size_long_fwhm[zz] = 0

                # try:
                size_long_std[zz] = std_moment(s, g.p_int[:, zz])
                # except:
                #     size_long_std[zz] = 0
            else:
                size_long_fwhm[zz] = None
                size_long_std[zz] = None

        ax_size_s.plot(g.z, size_long_fwhm, color='navy', linestyle='--', linewidth=1.0, label="fwhm")
        ax_size_s.plot(g.z, 2 * size_long_std, color='navy', linestyle=':', linewidth=1.0, label="std")
        ax_size_s.set_ylim(ymin=0)
        ax_size_s.set_ylabel(r'size$_{long}$ [$\mu$m]' + '\n' + r'(-- fwhm, $\cdots2\sigma$)')
        ax_size_s.grid(False)
        plt.yticks(plt.yticks()[0][0:-1])

        ax_size_t.yaxis.label.set_color('b')
        ax_size_t.tick_params(axis='y', which='both', colors=ax_size_t.yaxis.label.get_color())
        ax_size_t.yaxis.get_offset_text().set_color(ax_size_t.yaxis.label.get_color())

        ax_size_s.yaxis.label.set_color('navy')
        ax_size_s.tick_params(axis='y', which='both', colors=ax_size_s.yaxis.label.get_color())
        ax_size_s.yaxis.get_offset_text().set_color(ax_size_s.yaxis.label.get_color())

        if legend:
            ax_size_s.legend()


@if_plottable
def subfig_evo_rad_pow_sz(ax_power_evo, g, legend, norm=1, **kwargs):
    if g.nSlices > 1:

        y_units = kwargs.get('subfig_evo_rad_pow_sz_yunits', 'um')
        dgrid = kwargs.get('subfig_evo_rad_pow_sz_dgrid', True)
        z = g.z
        if y_units in ['um']:
            s = g.s * 1e6
            y_label = 's [$\mu$m]'
        elif y_units in ['fs']:
            s = g.s / speed_of_light * 1e15
            y_label = 't [fs]'
        power = g.p_int
        if norm == 1:
            max_power = np.max(power, 0)[np.newaxis, :]
            max_power[max_power == 0] = 1  # avoid division by zero
            power = power / max_power
            # power[isnan(power)]=0
        ax_power_evo.pcolormesh(z, s, power)
        ax_power_evo.set_xlabel('z [m]')
        ax_power_evo.set_ylabel(y_label)
        ax_power_evo.axis('tight')

        ax_power_evo.grid(dgrid)
    else:
        pass


@if_plottable
def subfig_evo_rad_spec_sz(ax_spectrum_evo, g, legend, norm=1, **kwargs):
    if not hasattr(g, 'spec'):
        g.calc_spec()

    if g.nSlices > 1:
        z = g.z

        y_units = kwargs.get('subfig_evo_rad_spec_sz_yunits', 'ev')
        dgrid = kwargs.get('subfig_evo_rad_spec_sz_dgrid', True)

        if y_units in ['ev', 'eV', 'phen']:
            l = g.freq_ev
            y_label = '$E_{photon}$ [eV]'
        else:
            l = g.freq_lamd
            y_label = '$\lambda$ [nm]'
        spectrum = g.spec
        if norm == 1:
            max_spectrum = np.max(spectrum, 0)[np.newaxis, :]
            max_spectrum[max_spectrum == 0] = 1  # avoid division by zero
            spectrum = spectrum / max_spectrum
            # spectrum[isnan(spectrum)]=0
        ax_spectrum_evo.pcolormesh(z, l, spectrum)
        ax_spectrum_evo.set_xlabel('z [m]')
        ax_spectrum_evo.set_ylabel(y_label)
        ax_spectrum_evo.axis('tight')
        ax_spectrum_evo.grid(dgrid)
    else:
        pass


@if_plottable
def plot_gen_out_scanned_z(g, figsize=(10, 14), legend=True, fig_name=None, z=inf, savefig=False):
    if g('itdp') is True:
        print('    plotting scan at ' + str(z) + ' [m]')
        print('!     Not implemented yet for time dependent, skipping')
        return

    if g('iscan') == 0 and g('scan') == 0:
        print('    plotting scan at ' + str(z) + ' [m]')
        print('!     Not a scan, skipping')
        return

    import matplotlib.ticker as ticker

    if z == inf:
        # print 'Showing profile parameters at the end of undulator'
        z = np.amax(g.z)

    elif z > np.amax(g.z):
        # print 'Z parameter too large, setting to the undulator end'
        z = np.amax(g.z)

    elif z < np.amin(g.z):
        # print 'Z parameter too small, setting to the undulator entrance'
        z = np.amin(g.z)

    zi = np.where(g.z >= z)[0][0]
    z = g.z[zi]

    print('    plotting scan at ' + str(z) + ' [m]')

    font_size = 1
    if fig_name is None:
        if g.fileName() == '':
            fig = plt.figure('Genesis scan at ' + str(z) + 'm')
        else:
            fig = plt.figure('Genesis scan at ' + str(z) + 'm ' + g.fileName())
    else:
        fig = plt.figure(fig_name)
    fig.set_size_inches(figsize, forward=True)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85

    plt.clf()

    ax_curr = fig.add_subplot(2, 1, 1)
    ax_curr.clear()
    ax_energy = fig.add_subplot(2, 1, 2, sharex=ax_curr)
    ax_energy.clear()
    # ax_phase=fig.add_subplot(4, 1, 3,sharex=ax_curr)
    # ax_phase.clear()
    # ax_spectrum=fig.add_subplot(4, 1, 4)
    # ax_spectrum.clear()
    for ax in [ax_curr]:  # , ax_energy: #ax_phase, ax_spectrum,
        for label in ax.get_xticklabels():
            label.set_visible(False)

    fig.subplots_adjust(hspace=0)

    s = g.scv  # scan value is written to current colunm

    ax_curr.plot(s, np.linspace(g('curpeak'), g('curpeak'), len(s)), 'k--')
    ax_curr.set_ylabel(r'I[kA]')

    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(s, g.p_int[:, zi], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim([0, np.amax(g.p_int[:, zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]

    # ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:, zi] * m_e_GeV, 'b-', s, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV,
                   'r--', s, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$\n[GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.get_xaxis().get_major_formatter().set_useOffset(False)
    ax_energy.get_xaxis().get_major_formatter().set_scientific(True)

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(s, g.bunching[:, zi], 'grey', linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.grid(False)

    # ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(s[0], s[-1])

    fig.subplots_adjust(top=0.95, bottom=0.2, right=0.85, left=0.15)

    # fig.set_size_inches((8,8),forward=True)

    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')

    ax_bunching.tick_params(axis='y', which='both', colors='grey')
    ax_bunching.yaxis.label.set_color('grey')

    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())

    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        fig.savefig(g.filePath + '_z_' + str(z) + 'm_scan.' + str(savefig), format=savefig)

    return fig


@if_plottable
def plot_gen_stat(proj_dir, run_inp=[], stage_inp=[], param_inp=[],
                  s_param_inp=['p_int', 'pulse_energy', 'r_size_weighted', 'spec', 'spec_phot_density', 'error'],
                  z_param_inp=['p_int', 'phi_mid_disp', 'spec', 'spec_phot_density', 'bunching', 'wigner'],
                  dfl_param_inp=['dfl_spec'], run_param_inp=['p_int', 'spec', 'spec_phot_density', 'pulse_energy'],
                  s_inp=['max'], z_inp=[0, 'end'], run_s_inp=['max'], run_z_inp=['end'], spec_pad=1, savefig=1,
                  saveval=1, showfig=0, debug=1):
    """
    The routine for plotting the statistical info of many GENESIS runs
    --- Will be rewritten and split in several separate modules ---
    
    proj_dir is the directory path in which \run_xxx folders are located.
    run_inp=[1,2,3] number of runs to be processed, default - all possible up to run 1000
    stage_inp=[1,2,3] stages to be processed, default - all possible up to stage 15
    s_param_inp=['p_int','energy'] parameters to be displayed at certain position along the beam as a function of undulator length
    z_param_inp=['p_int','phi_mid_disp','spec','bunching'] parameters to be displayed at certain position along the undulator length as a function of location along the beam.
    s_inp=[1e-6,'max','mean'] positions at s to be plotted as function of z, max value of s as a function of z, mean value of s as a function of z
    z_inp=[12,'end'] position of z at which radiation and spectrum parameters are plotted
    savefig=1 save figures to given file format into proj_dir/results folder. 1 corresponds to 'png'. accepts other values, such as 'eps'
    saveval=1, saves values being plotted to text files with the same names as the figures. first column - argument value (s[um],z[m],or lamd[nm]), second column - averaged parameters over shots, rest columns - single shot values.
    showfig=1 envokes plt.show() to display figures interactively. May be time- and processor-consuming

    dfl_power, dfl_spec, dfl_size, dfl_divergence
    """
    import copy
    rc('text', usetex=False)
    dict_name = {'p_int': 'radiation power', 'pulse_energy': 'radiation pulse energy',
                 'el_e_spread': 'el.beam energy spread', 'el_energy': 'el.beam energy average',
                 'bunching': 'el.beam bunching', 'spec': 'radiation on-axis spectral density',
                 'spec_phot_density': 'radiation spectral photon density',
                 'dfl_spec': 'total radiation photon spectral density (dfl)', 'r_size': 'radiation transv size',
                 'r_size_weighted': 'radiation transv size (weighted)', 'xrms': 'el.beam x size',
                 'yrms': 'el.beam y size', 'error': 'genesis simulation error', 'p_mid': 'radiation power on-axis',
                 'phi_mid': 'radiation phase on-axis', 'increment': 'radiation power increment'}
    dict_unit = {'p_int': '[W]', 'pulse_energy': '[J]', 'el_e_spread': '(gamma)', 'el_energy': '(gamma)',
                 'bunching': '', 'spec': '[arb.units]', 'spec_phot_density': '(estimation) [ph/eV]',
                 'dfl_spec': '[ph/eV]', 'r_size': '[m]', 'xrms': '[m]', 'yrms': '[m]', 'error': ''}

    figsize = (14, 7)
    figsize = (8, 6)

    _logger.info('statistical postprocessing started')
    start_time = time.time()

    if proj_dir[-1] != '/':
        proj_dir += '/'

    if stage_inp == []:
        stage_range = range(15)  # guess possible stages (0 to 100)
    else:
        stage_range = stage_inp

    for stage in stage_range:  # scan through stages

        outlist = [GenesisOutput() for i in range(1000)]

        if run_inp == []:
            run_range = range(1000)
        else:
            run_range = run_inp

        run_range_good = []

        for irun in run_range:
            out_file = proj_dir + 'run_' + str(irun) + '/run.' + str(irun) + '.s' + str(stage) + '.gout'
            if os.path.isfile(out_file):
                try:
                    outlist[irun] = read_out_file(out_file, read_level=2, debug=1)
                    outlist[irun].calc_spec(npad=spec_pad)
                    run_range_good.append(irun)
                except:
                    print('     could not read ' + out_file)

        run_range = run_range_good

        if len(run_range_good) == 0:
            continue

        nSlices = np.array([int(outlist[run].nSlices) for run in run_range])
        nZ = np.array([int(outlist[run].nZ) for run in run_range])

        f_nSlices = np.argmax(np.bincount(nSlices))  # most abundant number of slices
        f_nZ = np.argmax(np.bincount(nZ))  # most abundant number of z records

        index = np.where(np.logical_and(nSlices == f_nSlices, nZ == f_nZ))[0]
        run_range_good = [run_range[i] for i in index]
        if len(list(set(run_range) - set(run_range_good))) > 0:
            _logger.info('run_range {}'.format(run_range))
            _logger.info('run_range_good {}'.format(run_range_good))
            _logger.info('discarding runs {}'.format(list(set(run_range) - set(run_range_good))))

        run_range = run_range_good

        # if len(run_range)!=0 and debug>0:
        # print('stage = ', stage)

        # check if all gout have the same number of slices nSlice and history records nZ
        # for irun in run_range[1:]:
        #     if outlist[irun].nSlices != outlist[run_range[0]].nSlices or outlist[irun].nZ != outlist[run_range[0]].nZ:
        #         raise ValueError('Non-uniform out objects')

        if run_range == [] or len(run_range) == 1:
            continue

        _logger.info(ind_str + 'processing runs ' + str(run_range) + ' of stage ' + str(stage))

        # for irun in run_range:
        #     out_file=proj_dir+'run_'+str(irun)+'/run.'+str(irun)+'.s'+str(stage)+'.gout'
        #     outlist[irun] = read_out_file(out_file,read_level=1)
        #     print(outlist[irun].sliceKeys)

        # if param_inp==[]:
        #    if debug>1: print(outlist[run_range[0]].sliceKeys_used)
        #    param_range=outlist[run_range[0]].sliceKeys_used
        # else:
        param_range = param_inp

        if savefig != False or saveval != False:
            if savefig == True:
                savefig = 'png'
            saving_path = proj_dir + 'results/'
            if not os.path.isdir(saving_path):
                os.makedirs(saving_path)
            if debug > 1:
                print('      saving to ' + saving_path)

        # if s_param_inp==[]:
        #     s_param_range=param_range
        # else:

        s_param_range = s_param_inp
        if debug > 0:
            print('    processing S parameters ' + str(s_param_range))
        if debug > 1:
            print('      s_inp ' + str(s_inp))

        for param in s_param_range:
            for s_ind in s_inp:
                s_value = []
                s_fig_name = 'stage_' + str(stage) + '__Z__' + dict_name.get(param, param).replace(' ', '_').replace(
                    '.', '_') + '__' + str(s_ind)
                for irun in run_range:
                    if not hasattr(outlist[irun], param):
                        continue
                    else:
                        if debug > 1:
                            print('parameter = %s; s = %s; run = %s;' % (param, s_ind, irun))
                        param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                    if debug > 1:
                        print('shape param_matrix', np.shape(param_matrix))
                    if debug > 1:
                        print('length', len(param_matrix), len(outlist[irun].z))

                    if len(param_matrix) == len(outlist[irun].z):
                        s_value.append(param_matrix)
                    else:
                        if s_ind == 'max':
                            s_value.append(np.amax(param_matrix, axis=0))
                        elif s_ind == 'max_cur':
                            s_value.append(param_matrix[outlist[irun].sn_Imax, :])
                        elif s_ind == 'mean':
                            s_value.append(np.mean(param_matrix, axis=0))
                        else:
                            si = np.where(outlist[irun].s <= s_ind)[-1][-1]
                            s_value.append(param_matrix[si, :])
                if s_value != []:
                    fig = plt.figure(s_fig_name)
                    fig.clf()
                    fig.set_size_inches(figsize, forward=True)
                    if debug > 1:
                        print('plotting array shapes', np.shape(outlist[irun].z), np.shape(np.swapaxes(s_value, 0, 1)))
                    fig = plt.plot(outlist[irun].z, np.swapaxes(s_value, 0, 1), '0.8', linewidth=1)
                    fig = plt.plot(outlist[irun].z, s_value[0], '0.5', linewidth=1)
                    fig = plt.plot(outlist[irun].z, np.mean(s_value, 0), 'k', linewidth=2)
                    plt.xlim([np.min(outlist[irun].z), np.max(outlist[irun].z)])

                    # fig[0].axes.get_yaxis().get_major_formatter().set_scientific(True)
                    # plt.ticklabel_format(style='sci')
                    plt.xlabel('z [m]')
                    plt.ylabel(dict_name.get(param, param) + ' ' + dict_unit.get(param, ''))
                    if savefig != False:
                        if debug > 1:
                            print('      saving ' + s_fig_name + '.' + savefig)
                        plt.draw()
                        plt.savefig(saving_path + s_fig_name + '.' + savefig, format=savefig)
                    if saveval != False:
                        if debug > 1:
                            print('      saving ' + s_fig_name + '.txt')
                        np.savetxt(saving_path + s_fig_name + '.txt',
                                   np.vstack([outlist[irun].z, np.mean(s_value, 0), s_value]).T, fmt="%E", newline='\n',
                                   comments='')
                    if not showfig:
                        plt.close('all')
        # if z_param_inp==[]:
        #     z_param_range=param_range
        # else:
        z_param_range = z_param_inp
        if debug > 0:
            print('    processing Z parameters ' + str(z_param_range))
        if debug > 1:
            print('      z_inp ' + str(z_inp))

        if 'wigner' in z_param_range:
            wig_pad = 2
            if debug > 0:
                print('    processing Wigner')
            for z_ind in z_inp:
                print('      z=', z_ind)
                w = np.zeros((outlist[irun].nSlices * wig_pad, outlist[irun].nSlices * wig_pad))
                for irun in run_range:
                    out = outlist[irun]
                    W = wigner_out(out, z=z_ind, debug=0, pad=wig_pad)
                    w += W.wig
                W.wig = w / len(outlist)

                W.filePath = proj_dir + 'results' + os.path.sep + 'stage_' + str(stage) + '__WIG__' + str(z_ind) + '__m'
                wig_fig_name = 'stage_' + str(stage) + '__WIG__' + str(z_ind) + '__m'
                plot_wigner(W, z=z_ind, x_units='um', y_units='ev', fig_name=wig_fig_name, savefig=savefig,
                            showfig=showfig, debug=0)
                if saveval != False:
                    if debug > 1:
                        print('      saving ' + wig_fig_name + '.txt')
                    np.savetxt(saving_path + wig_fig_name + '.txt', W.wig, fmt='%E ', newline='\n')
                    np.savetxt(saving_path + wig_fig_name + '_sc.txt',
                               np.vstack([speed_of_light * h_eV_s * 1e9 / W.freq_lamd, W.s]).T, fmt='%E ', newline='\n',
                               header=' E[eV], s[m]')

        for param in z_param_range:
            for z_ind in z_inp:
                z_value = []
                z_fig_name = 'stage_' + str(stage) + '__S__' + dict_name.get(param, param).replace(' ', '_').replace(
                    '.', '_') + '__' + str(z_ind) + '__m'
                for irun in run_range:
                    if not hasattr(outlist[irun], param):
                        break
                    else:
                        if debug > 1:
                            print('parameter = %s; z = %s; run = %s;' % (param, z_ind, irun))
                        param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                    if debug > 1:
                        print('shape param_matrix', np.shape(param_matrix))
                    if debug > 1:
                        print('length', len(param_matrix), len(outlist[irun].z))

                    if len(param_matrix) == len(outlist[irun].z):  # case if the array is 1D (no s/z matrix presented)
                        break
                    else:
                        if z_ind == 'end' or z_ind == inf:
                            z_value.append(param_matrix[:, -1])  # after undulator
                        elif z_ind == 'start':
                            z_value.append(param_matrix[:, 0])  # before undulator
                        else:
                            zi = np.where(outlist[irun].z <= z_ind)[-1][-1]
                            z_value.append(param_matrix[:, zi])

                if z_value != []:
                    fig = plt.figure(z_fig_name)
                    fig.clf()
                    fig.set_size_inches(figsize, forward=True)
                    if param in ['spec', 'spec_phot_density']:
                        freq_scale = outlist[irun].freq_ev  # *1e9
                        if debug > 1:
                            print('plotting array shapes freq', np.shape(freq_scale),
                                  np.shape(np.swapaxes(z_value, 0, 1)))
                        fig = plt.plot(freq_scale, np.swapaxes(z_value, 0, 1), '0.8')
                        fig = plt.plot(freq_scale, z_value[0], '0.5', linewidth=1)
                        fig = plt.plot(freq_scale, np.mean(z_value, 0), 'k', linewidth=2)
                        plt.xlim([np.min(freq_scale), np.max(freq_scale)])
                        plt.xlabel('$E_{photon}$ [eV]')
                    else:
                        s_scale = outlist[irun].s * 1e6
                        if debug > 1:
                            print('plotting array shapes', np.shape(s_scale), np.shape(np.swapaxes(z_value, 0, 1)))
                        fig = plt.plot(s_scale, np.swapaxes(z_value, 0, 1), '0.8')
                        fig = plt.plot(s_scale, z_value[0], '0.5', linewidth=1)
                        fig = plt.plot(s_scale, np.mean(z_value, 0), 'k', linewidth=2)
                        plt.xlim([np.min(s_scale), np.max(s_scale)])
                        plt.xlabel('s [um]')
                    plt.ylabel(dict_name.get(param, param) + ' ' + dict_unit.get(param, ''))
                    if savefig != False:
                        if debug > 1:
                            print('      saving ' + z_fig_name + '.' + savefig)
                        plt.draw()
                        plt.savefig(saving_path + z_fig_name + '.' + savefig, format=savefig)
                    if saveval != False:
                        if debug > 1:
                            print('      saving ' + z_fig_name + '.txt')
                        if param in ['spec', 'spec_phot_density']:
                            np.savetxt(saving_path + z_fig_name + '.txt',
                                       np.vstack([outlist[irun].freq_ev, np.mean(z_value, 0), z_value]).T, fmt="%E",
                                       newline='\n', comments='')
                        else:
                            np.savetxt(saving_path + z_fig_name + '.txt',
                                       np.vstack([outlist[irun].s * 1e6, np.mean(z_value, 0), z_value]).T, fmt="%E",
                                       newline='\n', comments='')
                    if not showfig:
                        plt.close('all')
        # if run_param_inp==[]:
        #     run_param_range=[]
        # else:
        run_param_range = run_param_inp
        if debug > 0:
            print('    processing run parameters ' + str(run_param_range))
        if debug > 1:
            print('      run_s_inp ' + str(run_s_inp))
        if debug > 1:
            print('      run_z_inp ' + str(run_z_inp))

        for param in run_param_range:
            for z_ind in run_z_inp:
                for s_ind in run_s_inp:  # not optimal
                    run_value = []
                    run_value_arr = []
                    run_fig_name = 'stage_' + str(stage) + '__RUN__' + dict_name.get(param, param).replace(' ',
                                                                                                           '_').replace(
                        '.', '_') + '__' + str(s_ind) + '__um__' + str(z_ind) + '__m'
                    for irun in run_range:
                        if not hasattr(outlist[irun], param):
                            break
                        else:
                            if debug > 1:
                                print('parameter = %s; z = %s; s = %s; run = %s' % (param, z_ind, s_ind, irun))
                            param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                            if debug > 1:
                                print('shape param_matrix', np.shape(param_matrix))
                            if debug > 1:
                                print('length', len(param_matrix), len(outlist[irun].z))

                        if len(param_matrix) != len(
                                outlist[irun].z):  # case if the array is 1D (no s/z matrix presented)
                            if z_ind == 'end' or z_ind == inf:
                                run_value = param_matrix[:, -1]  # after undulator
                            elif z_ind == 'start':
                                run_value = param_matrix[:, 0]  # before undulator
                            else:
                                zi = np.where(outlist[irun].z <= z_ind)[-1][-1]
                                run_value = param_matrix[:, zi]
                        else:
                            run_value = param_matrix

                        if s_ind == 'max':
                            run_value = np.amax(run_value)
                        elif s_ind == 'max_cur':
                            run_value = run_value[outlist[irun].sn_Imax]
                        elif s_ind == 'mean':
                            run_value = np.amean(run_value)
                        else:
                            si = np.where(outlist[irun].s <= s_ind)[-1][-1]
                            run_value = run_value[si]

                        if debug > 1:
                            print('run_value ', run_value)

                        run_value_arr.append(run_value)

                    if run_value_arr != []:
                        fig = plt.figure(run_fig_name)
                        fig.clf()
                        fig.set_size_inches(figsize, forward=True)

                        fig = plt.plot(run_range, run_value_arr, 'k')
                        plt.xlim([np.min(run_range), np.max(run_range)])
                        plt.xlabel('run')

                        plt.ylabel(dict_name.get(param, param) + '  ' + dict_unit.get(param, '') + ' (' + str(
                            s_ind) + ' um, ' + str(z_ind) + ' m)')
                        if savefig is not False:
                            if debug > 1:
                                print('      saving ' + run_fig_name + '.' + savefig)
                            plt.draw()
                            plt.savefig(saving_path + run_fig_name + '.' + savefig, format=savefig)
                        if saveval != False:
                            if debug > 1:
                                print('      saving ' + run_fig_name + '.txt')
                            np.savetxt(saving_path + run_fig_name + '.txt', np.vstack([run_range, run_value_arr]).T,
                                       fmt="%E", newline='\n', comments='')
                        if not showfig:
                            plt.close('all')

        if dfl_param_inp != []:
            if debug > 0:
                print('    processing DFL parameters ' + str(dfl_param_inp))

        for param in dfl_param_inp:
            dfl_value = []
            dfl_fig_name = 'stage_' + str(stage) + '__DFL__' + param.replace(' ', '_').replace('.', '_') + '__end'
            for irun in run_range:
                dfl_filePath = proj_dir + 'run_' + str(irun) + '/run.' + str(irun) + '.s' + str(stage) + '.gout.dfl'
                dfl = read_dfl_file_out(outlist[irun], debug=debug)
                dfl_pad_z(dfl, spec_pad)
                # dfl=read_dfl_file(dfl_filePath, Nxy=outlist[irun]('ncar'),debug=debug)
                # read_dfl_file(filePath, Nxy=None, Lxy=None, Lz=None, zsep=None, xlamds=None, vartype=complex,debug=1):
                # dfl = dfl.fld
                if dfl.Nz() != 1:
                    freq_scale, spec = dfl.ph_sp_dens()
                    dfl_value.append(spec)
            if dfl_value != []:
                fig = plt.figure(dfl_fig_name)
                fig.clf()
                fig.set_size_inches(figsize, forward=True)
                if param == 'dfl_spec':
                    fig = plt.plot(freq_scale, np.swapaxes(dfl_value, 0, 1), '0.8')
                    fig = plt.plot(freq_scale, dfl_value[0], '0.5', linewidth=1)
                    fig = plt.plot(freq_scale, np.mean(dfl_value, 0), 'k', linewidth=2)
                    plt.xlabel('$E_{photon}$ [eV]')
                plt.ylabel(dict_name.get(param, param) + ' ' + dict_unit.get(param, ''))
                if savefig != False:
                    if debug > 1:
                        print('      saving ' + dfl_fig_name + '.' + savefig)
                    plt.draw()
                    plt.savefig(saving_path + dfl_fig_name + '.' + savefig, format=savefig)
                if saveval != False:
                    if debug > 1:
                        print('      saving ' + dfl_fig_name + '.txt')
                    if param == 'dfl_spec':
                        np.savetxt(saving_path + dfl_fig_name + '.txt',
                                   np.vstack([freq_scale, np.mean(dfl_value, 0), dfl_value]).T, fmt="%E", newline='\n',
                                   comments='')
                if not showfig:
                    plt.close('all')
    if showfig:
        plt.draw()
        plt.show()
    else:
        plt.close('all')

    if debug > 0:
        print('done in %.2f seconds' % (time.time() - start_time))


@if_plottable
def plot_gen_corr(proj_dir, run_inp=[], p1=(), p2=(), savefig=False, showfig=True, saveval=False):
    # param (parameter[str], stage[int], z_position[double], s_position [double or 'max'/'mean' stings])
    # e.g. ('p_int',1,inf,'max') , ('spec',1,inf,'max')

    figsize = (7, 7)

    if proj_dir[-1] != '/':
        proj_dir += '/'

    param_1, stage_1, z_1, s_1 = p1
    param_2, stage_2, z_2, s_2 = p2

    outlist_1 = [GenesisOutput() for i in range(1000)]
    outlist_2 = [GenesisOutput() for i in range(1000)]
    if run_inp == []:
        run_range = range(1000)
    else:
        run_range = run_inp

    run_range_good_1 = []
    run_range_good_2 = []

    if param_1 not in []:
        for irun in run_range:
            out_file_1 = proj_dir + 'run_' + str(irun) + '/run.' + str(irun) + '.s' + str(stage_1) + '.gout'
            if os.path.isfile(out_file_1):
                outlist_1[irun] = read_out_file(out_file_1, read_level=2)
                run_range_good_1.append(irun)

    if param_2 not in []:
        for irun in run_range:
            out_file_2 = proj_dir + 'run_' + str(irun) + '/run.' + str(irun) + '.s' + str(stage_2) + '.gout'
            if os.path.isfile(out_file_2):
                outlist_2[irun] = read_out_file(out_file_2, read_level=2)
                run_range_good_2.append(irun)

    run_range_good = [val for val in run_range_good_1 if val in run_range_good_2]

    if param_1 not in []:
        irun = run_range_good[0]
        if isinstance(s_1, (int, long, float)):
            index_s1 = np.where(outlist_1[irun].s <= s_1)[-1][-1]

        if isinstance(z_1, (int, long, float)):
            index_z1 = np.where(outlist_1[irun].z <= z_1)[-1][-1]

    if param_2 not in []:
        if isinstance(s_2, (int, long, float)):
            index_s2 = np.where(outlist_2[irun].s <= s_2)[-1][-1]

        if isinstance(z_2, (int, long, float)):
            index_z2 = np.where(outlist_2[irun].z <= z_2)[-1][-1]

    matrix_1 = []
    matrix_2 = []

    for i in run_range_good:
        matrix_1.append(getattr(outlist_1[i], param_1))
        matrix_2.append(getattr(outlist_2[i], param_2))

    matrix_1 = np.array(matrix_1)
    matrix_2 = np.array(matrix_2)

    if ndim(matrix_1) == 2:
        var_1 = matrix_1[:, index_z1]
    else:
        if s_1 == 'mean':
            var_1 = np.mean(matrix_1[:, :, index_z1], axis=1)
        elif s_1 == 'max':
            var_1 = np.amax(matrix_1[:, :, index_z1], axis=1)
        else:
            var_1 = matrix_1[:, index_s1, index_z1]

    if ndim(matrix_2) == 2:
        var_2 = matrix_2[:, index_z2]
    else:
        if s_2 == 'mean':
            var_2 = np.mean(matrix_2[:, :, index_z2], axis=1)
        elif s_2 == 'max':
            var_2 = np.amax(matrix_2[:, :, index_z2], axis=1)
        else:
            var_2 = matrix_2[:, index_s2, index_z2]

    corr_fig_name = 'corr_' + param_1 + '_s' + str(stage_1) + '_at' + str(z_1) + '_' + str(
        s_1) + '__' + param_2 + '_s' + str(stage_2) + '_at' + str(z_2) + '_' + str(s_2)

    fig = plt.figure(corr_fig_name)
    fig.clf()
    fig.set_size_inches(figsize, forward=True)
    fig = plt.scatter(var_1, var_2)

    label1 = param_1 + '_s' + str(stage_1) + '_z=' + str(z_1) + '_s=' + str(s_1)
    label2 = param_2 + '_s' + str(stage_2) + '_z=' + str(z_2) + '_s=' + str(s_2)
    label1 = label1.replace('_', ' ')
    label2 = label2.replace('_', ' ')
    plt.xlabel(label1)
    plt.ylabel(label2)

    plt.xlim(np.amin(var_1), np.amax(var_1))
    plt.ylim(np.amin(var_2), np.amax(var_2))

    plt.xlim(0, np.amax(var_1) * 1.05)
    plt.ylim(0, np.amax(var_2) * 1.05)

    saving_path = proj_dir + 'results/'

    plt.draw()
    if savefig is not False:
        print('      saving ' + corr_fig_name + '.' + savefig)
        plt.savefig(saving_path + corr_fig_name + '.' + savefig, format=savefig)
    if saveval is not False:
        print('      saving ' + corr_fig_name + '.txt')
        np.savetxt(saving_path + corr_fig_name + '.txt', np.vstack([var_1, var_2]).T, fmt="%E", newline='\n',
                   comments=param_1 + '_s' + str(stage_1) + '_at' + str(z_1) + '_' + str(
                       s_1) + ' ' + param_2 + '_s' + str(stage_2) + '_at' + str(z_2) + '_' + str(s_2))

    if showfig:
        plt.show()
    else:
        plt.close('all')

    return fig


def plot_dpa_bucket_out(out, dpa=None, slice_pos='max_I', repeat=1, GeV=1, figsize=4, cmap=def_cmap, scatter=True,
                        energy_mean=None, legend=True, fig_name=None, savefig=False, showfig=True, bins=[50, 50],
                        debug=1):
    if dpa is None:
        dpa = read_dpa_file_out(out)

    if out.nSlices > 1:
        if type(slice_pos) == str:
            if slice_pos == 'max_I':
                slice_num = np.argmax(out.I)
            elif slice_pos == 'max_P':
                slice_num = np.argmax(out.power)
            elif slice_pos == 'max_B':
                slice_num = np.argmax(out.bunching[:, -1])
            else:
                raise ValueError('slice_pos text should be "max_I" or "max_P" or "max_B"')
        else:
            if slice_pos < np.amin(out.s) or slice_pos > np.amax(out.s):
                raise ValueError('slice_pos outside out.s range')
            else:
                slice_num = np.where(out.s > slice_pos)[0][0]
        # return plot_dpa_bucket(dpa=dpa, slice_num=slice_num, repeat=repeat, GeV=GeV, figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)
    else:
        slice_num = 0
    slice_pos_act = out.s[slice_num]
    suffix = '_%.2fum_%2.fm' % (slice_pos_act * 1e6, np.amax(out.z))
    if scatter: suffix = '_scatter' + suffix
    return plot_dpa_bucket(dpa=dpa, slice_num=slice_num, repeat=repeat, GeV=GeV, figsize=figsize, cmap=cmap,
                           scatter=scatter, energy_mean=energy_mean, legend=legend, fig_name=fig_name, savefig=savefig,
                           showfig=showfig, suffix=suffix, bins=bins, debug=debug)


@if_plottable
def plot_dpa_bucket(dpa, slice_num=None, repeat=1, GeV=1, figsize=4, cmap=def_cmap, scatter=False, energy_mean=None,
                    legend=True, fig_name=None, savefig=False, showfig=True, suffix='', bins=(50, 50), debug=1,
                    return_mode_gamma=0):
    part_colors = ['darkred', 'orange', 'g', 'b', 'm', 'c', 'y']
    # cmap='BuPu'
    y_bins = bins[0]
    z_bins = bins[1]

    if showfig == False and savefig == False:
        return

    _logger.info('plotting dpa bucket')
    start_time = time.time()

    # if dpa.__class__ != GenesisParticlesDump:
    #     raise ValueError('wrong particle object: should be GenesisParticlesDump')

    if np.shape(dpa.ph)[0] == 1:
        slice_num = 0
    if slice_num is None:
        slice_num = int(np.shape(dpa.ph)[0] / 2)
        _logger.debug(
            ind_str + 'no slice number provided, using middle of the distribution - slice number {}'.format(slice_num))
    else:
        assert (slice_num <= np.shape(dpa.ph)[0]), 'slice_num larger than the dpa shape'

    if fig_name is None:
        fig_name = 'Electron phase space ' + dpa.fileName()
    fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches((5 * figsize, 3 * figsize), forward=True)

    left, width = 0.18, 0.57
    bottom, height = 0.14, 0.55
    left_h = left + width + 0.02 - 0.02
    bottom_h = bottom + height + 0.02 - 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.15, height]

    ax_main = plt.axes(rect_scatter)
    ax_z_hist = plt.axes(rect_histx, sharex=ax_main)
    ax_y_hist = plt.axes(rect_histy, sharey=ax_main)

    # ax_z_hist = plt.subplot2grid((4, 1), (0, 0), rowspan=1)
    # ax_y_hist = plt.subplot2grid((4, 1), (0, 0), rowspan=1)

    # ax_main = plt.subplot2grid((4, 1), (1, 0), rowspan=3, sharex=ax_z_hist)

    nbins = np.shape(dpa.ph)[1]
    phase = deepcopy(dpa.ph[slice_num, :, :])
    energy = deepcopy(dpa.e[slice_num, :, :])
    _logger.debug(ind_str + 'nbins =  {}'.format(nbins))

    if GeV:
        energy *= m_e_MeV
        if energy_mean is None:
            energy_mean = round(np.mean(energy), 0)
    else:
        if energy_mean is None:
            energy_mean = round(np.mean(energy), 1)
    energy -= energy_mean

    phase_flat = phase.flatten()
    energy_flat = energy.flatten()
    for irep in range(repeat - 1):
        phase_flat = np.append(phase_flat, phase.flatten() + 2 * np.pi * (irep + 1))
        energy_flat = np.append(energy_flat, energy.flatten())

    # phase_hist = np.ravel(phase)
    # for irep in range(repeat-1):
    #     phase_hist = np.concatenate((phase_hist, np.ravel(phase) + 2 * np.pi * (irep+1)))

    # hist, edges = np.histogram(phase_hist, bins=50 * repeat)  # calculate current histogram
    hist_z, edges_z = np.histogram(phase_flat, bins=z_bins * repeat)  # calculate current histogram
    edges_z = edges_z[0:-1]  # remove the last bin edge to save equal number of points
    ax_z_hist.bar(edges_z, hist_z, width=edges_z[1] - edges_z[0], color='silver')
    ax_z_hist.set_ylabel('counts')

    hist_y, edges_y = np.histogram(energy_flat, bins=y_bins)  # calculate current histogram
    edges_y = edges_y[0:-1]  # remove the last bin edge to save equal number of points
    ax_y_hist.barh(edges_y, hist_y, height=edges_y[1] - edges_y[0], color='silver')
    ax_y_hist.set_xlabel('counts')

    for label in ax_z_hist.get_xticklabels():
        label.set_visible(False)

    for label in ax_y_hist.get_yticklabels():
        label.set_visible(False)

    if scatter is True:
        for irep in range(repeat):
            for ibin in range(nbins):
                ax_main.scatter(phase[ibin, :] + 2 * np.pi * (irep), energy[ibin, :], color=part_colors[ibin],
                                marker='.')

        # ax_z_hist.set_xlim([edges[0], edges[-1]])

    elif scatter is False:
        ax_main.hist2d(phase_flat, energy_flat, bins=[z_bins * repeat, y_bins], cmin=0, cmap=cmap)

    ax_main.set_xlabel('$\phi$ [rad]')
    if GeV:
        ax_main.set_ylabel('E [MeV] + ' + str(energy_mean / 1000) + ' [GeV]')
    else:
        ax_main.set_ylabel('$\gamma$ + ' + str(energy_mean))

    plt.draw()
    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving to {}'.format(dpa.fileName() + suffix + '.' + savefig))
        plt.savefig(dpa.filePath + suffix + '.' + savefig, format=savefig)

    if showfig:
        rcParams["savefig.directory"] = os.path.dirname(dpa.filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


@if_plottable
def plot_edist(edist, figsize=4, fig_name=None, savefig=False, showfig=True, scatter=False, plot_x_y=True,
               plot_xy_s=True, bins=(50, 50, 50, 50), flip_t=False, x_units='um', y_units='ev', cmin=0, y_offset=None,
               cmap=def_cmap, debug=1):
    if showfig is False and savefig is False:
        return
    _logger.info('plotting edist file')
    start_time = time.time()
    # suffix=''
    # if edist.__class__ != GenesisElectronDist:
    #     raise ValueError('wrong distribution object: should be GenesisElectronDist')

    if np.size(bins) == 1:
        bins = (bins, bins, bins, bins)  # x,y,t,e

    if fig_name is None:
        fig_name = 'Electron distribution ' + edist.fileName()
    fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches(((3 + plot_x_y + plot_xy_s) * figsize, 3 * figsize), forward=True)

    if x_units == 'fs':
        mult = 1e15
        s_label = 't [fs]'
    elif x_units == 'um':
        mult = speed_of_light * 1e6
        s_label = 's [$\mu$m]'
    s = edist.t * mult
    # if flip_t:
    #     s = -edist.t * speed_of_light * 1e6
    # else:
    #     s = edist.t * speed_of_light * 1e6

    hist, edges = np.histogram(s, bins=bins[2])  # calculate current histogram
    edges = edges[0:-1]  # remove the last bin edge to save equal number of points
    hist_int = np.trapz(hist, edges) / mult  # normalize
    hist = np.rint(hist.astype(float) / (hist_int / float(edist.charge())))

    ax_curr = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 1)
    # ax_curr.hist(s, bins,color='b')
    ax_curr.plot(edges, hist / 1000, color='k', linewidth=2)
    ax_curr.set_xlabel(s_label)
    ax_curr.set_ylabel('I [kA]')

    ax_se = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 2 + plot_x_y + plot_xy_s, sharex=ax_curr)
    if y_units in ['ev', 'eV']:
        energy = edist.g * m_e_MeV
    else:  # elif beam_E_plot=='gamma':
        energy = edist.g

    if y_offset is None:
        y_offset = int(np.mean(energy))
    if scatter:
        ax_se.scatter(s, energy - y_offset, marker='.')
    else:
        ax_se.hist2d(s, energy - y_offset, [bins[2], bins[3]], cmin=cmin, cmap=cmap)

    ax_se.set_xlabel(s_label)
    if y_units in ['ev', 'eV']:
        ax_se.set_ylabel('E + ' + str(y_offset) + ' [MeV]')
    else:  # elif beam_E_plot=='gamma':
        ax_se.set_ylabel('$\gamma$ + ' + str(y_offset))

    if plot_xy_s:
        ax_xs = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 4 + plot_x_y, sharex=ax_curr)
        if scatter:
            ax_xs.scatter(s, 1e6 * edist.x, marker='.')
        else:
            ax_xs.hist2d(s, 1e6 * edist.x, [bins[2], bins[0]], cmin=cmin, cmap=cmap)
        ax_xs.set_xlabel(s_label)
        ax_xs.set_ylabel('x [$\mu$m]')

        ax_ys = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 2, sharex=ax_curr)
        if scatter:
            ax_ys.scatter(s, 1e6 * edist.y, marker='.')
        else:
            ax_ys.hist2d(s, 1e6 * edist.y, [bins[2], bins[1]], cmin=cmin, cmap=cmap)
        ax_ys.set_xlabel(s_label)
        ax_ys.set_ylabel('y [$\mu$m]')

    if plot_x_y:
        ax_xy = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 2 + plot_xy_s)
        if scatter:
            ax_xy.scatter(edist.x * 1e6, edist.y * 1e6, marker='.')
        else:
            ax_xy.hist2d(edist.x * 1e6, edist.y * 1e6, [bins[0], bins[1]], cmin=cmin, cmap=cmap)
        ax_xy.set_xlabel('x [$\mu$m]')
        ax_xy.set_ylabel('y [$\mu$m]')

        ax_pxpy = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 4 + 2 * plot_xy_s)
        if scatter:
            ax_pxpy.scatter(edist.xp * 1e6, edist.yp * 1e6, marker='.')
        else:
            ax_pxpy.hist2d(edist.xp * 1e6, edist.yp * 1e6, [bins[0], bins[1]], cmin=cmin, cmap=cmap)
        ax_pxpy.set_xlabel('xp [$\mu$rad]')
        ax_pxpy.set_ylabel('yp [$\mu$rad]')

    # if scatter:
    if flip_t:
        ax_curr.set_xlim([np.amax(s), np.amin(s)])
    else:
        ax_curr.set_xlim([np.amin(s), np.amax(s)])

    ax_curr.set_ylim(ymin=0)

    fig.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.draw()
    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        if savefig in ['png', 'eps', 'pdf', 'jpeg']:
            savepath = edist.filePath + '.' + savefig
        else:
            savepath = savefig
        _logger.debug(ind_str + 'saving to {}'.format(savepath))
        plt.savefig(savepath, format=savepath.split('.')[-1])

    if showfig:
        rcParams["savefig.directory"] = os.path.dirname(edist.filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)

    _logger.info(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
    # return fig


"""
tmp for HXRSS
"""


@if_plottable
def read_plot_dump_proj(exp_dir, stage, run_ids, plot_phase=1, showfig=True, savefig=0, debug=1):
    if showfig == 0 and savefig == 0:
        return None

    t_l_int_arr = []
    t_l_pha_arr = []
    f_l_int_arr = []
    for run_id in run_ids:
        array = np.loadtxt(exp_dir + 'run_' + str(run_id) + '/run.' + str(run_id) + '.s' + str(stage) + '.dfl.t.txt',
                           skiprows=1)
        array = np.rollaxis(array, 1)
        t_l_scale, t_l_int_a, t_l_pha_a = array[0], array[1], array[2]

        array = np.loadtxt(exp_dir + 'run_' + str(run_id) + '/run.' + str(run_id) + '.s' + str(stage) + '.dfl.f.txt',
                           skiprows=1)
        array = np.rollaxis(array, 1)
        f_l_scale, f_l_int_a, f_l_ftlt_abs, f_l_ftlt_ang = array[0], array[1], array[2], array[3]

        t_l_int_arr.append(t_l_int_a)
        t_l_pha_arr.append(t_l_pha_a)
        f_l_int_arr.append(f_l_int_a)

    t_l_scale *= 1e6
    t_l_int_arr = np.array(t_l_int_arr)
    t_l_pha_arr = np.array(t_l_pha_arr)
    f_l_int_arr = np.array(f_l_int_arr)
    t_l_int_arr = np.rollaxis(t_l_int_arr, 1)
    t_l_pha_arr = np.rollaxis(t_l_pha_arr, 1)
    f_l_int_arr = np.rollaxis(f_l_int_arr, 1)

    t_l_pha_arr = np.unwrap(t_l_pha_arr, axis=0)

    if len(run_ids) > 1:
        t_l_int_mean = np.mean(t_l_int_arr, axis=1)
        t_l_pha_mean = np.mean(t_l_pha_arr, axis=1)
        f_l_int_mean = np.mean(f_l_int_arr, axis=1)
    else:
        t_l_int_mean = t_l_int_arr[:, 0]
        t_l_pha_mean = t_l_pha_arr[:, 0]
        f_l_int_mean = f_l_int_arr[:, 0]
    # t_domain,t_norm=plt.figure('t_domain_filtered')
    fig_name = 'stage_' + str(stage) + '__FILT__power'
    t_domain = plt.figure(fig_name, figsize=(15, 7))
    ax1 = t_domain.add_subplot(2 + plot_phase, 1, 1)
    pulse_average_pos = np.sum(t_l_scale * t_l_int_mean) / np.sum(t_l_int_mean)
    ax1.plot(t_l_scale - pulse_average_pos, t_l_int_arr, '0.5')
    ax1.plot(t_l_scale - pulse_average_pos, t_l_int_mean, 'k', linewidth=1.5)
    ax1.plot([0, 0], [0, np.max(t_l_int_arr)], 'r')
    ax1.grid(True)
    plt.ylabel(r'$P$ [W]')

    ax2 = t_domain.add_subplot(2 + plot_phase, 1, 2, sharex=ax1)
    ax2.semilogy(t_l_scale - pulse_average_pos, t_l_int_arr, '0.5')
    ax2.semilogy(t_l_scale - pulse_average_pos, t_l_int_mean, 'k', linewidth=1.5)
    ax2.plot([0, 0], [np.min(t_l_int_arr), np.max(t_l_int_arr)], 'r')
    ax2.grid(True)
    plt.ylabel(r'$P$ [W]')

    if plot_phase:
        ax3 = t_domain.add_subplot(2 + plot_phase, 1, 3, sharex=ax1)
        ax3.plot(t_l_scale - pulse_average_pos, t_l_pha_arr, '0.5')
        ax3.plot(t_l_scale - pulse_average_pos, t_l_pha_mean, 'k', linewidth=1.5)
        plt.ylabel(r'$\phi [rad]$')
    plt.xlabel(r'$S [\mu m]$')

    plt.draw()

    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        if debug > 1:
            print('      saving ' + fig_name + '.' + savefig)
        plt.savefig(exp_dir + 'results/' + fig_name + '.' + savefig, format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')

    fig_name = 'stage_' + str(stage) + '__FILT__spectrum'
    f_domain = plt.figure(fig_name, figsize=(15, 7))
    ax1 = f_domain.add_subplot(2, 1, 1)
    ax1.plot(h_eV_s * speed_of_light / f_l_scale, f_l_int_arr, '0.5')
    ax1.plot(h_eV_s * speed_of_light / f_l_scale, f_l_int_mean, 'k', linewidth=1.5)
    ax1.grid(True)

    plt.ylabel(r'$P(\lambda)$ [a.u.]')
    ax2 = f_domain.add_subplot(2, 1, 2, sharex=ax1)
    # plt.xlabel(r'$S [\mu m]$')

    ax2.plot(h_eV_s * speed_of_light / f_l_scale, f_l_ftlt_abs ** 2, 'r')
    plt.xlabel(r'$E$ [eV]')
    plt.ylabel(r'Filter (ampl$^2$ & phase )')
    ax2_phase = ax2.twinx()
    ax2_phase.plot(h_eV_s * speed_of_light / f_l_scale, f_l_ftlt_ang, 'r--')
    ax2.grid(True)

    # plt.ylabel(r'$Transm$')
    # ax[1].xlabel(r'$E$ [eV]')
    # ax[0].xlabel(r'$P(\lambda)$ [a.u.]')
    # ax[1].xlabel(r'$abs(TrF)$')

    plt.draw()

    if savefig is not False:
        if savefig is True:
            savefig = 'png'
        if debug > 1:
            print('      saving ' + fig_name + '.' + savefig)
        plt.savefig(exp_dir + 'results/' + fig_name + '.' + savefig, format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')


"""
    scheduled for removal
"""


def show_output(g, show_field=False, show_slice=0):
    print('plotting slice', show_slice)

    h = 4.135667516e-15
    c = 299792458.0
    xrms = np.array(g.sliceValues[g.sliceValues.keys()[show_slice]]['xrms'])
    yrms = np.array(g.sliceValues[g.sliceValues.keys()[show_slice]]['yrms'])

    f = plt.figure()
    f.add_subplot(131), plt.plot(g.z, xrms, lw=3), plt.plot(g.z, yrms, lw=3), plt.grid(True)
    f.add_subplot(132), plt.plot(g.z, g.power_z, lw=3), plt.grid(True)
    t = 1.0e+15 * float(g('zsep')) * float(g('xlamds')) * np.arange(0, len(g.I)) / c

    f.add_subplot(133)
    plt.plot(g.t, g.power_int, lw=3)
    plt.plot(t, g.I * np.max(g.power_int) / np.max(g.I), lw=3)
    plt.grid(True)

    npoints = g('ncar')
    zstop = g('zstop')
    delz = g('delz')
    xlamd = g('xlamd')
    xlamds = g('xlamds')
    nslice = g('nslice')
    zsep = g('zsep')
    dgrid = g('dgrid')

    smax = nslice * zsep * xlamds

    print('wavelength ', xlamds)

    if show_field:
        # from mpi4py import MPI

        # comm = MPI.COMM_WORLD
        # slices = readRadiationFile_mpi(comm=comm, fileName=file+'.dfl', npoints=npoints)
        slices = readRadiationFile(fileName=g.path + '.dfl', npoints=npoints)
        print('slices:', slices.shape)

        E = np.zeros_like(slices[0, :, :])
        for i in range(slices.shape[0]):
            E += np.multiply(slices[i, :, :], slices[i, :, :].conjugate())

        fig = plt.figure()
        fig.add_subplot(131)
        m = plt.imshow(abs(E), cmap='YlOrRd')
        z = abs(slices[100, :, :])

        fig.add_subplot(132)
        P = np.zeros_like(slices[:, 0, 0])
        for i in range(len(P)):
            s = sum(np.abs(np.multiply(slices[i, :, :], slices[i, :, :])))
            P[i] = abs(s * s.conjugate()) * (dgrid ** 2 / npoints) ** 2

        t = 1.0e+15 * float(g('zsep')) * float(g('xlamds')) * np.arange(0, len(P)) / c
        plt.plot(t, P)
        plt.title('Pulse/axis')

        fig.add_subplot(133)
        spec = np.abs(np.fft.fft(slices[:, int(npoints / 2), int(npoints / 2)])) ** 2
        freq_ev = h * fftfreq(len(spec), d=zsep * xlamds / c)
        plt.plot(freq_ev, spec)
        plt.title('Spectrum/axis')


def show_plots(displays, fig):
    """
    putting arbitrarily many plots on single figure
    """
    n1 = (len(displays) - 1) / 2 + 1
    n2 = (len(displays) - 1) / n1 + 1
    # print n1, n2
    fmt = str(n1) + str(n2)
    print(fmt)

    for i in range(len(displays)):
        ax = fig.add_subplot(fmt + str(i + 1))
        ax.grid(True)
        for f in displays[i].data:
            x, y = f(x=np.linspace(-10, 10, 100))
            ax.plot(x, y, '.')

    show()


class Display:

    def __init__(self, data=lambda x: (x, 0 * x), xlabel='', ylabel=''):
        self.data = (data,)
        self.xlabel = xlabel
        self.ylabel = ylabel
