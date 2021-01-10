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

# check if Xserver is connected
# havedisplay = "DISPLAY" in os.environ
# if not havedisplay:
# # re-check
# exitval = os.system('python -c "import matplotlib.pyplot as plt; plt.figure()"')
# havedisplay = (exitval == 0)
# if not havedisplay:
# # force matplotlib not ot use Xwindows backend. plots may still be plotted into e.g. *.png
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import logging
from ocelot.adaptors.genesis4 import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.common.math_op import *  # import of mathematical functions
from ocelot.utils.xfel_utils import *
from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.optics.wave import *
from copy import deepcopy
from ocelot.gui.settings_plot import *
from ocelot.gui.dfl_plot import plot_dfl, plot_wigner
from ocelot.gui.genesis_plot import plot_edist

# in order to run decorators properly
import functools

_logger = logging.getLogger(__name__)


# def plot_gen_out_all_paral(exp_dir, stage=1, savefig='png', debug=1):
#     print('start')
#     from ocelot.utils.xfel_utils import background
#     i = 0
#     dir = exp_dir + 'run_' + str(i) + '/'
#
#     while(os.path.exists(dir)):
#         print(i)
#         file = dir + 'run.' + str(i) + '.s'+str(stage)+'.gout'
#         if(file):
#             print('good',i)
#             background("""plot_gen_out_all(""""+file+"""", choice=(1,1,1,1,0,0,0,0,0,0,0),debug="""+str(debug)+""")""")
#
#         i += 1
#         dir = exp_dir + 'run_' + str(i) + '/'
#         print(dir)
#
#     return

@if_plottable
def plot_gen4_out_all(handle=None, savefig='png', showfig=False, choice=(1, 1, 1, 1, 10, 1, 0, 0, 0, 0, 0, 10, 1),
                      vartype_dfl=complex128, *args, **kwargs):
    debug = 1
    """
    plots all possible output from the genesis output
    handle is either:
        genesis output object
        path to genesis output file
        path to folders with genesis output files
    choice=(1,1,1,1,[],1,0,0,0,0,0, 0, 0)
            0 1 2 3 4  5 6 7 8 9 10,11,12
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
        11 - wigner distribution every m meters,
        12 - ebeam bucket at max power
    picks as an input "GenesisOutput" object, file path of directory as strings.
    plots e-beam evolution, radiation evolution, initial and final simulation window
    If folder path is provided, all *.gout and *.out files are plotted
    """
    
    # if debug > 0:
    #     print('  plotting genesis output')
    plotting_time = time.time()

    # # plt.ioff()

    if savefig == True:
        savefig = 'png'

    if choice == 'all':
        choice = (1, 1, 1, 1, 10, 1, 1, 1, 1, 1, 0, 10, 0)
    elif choice == 'gen':
        choice = (1, 1, 1, 1, 10, 0, 0, 0, 0, 0, 0, 0, 0)

    if len(choice) > 13:
        choice = choice[:13]
    elif len(choice) < 13:
        choice += tuple((np.zeros(13 - len(choice)).astype(int)))

    if os.path.isdir(str(handle)):
        handles = []
        for root, dirs, files in os.walk(handle):
            for name in files:
                if name.endswith('out.h5'):
                    handles.append(os.path.join(root, name))
        _logger.info('\n  plotting all files in {}'.format(str(handle)))
    else:
        handles = [handle]

    for handle in handles:

        if os.path.isfile(str(handle)):
            _logger.info('plotting ' + str(handle))
            try:
                handle = read_gout4(handle)
            except (IOError, ValueError):
                continue

        if isinstance(handle, Genesis4Output):
            if choice[0]:
                f0 = plot_gen4_out_e(handle, showfig=showfig, savefig=savefig, debug=debug)
            if choice[1]:
                f1 = plot_gen4_out_ph(handle, showfig=showfig, savefig=savefig, debug=debug)
            if choice[2]:
                f2 = plot_gen4_out_z(handle, z=0, showfig=showfig, savefig=savefig, debug=debug)
            if choice[3]:
                f3 = plot_gen4_out_z(handle, z=np.inf, showfig=showfig, savefig=savefig, debug=debug)
            if choice[11] != 0:
                if choice[11] == -1:
                    # try:
                    W = wigner_out(handle, pad=2)
                    plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    # except:
                    # _logger.warning('could not plot wigner')
                else:
                    if choice[11] == 1:
                        _logger.warning(
                            'choice[11] in plot_gen_out_all defines interval of Wigner plotting. To plot at the end set to "-1"')
                    # try:
                    for z in np.arange(0, np.amax(handle.z), choice[11]):
                        W = wigner_out(handle, z=z, pad=2)
                        plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    W = wigner_out(handle, z=np.inf, pad=2)
                    plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    # except:
                    # _logger.warning('could not plot wigner')
            # if choice[4] != 0:
            # for z in np.arange(choice[4], np.amax(handle.z), choice[4]):
            # plot_gen4_out_z(handle, z=z, showfig=showfig, savefig=savefig, debug=debug)
            if choice[4] != 0 and choice[4] != []:
                for z in np.arange(choice[4], np.amax(handle.z), choice[4]):
                    plot_gen4_out_z(handle, z=z, showfig=showfig, savefig=savefig, debug=debug, *args, **kwargs)

            if choice[12]:
                pass
                # try:
                # plot_dpa_bucket_out(handle,scatter=0,slice_pos='max_P',repeat=3, showfig=showfig,
                # savefig=savefig, cmap=def_cmap)
                # except IOError:
                # pass

        if os.path.isfile(handle.filePath.replace('.out.h5', '.fld.h5')) and any(choice[5:8]):
            dfl = read_dfl4(handle.filePath.replace('.out.h5', '.fld.h5'))
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

        if os.path.isfile(handle.filePath.replace('.out.h5', '.par.h5')) and (choice[9] or choice[10]):
            dpa = read_dpa4(handle.filePath.replace('.out.h5', '.par.h5'))
            if choice[9]:
                try:
                    edist = dpa42edist(dpa, n_part=5e4, fill_gaps=1)
                    f9 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig, bins=100,
                                    debug=debug)
                except:
                    _logger.warning('could not plot smeared edist')
            if choice[10]:
                try:
                    edist = dpa42edist(dpa, n_part=5e4, fill_gaps=0)
                    f10 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig,
                                     bins=(50, 50, 300, 300), debug=debug)
                except:
                    _logger.warning('could not plot unsmeared edist')

    if savefig != False:
        if debug > 0:
            _logger.info('{}plots recorded to *. {} files'.format(ind_str, savefig))

    if showfig:
        if debug > 0:
            _logger.info(ind_str + 'showing plots, close all to proceed')
        plt.show()
    # else:
    # plt.close('all')

    _logger.info(ind_str + 'total plotting time {:.2f} seconds'.format(time.time() - plotting_time))


@if_plottable
def plot_gen4_out_z(out, z=np.inf, params=['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'rad_spec'],
                    figsize=3, x_units='um', y_units='ev', legend=False, fig_name=None, savefig=False, showfig=True,
                    debug=1, *args, **kwargs):
    """
    radiation parameters at distance z
    out/out = GenesisOutput() object
    z distance along undulator [m]
    params = parameters of interest:
        'rad_power+el_current' - radiation power and electron beam current
        'el_energy+el_espread+el_bunching' - electron beam energy +/- spread and bunching
        'rad_phase' - phase of radiation
        'rad_spec' - on-axis spectrum
    figsize - np.size of figure (unit-less)
    x_units - units of time domain ('um' of 'fs')
    y_units - units of frequency domain ('nm' of 'ev')
    legend - plot legend - tbd
    fig_name - override figure name
    savefig - save figure
    showfig - show figure
    """
    import matplotlib.ticker as ticker

    if showfig == False and savefig == False:
        return

    t_domain = ['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'el_energy+el_espread', 'rad_phase']
    f_domain = ['rad_spec']
    # t_domain_i = list(set(t_domain).intersection(params))
    # f_domain_i = list(set(f_domain).intersection(params))
    # t_domain_n = len(t_domain_i)
    # f_domain_n = len(f_domain_i)
    # add sorting of f_domain to the end params += [params.pop(i)]
    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')

    if z == np.inf:
        # print ('Showing profile parameters at the end of undulator')
        z = np.amax(out.z)

    elif z > np.amax(out.z):
        # print ('Z parameter too large, setting to the undulator end')
        z = np.amax(out.z)

    elif z < np.amin(out.z):
        # print ('Z parameter too small, setting to the undulator entrance')
        z = np.amin(out.z)

    zi = np.where(out.z >= z)[0][0]
    z = out.z[zi]

    # zi = np.abs(out.z-z).argmin()
    # z = out.z[zi]

    if fig_name is None:
        if out.fileName() == '':
            fig = plt.figure('Bunch profile at ' + str(z) + 'm')
        else:
            fig = plt.figure('Bunch profile at ' + str(z) + 'm ' + out.fileName())
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

    if not out.tdp:
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
            subfig_z_power_curr(ax[-1], out, zi=zi, x_units=x_units, legend=legend)
        elif param == 'el_energy+el_espread+el_bunching':
            subfig_z_energy_espread_bunching(ax[-1], out, zi=zi, x_units=x_units, legend=legend)
        elif param == 'el_energy+el_espread':
            subfig_z_energy_espread(ax[-1], out, zi=zi, x_units=x_units, legend=legend)
        elif param == 'rad_phase':
            subfig_z_phase(ax[-1], out, zi=zi, x_units=x_units, legend=legend)
        # elif param == 'el_energy':
        # subfig_evo_el_energy(ax[-1], out, legend)
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
            subfig_z_spec(ax[-1], out, zi=zi, y_units=y_units, estimate_ph_sp_dens=True, legend=legend)
        else:
            _logger.warning(ind_str + 'wrong parameter ' + param)
    axf = len(ax) - axt
    # ax[0].set_xlim(out.z[0], out.z[-1])
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
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        fig.savefig(out.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')


@if_plottable
def subfig_z_power_curr(ax_curr, out, zi=None, x_units='um', legend=False):
    ax_curr.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_curr.set_xlabel(r's [$\mu$m]')
        x = out.t * speed_of_light * 1e6
    elif x_units == 'fs':
        ax_curr.set_xlabel(r't [fs]')
        x = out.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    ax_curr.plot(x, out.I / 1e3, 'k--')
    ax_curr.set_ylabel(r'I [kA]')
    ax_curr.set_ylim(ymin=0)
    ax_curr.text(0.02, 0.98, "Q= %.2f pC" % (out.beam_charge * 1e12), fontsize=12, horizontalalignment='left',
                 verticalalignment='top', transform=ax_curr.transAxes,
                 color='black')  # horizontalalignment='center', verticalalignment='center',
    ax_curr.grid(True)

    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(x, out.rad_power[zi, :], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(out.rad_power[:,zi])!=np.amin(out.rad_power[:,zi]):
    # ax_power.set_ylim([0, np.amax(out.rad_power[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    if 'n_photons' in dir(out):
        ax_curr.text(0.98, 0.98, "E= %.2e J\nN$_{phot}$= %.2e" % (out.rad_energy[zi], out.n_photons[zi]), fontsize=12,
                     horizontalalignment='right', verticalalignment='top', transform=ax_curr.transAxes,
                     color='green')  # horizontalalignment='center', verticalalignment='center',
    else:
        ax_curr.text(0.98, 0.98, "E= %.2e J" % (out.rad_energy[zi]), fontsize=12, horizontalalignment='right',
                     verticalalignment='top', transform=ax_curr.transAxes,
                     color='green')  # horizontalalignment='center', verticalalignment='center',

    ax_curr.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())

    ax_power.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_energy_espread_bunching(ax_energy, out, zi=None, x_units='um', legend=False):
    ax_energy.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_energy.set_xlabel(r's [$\mu$m]')
        x = out.t * speed_of_light * 1e6
    elif x_units == 'fs':
        ax_energy.set_xlabel(r't [fs]')
        x = out.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    ax_energy.plot(x, out.h5['Beam/energy'][zi, :] * m_e_GeV, 'b-', x,
                   (out.h5['Beam/energy'][zi, :] + out.h5['Beam/energyspread'][zi, :]) * m_e_GeV, 'r--', x,
                   (out.h5['Beam/energy'][zi, :] - out.h5['Beam/energyspread'][zi, :]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(True)
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(x, out.h5['Beam/bunching'][zi, :], 'grey', linewidth=0.5)
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
def subfig_z_energy_espread(ax_energy, out, zi=None, x_units='um', legend=False):
    ax_energy.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_energy.set_xlabel(r's [$\mu$m]')
        x = out.t * speed_of_light * 1.0e-15 * 1e6
    elif x_units == 'fs':
        ax_energy.set_xlabel(r't [fs]')
        x = out.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    ax_energy.plot(x, out.h5['Beam/energy'][zi, :] * m_e_GeV, 'b-', x,
                   (out.h5['Beam/energy'][zi, :] + out.h5['Beam/energyspread'][zi, :]) * m_e_GeV, 'r--', x,
                   (out.h5['Beam/energy'][zi, :] - out.h5['Beam/energyspread'][zi, :]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(True)
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')

    ax_energy.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_phase(ax_phase, out, zi=None, x_units='um', legend=False, rewrap=False):
    ax_phase.clear()
    number_ticks = 6

    if x_units == 'um':
        ax_phase.set_xlabel(r's [$\mu$m]')
        x = out.t * speed_of_light * 1e6
    elif x_units == 'fs':
        ax_phase.set_xlabel(r't [fs]')
        x = out.t
    else:
        raise ValueError('Unknown parameter x_units (should be um or fs)')

    if zi == None:
        zi = -1

    #    if rewrap:
    #        phase = unwrap(out.h5['Beam/phase-nearfield'][zi, :])
    #        phase_cor = np.arange(out.nSlices) * (maxspectrum_wavelength - out.lambdaref) / out.lambdaref * out('zsep') * 2 * pi
    #        phase_fixed = phase + phase_cor
    #        phase_fixed -= power[maxspower_index, zi]
    #        n = 1
    #        phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
    #    else:
    phase_fixed = out.h5['Field/phase-nearfield'][zi, :]
    ax_phase.plot(x, phase_fixed, 'k-', linewidth=0.5)
    ax_phase.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top',
                  transform=ax_phase.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_phase.set_ylabel(r'$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])
    ax_phase.grid(True)

    ax_phase.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_phase.set_xlim([x[0], x[-1]])


@if_plottable
def subfig_z_spec(ax_spectrum, out, zi=None, loc='near', y_units='ev', estimate_ph_sp_dens=True, legend=False):
    number_ticks = 6

    if zi == None:
        zi = -1

    scale_ev, spec = out.calc_spec(zi=zi, loc=loc, estimate_ph_sp_dens=True)
    #    if 'spec' not in dir(out):
    #        out.calc_spec()

    if y_units == 'nm':
        x = h_eV_s * speed_of_light * 1e9 / scale_ev
    elif y_units in ['ev', 'eV']:
        x = scale_ev

    if estimate_ph_sp_dens:
        y_units = 'ev'
        # spec = out.spec_phot_density[:, zi]
        # # spec = calc_ph_sp_dens(out.spec[:, zi], out.freq_ev, out.n_photons[zi])
    # else:
    # spec = out.spec[:, zi]

    # power = np.pad(out.p_mid, [(int(out.nSlices / 2) * n_pad, (out.nSlices - (int(out.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')
    # phase = np.pad(out.phi_mid, [(int(out.nSlices / 2) * n_pad, (out.nSlices - (int(out.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')  # not supported by the numpy 1.6.2

    ax_spectrum.plot(x, spec, 'r-')
    ax_spectrum.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top',
                     transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',

    ax_spectrum.set_ylim(ymin=0)
    ax_spectrum.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_spectrum.get_yaxis().get_major_formatter().set_scientific(True)
    ax_spectrum.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    ax_spectrum.grid(True)
    if np.amin(x) != np.amax(x):
        ax_spectrum.set_xlim([np.amin(x), np.amax(x)])

    maxspectrum_index = np.argmax(spec)
    # maxspower_index = np.argmax(power[:, zi])
    maxspectrum_value = x[maxspectrum_index]

    spec_width = None
    if np.sum(spec) != 0:
        pos, width, arr = fwhm3(spec)
        if width != None:
            if arr[0] == arr[-1]:
                dx = abs(x[pos] - x[pos - 1])
            else:
                dx = abs((x[arr[0]] - x[arr[-1]]) / (arr[0] - arr[-1]))
            spec_width = dx * width / x[
                pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

    if spec_width is not None and maxspectrum_value is not None:
        if y_units == 'nm':
            ax_spectrum.text(0.02, 0.98, r"$\lambda^{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (
            maxspectrum_value * 1e-9, spec_width), fontsize=12, horizontalalignment='left', verticalalignment='top',
                             transform=ax_spectrum.transAxes,
                             color='red')  # horizontalalignment='center', verticalalignment='center',
            if estimate_ph_sp_dens:
                ax_spectrum.set_ylabel(r'[$N_{phot}$/nm](estim)')
            else:
                ax_spectrum.set_ylabel(r'P($\lambda$) [a.u.]')
            ax_spectrum.set_xlabel(r'$\lambda$ [nm]')
        elif y_units in ['ev', 'eV']:
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
def plot_gen4_out_e(out, legend=False, figsize=3, fig_name='Electrons', savefig=False, showfig=True, debug=1):
    fig = plot_gen4_out_evo(out, params=['und_quad', 'el_size', 'el_energy', 'el_bunching'], figsize=figsize,
                            legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)


@if_plottable
def plot_gen4_out_ph(out, legend=False, figsize=3, fig_name='Radiation', savefig=False, showfig=True, debug=1):
    if out.tdp:
        fig = plot_gen4_out_evo(out, params=['rad_pow_en_log', 'rad_pow_en_lin', 'rad_spec_log', 'rad_size'],
                                figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig,
                                debug=debug)
    else:
        fig = plot_gen4_out_evo(out, params=['rad_pow_log', 'rad_size'], figsize=figsize, legend=legend,
                                fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)


@if_plottable
def plot_gen4_out_evo(out, params=['und_quad', 'el_size', 'el_pos', 'el_energy', 'el_bunching', 'rad_pow_en_log',
                                   'rad_pow_en_lin', 'rad_spec_log', 'rad_size', 'rad_spec_evo_n', 'rad_pow_evo_n'],
                      figsize=3, legend=False, fig_name=None, savefig=False, showfig=True, debug=1):
    """
    plots evolution of given parameters from genesis output with undulator length
    """
    import matplotlib.ticker as ticker

    if showfig == False and savefig == False:
        return

    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')

    if os.path.isfile(str(out)):
        out = read_out_file(out, read_level=2)
    # add check for output object
    if fig_name is None:
        if out.fileName() == '':
            fig = plt.figure(params_str)
            _logger.info('plotting ' + params_str)
        else:
            fig = plt.figure(out.fileName() + '_' + params_str)
            _logger.info('plotting ' + out.fileName() + '_' + params_str)
    else:
        fig = plt.figure(fig_name)
        if debug > 0:
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
    is_tdp = out.tdp
    for index, param in enumerate(params):
        if len(ax) == 0:
            ax.append(fig.add_subplot(len(params), 1, index + 1))
        else:
            ax.append(fig.add_subplot(len(params), 1, index + 1, sharex=ax[0]))
        # ax[-1]
        if param == 'und_quad':
            subfig_evo_und_quad(ax[-1], out, legend)
        elif param == 'und':
            subfig_evo_und(ax[-1], out, legend)
        elif param == 'el_size':
            subfig_evo_el_size(ax[-1], out, legend)
        elif param == 'el_pos':
            subfig_evo_el_pos(ax[-1], out, legend)
        elif param == 'el_energy':
            subfig_evo_el_energy(ax[-1], out, legend)
        elif param == 'el_bunching':
            subfig_evo_el_bunching(ax[-1], out, legend)
        elif param == 'rad_pow_en_log':
            if not out.tdp:
                subfig_evo_rad_pow(ax[-1], out, legend)
            else:
                subfig_evo_rad_pow_en(ax[-1], out, legend)
        elif param == 'rad_pow_en_lin':
            if not out.tdp:
                subfig_evo_rad_pow(ax[-1], out, legend, log=0)
            else:
                subfig_evo_rad_pow_en(ax[-1], out, legend, log=0)
        elif param == 'rad_pow_log':
            subfig_evo_rad_pow(ax[-1], out, legend)
        elif param == 'rad_pow_lin':
            subfig_evo_rad_pow(ax[-1], out, legend, log=0)
        elif param == 'rad_size':
            subfig_rad_size(ax[-1], out, legend)
        elif param == 'rad_spec_log':
            if out.tdp:
                subfig_evo_rad_spec(ax[-1], out, legend)
        elif param == 'rad_spec_nolog':
            if out.tdp:
                subfig_evo_rad_spec(ax[-1], out, legend, log=0)
        elif param == 'rad_spec_evo_n':
            if out.tdp:
                subfig_evo_rad_spec_sz(ax[-1], out, legend, norm=1)
        elif param == 'rad_pow_evo_n':
            if out.tdp:
                subfig_evo_rad_pow_sz(ax[-1], out, legend, norm=1)
        elif param == 'rad_spec_evo':
            if out.tdp:
                subfig_evo_rad_spec_sz(ax[-1], out, legend, norm=0)
        elif param == 'rad_pow_evo':
            if out.tdp:
                subfig_evo_rad_pow_sz(ax[-1], out, legend, norm=0)
        else:
            _logger.warning('wrong parameter ' + param)

    ax[0].set_xlim(out.z[0], out.z[-1])
    ax[-1].set_xlabel('z [m]')
    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.8, left=0.15)

    for axi in ax[0:-1]:
        for label in axi.get_xticklabels():
            label.set_visible(False)

    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if fig_name == 'Electrons':
            savepath = out.filePath + '_elec.' + str(savefig)
        elif fig_name == 'Radiation':
            savepath = out.filePath + '_rad.' + str(savefig)
        elif fig_name == '':
            savepath = out.filePath + '_' + params_str + '.' + str(savefig)
        else:
            savepath = out.filePath + '_' + fig_name + '.' + str(savefig)
        _logger.debug('saving figure to {}'.format(savepath))
        fig.savefig(out.filePath + '_' + params_str + '.' + str(savefig), format=savefig)
    plt.draw()
    if showfig == True:
        dir_lst = out.filePath.split(os.path.sep)
        dir = os.path.sep.join(dir_lst[0:-1]) + os.path.sep
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        plt.close('all')


@if_plottable
def subfig_evo_und_quad(ax_und, out, legend):
    number_ticks = 6
    aw = out.h5['Lattice/aw']
    qf = out.h5['Lattice/qf']
    z = out.h5['Lattice/z']

    ax_und.step(z, aw, 'b-', where='post', linewidth=1.5)
    # ax_und.scatter(z, aw)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(True)

    ax_quad = ax_und.twinx()
    ax_quad.step(z, qf, 'r-', where='post', linewidth=1.5)
    # ax_quad.scatter(z, qf)
    ax_quad.set_ylabel('Quad')
    ax_quad.grid(False)

    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_quad.yaxis.major.locator.set_params(nbins=number_ticks)

    if np.amax(aw) != 0:
        aw_tmp = np.array(aw)[np.array(aw) != 0]
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
def subfig_evo_und(ax_und, out, legend):
    number_ticks = 6
    aw = out.h5['Lattice/aw']
    qf = out.h5['Lattice/qf']
    z = out.h5['Lattice/z']

    ax_und.step(z, aw, 'b-', where='post', linewidth=1.5)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(True)

    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)

    if np.amax(aw) != 0:
        aw_tmp = np.array(aw)[np.array(aw) != 0]
        if np.amax(aw_tmp) != np.amin(aw_tmp):
            diff = np.amax(aw_tmp) - np.amin(aw_tmp)
            ax_und.set_ylim([np.amin(aw_tmp) - diff / 10, np.amax(aw_tmp) + diff / 10])
    else:
        ax_und.set_ylim([0, 1])
    ax_und.tick_params(axis='y', which='both', colors='b')
    ax_und.yaxis.label.set_color('b')


@if_plottable
def subfig_evo_el_size(ax_size_tsize, out, legend, which='both'):
    number_ticks = 6

    xrms = out.h5['Beam/xsize']
    yrms = out.h5['Beam/ysize']

    # x = out.h5['Beam/xsize']
    # y = out.h5['Beam/ysize']
    z = out.h5['Lattice/zplot']

    if np.sum(out.I) == 0:
        weights = None
        _logger.warning('charge=0, no weighting')
    else:
        weights = out.I

    if which == 'both' or which == 'averaged':
        ax_size_tsize.plot(z, np.average(xrms, axis=1, weights=weights) * 1e6, 'g-', z,
                           np.average(yrms, axis=1, weights=weights) * 1e6, 'b-')
    if which == 'both' or which == 'peak_curr':
        if weights is None:
            idx_pk = int(out.I.size / 2)
        else:
            idx_pk = np.where(out.I == np.amax(out.I))[0][0]
        ax_size_tsize.plot(z, xrms[:, idx_pk] * 1e6, 'g--', z, yrms[:, idx_pk] * 1e6, 'b--')
    ax_size_tsize.set_ylabel(r'$\sigma_{x,y}$ [$\mu$m]')

    ax_size_tsize.set_ylim(ymin=0)
    ax_size_tsize.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_size_tsize.grid(True)


@if_plottable
def subfig_evo_el_pos(ax_size_tpos, out, legend, which='both'):
    number_ticks = 6

    x = out.h5['Beam/xposition']
    y = out.h5['Beam/yposition']
    z = out.h5['Lattice/zplot']

    # if hasattr(out,'x') and hasattr(out,'y'):
    if which == 'both' or which == 'averaged':
        ax_size_tpos.plot(z, np.average(x, axis=1, weights=out.I) * 1e6, 'g-', z,
                          np.average(y, axis=1, weights=out.I) * 1e6, 'b-')
    if which == 'both' or which == 'peak_curr':
        idx_pk = np.where(out.I == np.amax(out.I))[0][0]
        ax_size_tpos.plot(z, x[:, idx_pk] * 1e6, 'g--', z, y[:, idx_pk] * 1e6, 'b--')
    ax_size_tpos.set_ylabel(r'$x,y$ [$\mu$m]')


@if_plottable
def subfig_evo_el_energy(ax_energy, out, legend):
    number_ticks = 6

    el_energy = out.h5['Beam/energy'][:] * m_e_MeV
    el_energy_av = int(np.nanmean(el_energy))
    z = out.h5['Lattice/zplot']
    el_energy_spread = out.h5['Beam/energyspread'][:]

    ax_energy.plot(z, np.average(el_energy - el_energy_av, axis=1), 'b-', linewidth=1.5)
    ax_energy.set_ylabel('<E> + ' + str(el_energy_av) + '[MeV]')
    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.grid(True)

    ax_spread = ax_energy.twinx()
    ax_spread.plot(z, np.average(el_energy_spread * m_e_MeV, weights=out.I, axis=1), 'm--', out.z,
                   np.amax(el_energy_spread * m_e_GeV * 1000, axis=1), 'r--', linewidth=1.5)
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
def subfig_evo_el_bunching(ax_bunching, out, legend):
    number_ticks = 6

    z = out.h5['Lattice/zplot']
    b = out.h5['Beam/bunching']

    ax_bunching.plot(z, np.average(b, weights=out.I, axis=1), 'k-', out.z, np.amax(b, axis=1), 'grey', linewidth=1.5)
    # ax_bunching.plot(out.z, np.amax(out.bunching, axis=0), 'grey',linewidth=1.5) #only max
    ax_bunching.set_ylabel(r'Bunching')
    ax_bunching.set_ylim(ymin=0)
    # ax_bunching.set_ylim([0,0.8])
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.grid(True)


@if_plottable
def subfig_evo_rad_pow_en(ax_rad_pow, out, legend, log=1):
    if log:
        e = out.rad_energy
        e[e == 0] = e[e != 0].min() / 10
        growth = np.divide(np.roll(e, -1), e)
        idx = growth < 2
    else:
        idx = np.ones_like(out.z).astype(bool)

    ax_rad_pow.plot(out.z[idx], np.amax(out.rad_power[idx, :], axis=1), 'g-', linewidth=1.5)
    ax_rad_pow.set_ylabel(r'P [W]')
    ax_rad_pow.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_pow.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(out.rad_power) > 0 and log:
        ax_rad_pow.set_yscale('log')
    if not log:
        ax_rad_pow.set_ylim(ymin=0)
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_en = ax_rad_pow.twinx()
    ax_rad_en.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_en.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(out.rad_power) > 0 and log:
        ax_rad_en.plot(out.z[idx], out.rad_energy[idx], 'k--', linewidth=1.5)
        ax_rad_en.set_ylabel(r'E [J]')
        ax_rad_en.set_yscale('log')
    if not log:
        if np.amax(out.rad_energy) < 1e-4:
            ax_rad_en.plot(out.z[idx], out.rad_energy[idx] * 1e6, 'k--', linewidth=1.5)
            ax_rad_en.set_ylabel(r'E [$\mu$J]')
        else:
            ax_rad_en.plot(out.z[idx], out.rad_energy[idx] * 1e3, 'k--', linewidth=1.5)
            ax_rad_en.set_ylabel(r'E [mJ]')
        ax_rad_en.set_ylim(ymin=0)
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_pow.grid(True)  # , which='minor')
    # ax_rad_pow.grid(False, which="minor")
    ax_rad_pow.tick_params(axis='y', which='both', colors='g')
    ax_rad_pow.yaxis.label.set_color('g')
    ax_rad_en.tick_params(axis='y', which='both', colors='k')
    ax_rad_en.yaxis.label.set_color('k')
    ax_rad_en.grid(False)
    # ax_rad_en.grid(False, which='minor')
    ax_rad_pow.yaxis.get_offset_text().set_color(ax_rad_pow.yaxis.label.get_color())
    ax_rad_en.yaxis.get_offset_text().set_color(ax_rad_en.yaxis.label.get_color())

    ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W ' '\n' r'$E_{end}$= %.2e J' % (
    np.amax(out.rad_power[-1, :]), out.rad_energy[-1]), fontsize=12, horizontalalignment='right',
                    verticalalignment='bottom', transform=ax_rad_pow.transAxes)


@if_plottable
def subfig_evo_rad_pow(ax_rad_pow, out, legend, log=1):
    ax_rad_pow.plot(out.z, np.amax(out.rad_power, axis=1), 'g-', linewidth=1.5)
    ax_rad_pow.set_ylabel('P [W]')
    ax_rad_pow.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_rad_pow.get_yaxis().get_major_formatter().set_scientific(True)
    if np.amax(out.rad_power) > 0 and log:
        ax_rad_pow.set_yscale('log')
    plt.yticks(plt.yticks()[0][0:-1])

    ax_rad_pow.grid(False)  # , which='minor')
    ax_rad_pow.tick_params(axis='y', which='both', colors='g')
    ax_rad_pow.yaxis.label.set_color('g')
    ax_rad_pow.yaxis.get_offset_text().set_color(ax_rad_pow.yaxis.label.get_color())
    ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W' % (np.amax(out.rad_power[-1, :])), fontsize=12,
                    horizontalalignment='right', verticalalignment='bottom', transform=ax_rad_pow.transAxes)


@if_plottable
def subfig_evo_rad_spec(ax_spectrum, out, legend, log=1):
    scale_ev, spec = out.calc_spec()

    ax_spectrum.plot(out.z, np.amax(spec, axis=1), 'r-', linewidth=1.5)
    ax_spectrum.text(0.5, 0.98, r"(on axis)", fontsize=10, horizontalalignment='center', verticalalignment='top',
                     transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel(r'P$(\lambda)_{max}$ [a.u.]')
    plt.yticks(plt.yticks()[0][0:-1])

    if np.amax(np.amax(spec, axis=0)) > 0 and log:
        ax_spectrum.set_yscale('log')
    ax_spectrum.grid(True)

    spectrum_lamdwidth_fwhm = np.zeros_like(out.z)
    spectrum_lamdwidth_std = np.zeros_like(out.z)

    for zz in range(out.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(spec[zz, :]) != 0:
            pos, width, arr = fwhm3(spec[zz, :])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(scale_ev[pos] - scale_ev[pos - 1])
                else:
                    dlambda = abs((scale_ev[arr[0]] - scale_ev[arr[-1]]) / (arr[0] - arr[-1]))
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / scale_ev[pos]
                # spectrum_lamdwidth_fwhm[zz] = abs(out.freq_lamd[arr[0]] - out.freq_lamd[arr[-1]]) / out.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

            spectrum_lamdwidth_std[zz] = std_moment(scale_ev, spec[zz, :]) / n_moment(scale_ev, spec[zz, :], 0, 1)

        # try:
        #     peak = fwhm3(out.spec[:, zz])
        #     spectrum_lamdwidth_fwhm[zz] = abs(out.freq_lamd[0] - out.freq_lamd[1]) * peak[1] / out.freq_lamd[peak[0]]  # the FWHM of spectral line (error when paekpos is at the edge of lamdscale)
        # except:
        #     spectrum_lamdwidth_fwhm[zz] = 0
        #
        # try:
        #     spectrum_lamdwidth_std[zz] = std_moment(out.freq_lamd, out.spec[:, zz]) / n_moment(out.freq_lamd, out.spec[:, zz], 0, 1)
        # except:
        #     spectrum_lamdwidth_std[zz] = 0

    ax_spec_bandw = ax_spectrum.twinx()
    ax_spec_bandw.plot(out.z, spectrum_lamdwidth_fwhm * 100, 'm--', label="fwhm")
    ax_spec_bandw.plot(out.z, 2 * spectrum_lamdwidth_std * 100, 'm:', label="std")
    ax_spec_bandw.grid(False)
    plt.yticks(plt.yticks()[0][0:-1])
    ax_spec_bandw.set_ylim(ymin=0)

    if legend:
        ax_spec_bandw.legend()
    ax_spec_bandw.set_ylabel(r'$\Delta\lambda/\lambda, \%$' + '\n' + r'(-- fwhm, $\cdots2\sigma$)')
    if spectrum_lamdwidth_fwhm[-1] != None:
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
def subfig_rad_size(ax_size_t, out, legend):
    x_size = out.h5['Field/xsize'][:]
    y_size = out.h5['Field/ysize'][:]
    r_size = np.sqrt(x_size ** 2 + y_size ** 2)

    if out.nSlices == 1:
        ax_size_t.plot(out.z, r_size * 2 * 1e6, 'b-', linewidth=1.5)
        #        ax_size_t.plot([np.amin(out.z), np.amax(out.z)], [out.leng * 1e6, out.leng * 1e6], 'b-', linewidth=1.0)
        ax_size_t.set_ylabel('transverse $[\mu m]$')
    else:
        if hasattr(out, r'rad_t_size_weighted'):
            ax_size_t.plot(out.z, out.rad_t_size_weighted, 'b-', linewidth=1.5)
        else:
            if np.amax(out.rad_power) > 0:
                # idx = out.rad_energy != 0
                weight = out.rad_power + np.amin(out.rad_power[out.rad_power != 0]) / 1e6
                weight[out.rad_energy == 0, :] = 1
                r_size = np.average(r_size * 2 * 1e6, weights=weight, axis=1)
            else:
                r_size = np.zeros_like(out.rad_power)
                # weight = np.ones_like(out.rad_power)

            ax_size_t.plot(out.z, r_size, 'b-', linewidth=1.5)

    ax_size_t.set_ylim(ymin=0)
    # ax_size_t.set_ylabel(r'$\sim$size$_{transv}$ [$\mu$m]'+'\n'+r'($2\sigma$)')
    ax_size_t.set_ylabel(r'$\sim$size$_{transv}$ [$\mu$m]' + '\n' + u'(\u2014 $2\sigma$)')
    ax_size_t.grid(True)
    plt.yticks(plt.yticks()[0][0:-1])

    if out.nSlices > 1:
        ax_size_s = ax_size_t.twinx()
        size_long_fwhm = np.zeros_like(out.z)
        size_long_std = np.zeros_like(out.z)
        s = out.t * speed_of_light * 1.0e-15 * 1e6
        delta_s = (s[1] - s[0])
        for zz in range(out.nZ):
            # size_long_fwhm[zz] = fwhm(out.s,out.rad_power[:, zz])
            if np.sum(out.rad_power[zz, :]) != 0:
                # try:
                # _, width, _ = fwhm3(out.rad_power[zz, :])
                # if width != None:
                #     size_long_fwhm[zz] = abs(delta_s) * width
                # else:
                #     size_long_fwhm[zz] = None
                # except:
                #     size_long_fwhm[zz] = 0

                # try:
                size_long_std[zz] = std_moment(s, out.rad_power[zz, :])
                # except:
                #     size_long_std[zz] = 0
            else:
                # size_long_fwhm[zz] = None
                size_long_std[zz] = None

        ax_size_s.plot(out.z, size_long_fwhm, color='navy', linestyle='--', linewidth=1.0, label="fwhm")
        ax_size_s.plot(out.z, 2 * size_long_std, color='navy', linestyle=':', linewidth=1.0, label="std")
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


#        plt.legend('fwhm','std')

@if_plottable
def subfig_evo_rad_pow_sz(ax_power_evo, out, legend, norm=1, **kwargs):
    if out.nSlices > 1:
        z = out.z
        s = out.s
        power = out.rad_power
        if norm == 1:
            max_power = np.nanmax(power, 1)[:, np.newaxis]
            max_power[max_power == 0] = 1  # avoid division by zero
            power = power / max_power
            # power[isnan(power)]=0
        ax_power_evo.pcolormesh(z, s * 1e6, power.T)
        ax_power_evo.set_xlabel('z [m]')
        ax_power_evo.set_ylabel('s [$\mu$m]')
        ax_power_evo.axis('tight')
        ax_power_evo.grid(True)
    else:
        pass


@if_plottable
def subfig_evo_rad_spec_sz(ax_spectrum_evo, out, legend, norm=1):
    if out.nSlices > 1:
        z = out.z
        l, spectrum = out.calc_spec()
        #        spectrum = out.spec
        if norm == 1:
            max_spectrum = np.nanmax(spectrum, 1)[:, np.newaxis]
            max_spectrum[max_spectrum == 0] = 1  # avoid division by zero
            spectrum = spectrum / max_spectrum
            # spectrum[isnan(spectrum)]=0
        ax_spectrum_evo.pcolormesh(z, l, spectrum.T)
        ax_spectrum_evo.set_xlabel('z [m]')
        ax_spectrum_evo.set_ylabel('[eV]')
        ax_spectrum_evo.axis('tight')
        ax_spectrum_evo.grid(True)
    else:
        pass
