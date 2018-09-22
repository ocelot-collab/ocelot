'''
user interface for viewing genesis simulation results
'''

'''
MEMO

plt.gcf() to get current figure
plt.gca() to get current axis

ax.set_xlabel('')
ax.get_xlim()
ax.set_xlim([0, 1])
ax.set_ylim(ymin=0)

'''

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

from ocelot.adaptors.genesis import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.common.math_op import *  # import of mathematical functions
from ocelot.utils.xfel_utils import *
from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.optics.wave import *

from ocelot.gui.colormaps2d.colormap2d import *

# from pylab import rc, rcParams #tmp
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

_logger = logging.getLogger('ocelot.genesis_plot')

def_cmap = 'viridis'
# def_cmap = 'Greys'

fntsz = 4
params = {'image.cmap': def_cmap, 'backend': 'ps', 'axes.labelsize': 3 * fntsz, 'font.size': 3 * fntsz, 'legend.fontsize': 4 * fntsz, 'xtick.labelsize': 4 * fntsz,  'ytick.labelsize': 4 * fntsz, 'text.usetex': False}
rcParams.update(params)
# plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
# rcParams["savefig.directory"] = os.chdir(os.path.dirname(__file__)) but __file__ appears to be genesis_plot

plt.ioff() #turn off interactive mode

def plot_gen_out_all_paral(exp_dir, stage=1, savefig='png', debug=1):
    print('start')
    from ocelot.utils.xfel_utils import background
    i = 0
    dir = exp_dir + 'run_' + str(i) + '/'

    while(os.path.exists(dir)):
        print(i)
        file = dir + 'run.' + str(i) + '.s'+str(stage)+'.gout'
        if(file): 
            print('good',i)
            background('''plot_gen_out_all("'''+file+'''", choice=(1,1,1,1,0,0,0,0,0,0,0),debug='''+str(debug)+''')''')
        
        i += 1
        dir = exp_dir + 'run_' + str(i) + '/'
        print(dir)
    
    return

# plot_gen_stat(proj_dir=exp_dir, run_inp=[], stage_inp=[], param_inp=[], s_param_inp=['p_int','energy','r_size_weighted'], z_param_inp=[], dfl_param_inp=[], s_inp=['max'], z_inp=[0,'end'], savefig=1, saveval=1, showfig=0, debug=0)

def plot_gen_out_all(handle=None, savefig='png', showfig=False, choice='all', vartype_dfl=complex128, debug=1):
    '''
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
    '''


    _logger.info('plotting all genesis output')
    _logger.debug('choice = ' + str(choice))
    plotting_time = time.time()

    # plt.ioff()

    if savefig == True:
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
        _logger.info('\n  plotting all files in ' + str(handle))
    else:
        handles = [handle]

    for handle in handles:

        if os.path.isfile(str(handle)):
            _logger.info('plotting '+str(handle))
            handle = read_out_file(handle, read_level=2, debug=debug)

        if isinstance(handle, GenesisOutput):
            if choice[0]:
                f0 = plot_gen_out_e(handle, showfig=showfig, savefig=savefig, debug=debug)
            if choice[1]:
                f1 = plot_gen_out_ph(handle, showfig=showfig, savefig=savefig, debug=debug)
            if choice[2]:
                f2 = plot_gen_out_z(handle, z=0, showfig=showfig, savefig=savefig, debug=debug)
            if choice[3]:
                f3 = plot_gen_out_z(handle, z=inf, showfig=showfig, savefig=savefig, debug=debug)
            if choice[11] != 0:
                if choice[11] == -1:
                    try:
                        W=wigner_out(handle, pad=2)
                        plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    except:
                        _logger.warning('could not plot wigner')
                else:
                    if choice[11] == 1:
                        _logger.warning('choice[11] in plot_gen_out_all defines interval of Wigner plotting. To plot at the end set to "-1"')
                    try:
                        for z in np.arange(0, np.amax(handle.z), choice[11]):
                            W=wigner_out(handle, z=z, pad=2)
                            plot_wigner(W, showfig=showfig, savefig=savefig, debug=debug, downsample=2)
                    except:
                        _logger.warning('could not plot wigner')
            if choice[4] != 0 and  choice[4] != []:
                for z in np.arange(choice[4], np.amax(handle.z), choice[4]):
                    plot_gen_out_z(handle, z=z, showfig=showfig, savefig=savefig, debug=debug)
            if choice[12]:
                try:
                    plot_dpa_bucket_out(handle,scatter=0,slice_pos='max_P',repeat=3, showfig=showfig, savefig=savefig, cmap=def_cmap)
                except:
                    _logger.warning('could not plot particle buckets')

            
                
        if os.path.isfile(handle.filePath + '.dfl') and any(choice[5:8]):
            dfl = read_dfl_file_out(handle, debug=debug)
            if dfl.Nz()==0:
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
            if np.size(dpa.ph)==0:
                _logger.warning('empty dpa, skipping')
            else:
                if choice[9]:
                    edist = dpa2edist(handle, dpa, num_part=5e4, smear=1, debug=debug)
                    f9 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig, bins=100, debug=debug)
                if choice[10]:
                    edist = dpa2edist(handle, dpa, num_part=5e4, smear=0, debug=debug)
                    f10 = plot_edist(edist, figsize=3, fig_name=None, savefig=savefig, showfig=showfig, bins=(100, 100, 300, 200), debug=debug)
        
    if savefig != False:
        _logger.info(ind_str + 'plots recorded to *.' + str(savefig) + ' files')

    if showfig:
        _logger.info(ind_str + 'showing plots, close all to proceed')
        plt.show()
    # else:
        # plt.close('all')

    _logger.info(ind_str + 'total plotting time {:.2f} seconds'.format(time.time() - plotting_time))


def plot_gen_out_z(g, z=inf, params=['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'rad_spec'], figsize=3.5, x_units='um', y_units='ev', legend=False, fig_name=None, savefig=False, showfig=True, debug=1):
    '''
    radiation parameters at distance z
    g/out = GenesisOutput() object
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
    '''
    import matplotlib.ticker as ticker
    
    
    if showfig == False and savefig == False:
        return

    t_domain = ['rad_power+el_current', 'el_energy+el_espread+el_bunching', 'el_energy+el_espread', 'rad_phase']
    f_domain = ['rad_spec']
    # t_domain_i = list(set(t_domain).intersection(params))
    # f_domain_i = list(set(f_domain).intersection(params))
    # t_domain_n = len(t_domain_i)
    # f_domain_n = len(f_domain_i)
    #add sorting of f_domain to the end params += [params.pop(i)]
    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')
    
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
    
    if os.path.isfile(str(g)):
        g = read_out_file(g, read_level=2)
    # add check for output object

    # if fig_name is None:
        # if g.fileName() is '':
            # fig = plt.figure(params_str)
            # if debug > 0:
                # print('    plotting ' + params_str)
        # else:
            # fig = plt.figure(g.fileName() + '_' + params_str)
            # if debug > 0:
                # print('    plotting ' + g.fileName() + '_' + params_str)
    # else:
        # fig = plt.figure(fig_name)
        # if debug > 0:
            # print('    plotting ' + fig_name)
    if fig_name is None:
        if g.fileName() is '':
            fig = plt.figure('Bunch profile at {:} [m]'.format(z))
        else:
            fig = plt.figure('Bunch profile at {:} [m]'.format(z) + g.fileName())
    else:
        fig = plt.figure(fig_name)
    _logger.info('plotting bunch profile at {:} [m]'.format(z))
        
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
            subfig_z_power_curr(ax[-1], g, zi=zi, x_units=x_units ,legend=legend)
        elif param == 'el_energy+el_espread+el_bunching':
            subfig_z_energy_espread_bunching(ax[-1], g, zi=zi, x_units=x_units ,legend=legend)
        elif param == 'el_energy+el_espread':
            subfig_z_energy_espread(ax[-1], g, zi=zi, x_units=x_units ,legend=legend)
        elif param == 'rad_phase':
            subfig_z_phase(ax[-1], g, zi=zi, x_units=x_units ,legend=legend)
        # elif param == 'el_energy':
            # subfig_evo_el_energy(ax[-1], g, legend)
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
        if len(ax)-axt == 0:
            ax.append(fig.add_subplot(len(params), 1, index + axt + 1))
        else:
            ax.append(fig.add_subplot(len(params), 1, index + axt + 1, sharex=ax[-1]))
        if param == 'rad_spec':
            subfig_z_spec(ax[-1], g, zi=zi, y_units=y_units, estimate_ph_sp_dens=True, legend=legend)
        else:
            print('! wrong parameter ' + param)
    axf = len(ax) - axt
    # ax[0].set_xlim(g.z[0], g.z[-1])
    # ax[-1].set_xlabel('z [m]')
    if axt is not 0 and axf is not 0:
        fig.subplots_adjust(top=0.95, bottom=0.2, right=0.8, left=0.15)
    else:
        fig.subplots_adjust(top=0.95, bottom=0.1, right=0.8, left=0.15)

    for axi in ax[axt:]:
        if axt is not 0:
            pos1 = axi.get_position()  # get the original position
            pos2 = [pos1.x0 + 0, pos1.y0 - 0.1,  pos1.width / 1.0, pos1.height / 1.0]
            axi.set_position(pos2)
            if axi is not ax[-1]:
                axi.set_xlabel('')
                for label in axi.get_xticklabels():
                    label.set_visible(False)

    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        fig.savefig(g.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)

    if showfig:
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


def subfig_z_power_curr(ax_curr, g, zi=None, x_units='um', legend=False):
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
        
    if zi == None:
        zi = -1
    
    ax_curr.plot(x, g.I / 1e3, 'k--')
    ax_curr.set_ylabel(r'I [kA]')
    ax_curr.set_ylim(ymin=0)
    ax_curr.text(0.02, 0.98, "Q= {:.2f} pC".format(g.beam_charge * 1e12), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_curr.transAxes, color='black')  # horizontalalignment='center', verticalalignment='center',
    ax_curr.grid(True)
    
    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(x, g.p_int[:, zi], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(g.p_int[:,zi])!=np.amin(g.p_int[:,zi]):
        # ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    if 'n_photons' in dir(g):
        ax_curr.text(0.98, 0.98, "E= {:.2e} J\nN$_{{phot}}$= {:.2e}".format(g.pulse_energy[zi], g.n_photons[zi]), fontsize=12, horizontalalignment='right', verticalalignment='top', transform=ax_curr.transAxes, color='green')  # horizontalalignment='center', verticalalignment='center',
    else:
        ax_curr.text(0.98, 0.98, "E= {:.2e} J".format(g.pulse_energy[zi]), fontsize=12, horizontalalignment='right', verticalalignment='top', transform=ax_curr.transAxes, color='green')  # horizontalalignment='center', verticalalignment='center',
    
    ax_curr.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
    
    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())
    
    ax_power.set_xlim([x[0],x[-1]])


def subfig_z_energy_espread_bunching(ax_energy, g, zi=None, x_units='um', legend=False):
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
        
    ax_energy.plot(x, g.el_energy[:, zi] * m_e_GeV, 'b-', x, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV, 'r--', x, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(True)
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
    
    ax_energy.set_xlim([x[0],x[-1]])


def subfig_z_energy_espread(ax_energy, g, zi=None, x_units='um', legend=False):
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
        
    ax_energy.plot(x, g.el_energy[:, zi] * m_e_GeV, 'b-', x, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV, 'r--', x, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(True)
    # plt.yticks(plt.yticks()[0][0:-1])
    
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')
    
    ax_energy.set_xlim([x[0],x[-1]])


def subfig_z_phase(ax_phase, g, zi=None, x_units='um', legend=False, **kwargs):
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
    
    if zi == None:
        zi = -1
    
    if "rewrap" in kwargs:
        _logger.warning(ind_str + '"rewrap" argument is obsolete')
    
    if hasattr(g, 'phi_mid_disp'):
        phase_disp = g.phi_mid_disp[:, zi]
    
    else:
        phase_disp = g.phi_mid[:, zi]
        # phase = unwrap(g.phi_mid[:, zi])
        # phase_cor = np.arange(g.nSlices) * (maxspectrum_wavelength - g('xlamds')) / g('xlamds') * g('zsep') * 2 * pi
        # phase_fixed = phase + phase_cor
        # phase_fixed -= power[maxspower_index, zi]
        # n = 1
        # phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
    # else:
        # phase_fixed = g.phi_mid[:, zi]
    ax_phase.plot(x, phase_disp, 'k-', linewidth=0.5)
    if hasattr(g, 'phi_mid_disp'):
        _txt = r'(on axis, rewrapped)'
    else:
        _txt = r'(on axis)'
        ax_phase.text(0.98, 0.98, _txt, fontsize=10, horizontalalignment='right', verticalalignment='top', transform=ax_phase.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_phase.set_ylabel(r'$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])
    ax_phase.grid(True)
    
    ax_phase.yaxis.major.locator.set_params(nbins=number_ticks)
    
    ax_phase.set_xlim([x[0],x[-1]])


def subfig_z_spec(ax_spectrum, g, zi=None, y_units='ev', estimate_ph_sp_dens=True, legend=False, mode='mid'):
    
    number_ticks = 6
    # n_pad = 1
    
    if zi == None:
        zi = -1
    
    if hasattr(g,'spec'):
        if g.spec_mode != mode:
            g.calc_spec(mode = mode)
    else:
        g.calc_spec(mode = mode)
    #TMP
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
    if mode == 'mid':
        ax_spectrum.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top', transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    else:
        ax_spectrum.text(0.98, 0.98, r'(integrated assuming on-axis phases)', fontsize=10, horizontalalignment='right', verticalalignment='top', transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    
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
    if np.sum(spec)!=0:
        pos, width, arr = fwhm3(spec)
        if width != None:
            if arr[0] == arr[-1]:
                dx = abs(x[pos] - x[pos-1])
            else:
                dx = abs( (x[arr[0]] - x[arr[-1]]) / (arr[0] - arr[-1]) )
            spec_width = dx * width / x[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

    if spec_width is not None and maxspectrum_value is not None:
            if y_units == 'nm':
                ax_spectrum.text(0.02, 0.98, r"$\lambda^{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (maxspectrum_value*1e-9, spec_width), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_spectrum.transAxes, color='red')  # horizontalalignment='center', verticalalignment='center',
                if estimate_ph_sp_dens:
                    ax_spectrum.set_ylabel(r'[$N_{phot}$/nm](estim)')
                else:
                    ax_spectrum.set_ylabel(r'P($\lambda$) [a.u.]')
                ax_spectrum.set_xlabel(r'$\lambda$ [nm]')
            elif y_units in ['ev', 'eV']:
                ax_spectrum.text(0.02, 0.98, r"$E_{ph}^{max}$= %.2f eV " "\n" "$(\Delta E/E)_{fwhm}$= %.2e" % (maxspectrum_value, spec_width), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_spectrum.transAxes, color='red')  # horizontalalignment='center', verticalalignment='center',
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
    
    

def plot_gen_out_z_old(g, figsize=(10, 14), x_units='um', y_units='ev', legend=True, fig_name=None, z=inf, savefig=False, showfig=1, debug=1):
    print('soon this function will be replaced by plot_gen_out_z_new (currently being tested)')
    number_ticks = 6

    if showfig == False and savefig == False:
        return

    if g('itdp') == False:
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
        if g.fileName() is '':
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
    ax_curr.text(0.02, 0.98, "Q= %.2f pC" % (g.beam_charge * 1e12), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_curr.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_curr.text(0.98, 0.98, "E= %.2e J" % (g.pulse_energy[zi]), fontsize=12, horizontalalignment='right', verticalalignment='top', transform=ax_curr.transAxes, color='green')  # horizontalalignment='center', verticalalignment='center',
    ax_curr.grid(True)

    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(s, g.p_int[:, zi], 'g-', linewidth=1.5)
    ax_power.set_ylabel(r'Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(g.p_int[:,zi])!=np.amin(g.p_int[:,zi]):
        # ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]

    # ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:, zi] * m_e_GeV, 'b-', s, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV, 'r--', s, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
    ax_energy.set_ylabel(r'$E\pm\sigma_E$ [GeV]')
    # ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')
    ax_energy.grid(True)
    # plt.yticks(plt.yticks()[0][0:-1])

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(s, g.bunching[:, zi], 'grey', linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.set_ylim(ymin=0)
    ax_bunching.grid(False)

    n_pad = 1
    power = np.pad(g.p_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')
    phase = np.pad(g.phi_mid, [(int(g.nSlices / 2) * n_pad, (g.nSlices - (int(g.nSlices / 2)))) * n_pad, (0, 0)], mode='constant')  # not supported by the numpy 1.6.2

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
    ax_spectrum.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top', transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel(r'P($\lambda$) [a.u.]')
    ax_spectrum.set_xlabel(r'$\lambda$ [nm]')
    ax_spectrum.set_ylim(ymin=0)
    ax_spectrum.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_spectrum.get_yaxis().get_major_formatter().set_scientific(True)
    ax_spectrum.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))  # [:,75,75]
    ax_spectrum.grid(True)
    if np.amin(g.freq_lamd) != np.amax(g.freq_lamd):
        ax_spectrum.set_xlim([np.amin(g.freq_lamd), np.amax(g.freq_lamd)])
    ax_phase.set_xlabel(r's [$\mu$m]')

    maxspectrum_index = np.argmax(g.spec[:, zi])
    maxspower_index = np.argmax(power[:, zi])
    maxspectrum_wavelength = g.freq_lamd[maxspectrum_index] * 1e-9
    
    spectrum_lamdwidth_fwhm = None
    if np.sum(g.spec[:,zi])!=0:
        pos, width, arr = fwhm3(g.spec[:, zi])
        if width != None:
            if arr[0] == arr[-1]:
                dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos-1])
            else:
                dlambda = abs( (g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]) )
            spectrum_lamdwidth_fwhm = dlambda * width / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)

    if spectrum_lamdwidth_fwhm is not None and maxspectrum_wavelength is not None:
        ax_spectrum.text(0.02, 0.98, r"$\lambda_{max}$= %.4e m " "\n" "$(\Delta\lambda/\lambda)_{fwhm}$= %.2e" % (maxspectrum_wavelength, spectrum_lamdwidth_fwhm), fontsize=12, horizontalalignment='left', verticalalignment='top', transform=ax_spectrum.transAxes, color='red')  # horizontalalignment='center', verticalalignment='center',

    phase = unwrap(g.phi_mid[:, zi])

    phase_cor = np.arange(g.nSlices) * (maxspectrum_wavelength - g('xlamds')) / g('xlamds') * g('zsep') * 2 * pi
    phase_fixed = phase + phase_cor
    phase_fixed -= power[maxspower_index, zi]
    n = 1
    phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
    ax_phase.plot(s, phase_fixed, 'k-', linewidth=0.5)
    ax_phase.text(0.98, 0.98, r'(on axis)', fontsize=10, horizontalalignment='right', verticalalignment='top', transform=ax_phase.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_phase.set_ylabel(r'$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])
    ax_phase.grid(True)
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
    pos2 = [pos1.x0 + 0, pos1.y0 - 0.1,  pos1.width / 1.0, pos1.height / 0.9]
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
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        fig.savefig(g.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')



def plot_gen_out_e(g, legend=False, figsize=4, fig_name='Electrons', savefig=False, showfig=True, debug=1):
    fig = plot_gen_out_evo(g, params=['und_quad', 'el_size', 'el_energy', 'el_bunching'], figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)


def plot_gen_out_ph(g, legend=False, figsize=4, fig_name='Radiation', savefig=False, showfig=True, debug=1):
    if g('itdp'):
        fig = plot_gen_out_evo(g, params=['rad_pow_en_log', 'rad_pow_en_lin', 'rad_spec_log', 'rad_size'], figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)
    else:
        fig = plot_gen_out_evo(g, params=['rad_pow_log', 'rad_size'], figsize=figsize, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, debug=debug)

def plot_gen_out_evo(g, params=['und_quad', 'el_size', 'el_pos', 'el_energy', 'el_bunching', 'rad_pow_en_log', 'rad_pow_en_lin', 'rad_spec_log', 'rad_size', 'rad_spec_evo_n', 'rad_pow_evo_n'], figsize=4, legend=False, fig_name=None, savefig=False, showfig=True, debug=1):
    '''
    plots evolution of given parameters from genesis output with undulator length
    '''
    import matplotlib.ticker as ticker
    
    if showfig == False and savefig == False:
        return
    
    params_str = str(params).replace("'", '').replace('[', '').replace(']', '').replace(' ', '').replace(',', '--')
    
    if os.path.isfile(str(g)):
        g = read_out_file(g, read_level=2)
    # add check for output object
    if fig_name is None:
        if g.fileName() is '':
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
            subfig_evo_und_quad(ax[-1], g, legend)
        elif param == 'und':
            subfig_evo_und(ax[-1], g, legend)
        elif param == 'el_size':
            subfig_evo_el_size(ax[-1], g, legend)
        elif param == 'el_pos':
            subfig_evo_el_pos(ax[-1], g, legend)
        elif param == 'el_energy':
            subfig_evo_el_energy(ax[-1], g, legend)
        elif param == 'el_bunching':
            subfig_evo_el_bunching(ax[-1], g, legend)
        elif param == 'rad_pow_en_log':
            if not is_tdp:
                subfig_evo_rad_pow(ax[-1], g, legend)
            else:
                subfig_evo_rad_pow_en(ax[-1], g, legend)
        elif param == 'rad_pow_en_lin':
            if not is_tdp:
                subfig_evo_rad_pow(ax[-1], g, legend, log=0)
            else:
                subfig_evo_rad_pow_en(ax[-1], g, legend, log=0)
        elif param == 'rad_pow_log':
            subfig_evo_rad_pow(ax[-1], g, legend)
        elif param == 'rad_pow_lin':
            subfig_evo_rad_pow(ax[-1], g, legend, log=0)
        elif param == 'rad_size':
            subfig_rad_size(ax[-1], g, legend)
        elif param == 'rad_spec_log':
            if is_tdp:
                subfig_evo_rad_spec(ax[-1], g, legend)
        elif param == 'rad_spec_nolog':
            if is_tdp:
                subfig_evo_rad_spec(ax[-1], g, legend, log=0)
        elif param == 'rad_spec_evo_n':
            if is_tdp:
                subfig_evo_rad_spec_sz(ax[-1], g, legend, norm=1)
        elif param == 'rad_pow_evo_n':
            if is_tdp:
                subfig_evo_rad_pow_sz(ax[-1], g, legend, norm=1)
        elif param == 'rad_spec_evo':
            if is_tdp:
                subfig_evo_rad_spec_sz(ax[-1], g, legend, norm=0)
        elif param == 'rad_pow_evo':
            if is_tdp:
                subfig_evo_rad_pow_sz(ax[-1], g, legend, norm=0)
        else:
            print('! wrong parameter ' + param)

    ax[0].set_xlim(g.z[0], g.z[-1])
    ax[-1].set_xlabel('z [m]')
    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.8, left=0.15)

    for axi in ax[0:-1]:
        for label in axi.get_xticklabels():
            label.set_visible(False)
    
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if fig_name == 'Electrons':
            fig.savefig(g.filePath + '_elec.' + str(savefig), format=savefig)
        elif fig_name == 'Radiation':
            fig.savefig(g.filePath + '_rad.' + str(savefig), format=savefig)
        else:
            fig.savefig(g.filePath + '_' + params_str + '.' + str(savefig), format=savefig)

    plt.draw()
    if showfig == True:
        dir_lst = g.filePath.split(os.path.sep)
        dir = os.path.sep.join(dir_lst[0:-1]) + os.path.sep
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


def subfig_evo_und_quad(ax_und, g, legend):
    number_ticks = 6
    
    ax_und.plot(g.z, g.aw, 'b-', linewidth=1.5)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(True)
    
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


def subfig_evo_und(ax_und, g, legend):
    number_ticks = 6
    
    ax_und.plot(g.z, g.aw, 'b-', linewidth=1.5)
    ax_und.set_ylabel('K (rms)')
    ax_und.grid(True)

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


def subfig_evo_el_size(ax_size_tsize, g, legend, which='both'):
    number_ticks = 6
    
    if which == 'both' or which == 'averaged':
        ax_size_tsize.plot(g.z, np.average(g.xrms, axis=0, weights=g.I) * 1e6, 'g-', g.z, np.average(g.yrms, axis=0, weights=g.I) * 1e6, 'b-')
    if which == 'both' or which == 'peak_curr':
        idx_pk = np.where(g.I == np.amax(g.I))[0][0]
        ax_size_tsize.plot(g.z, g.xrms[idx_pk, :] * 1e6, 'g--', g.z, g.yrms[idx_pk, :] * 1e6, 'b--')
    ax_size_tsize.set_ylabel(r'$\sigma_{x,y}$ [$\mu$m]')

    ax_size_tsize.set_ylim(ymin=0)
    ax_size_tsize.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_size_tsize.grid(True)
    
def subfig_evo_el_pos(ax_size_tpos, g, legend, which='both'):
    number_ticks = 6
    
    if hasattr(g,'x') and hasattr(g,'y'):
        if which == 'both' or which == 'averaged':
            ax_size_tpos.plot(g.z, np.average(g.x, axis=0, weights=g.I) * 1e6, 'g-', g.z, np.average(g.y, axis=0, weights=g.I) * 1e6, 'b-')
        if which == 'both' or which == 'peak_curr':
            idx_pk = np.where(g.I == np.amax(g.I))[0][0]
            ax_size_tpos.plot(g.z, g.x[idx_pk, :] * 1e6, 'g--', g.z, g.y[idx_pk, :] * 1e6, 'b--')
        ax_size_tpos.set_ylabel(r'$x,y$ [$\mu$m]')

def subfig_evo_el_energy(ax_energy, g, legend):
    number_ticks = 6
    
    el_energy = g.el_energy * m_e_MeV
    el_energy_av = int(np.mean(el_energy))
    ax_energy.plot(g.z, np.average(el_energy - el_energy_av, axis=0), 'b-', linewidth=1.5)
    ax_energy.set_ylabel('E + ' + str(el_energy_av) + '[MeV]')
    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.grid(True)

    ax_spread = ax_energy.twinx()
    ax_spread.plot(g.z, np.average(g.el_e_spread * m_e_GeV * 1000, weights=g.I, axis=0), 'm--', g.z, np.amax(g.el_e_spread * m_e_GeV * 1000, axis=0), 'r--', linewidth=1.5)
    ax_spread.set_ylabel(r'$\sigma_E$ [MeV]')
    ax_spread.grid(False)
    ax_spread.set_ylim(ymin=0)

    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spread.yaxis.major.locator.set_params(nbins=number_ticks)

    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')
    ax_spread.tick_params(axis='y', which='both', colors='r')
    ax_spread.yaxis.label.set_color('r')


def subfig_evo_el_bunching(ax_bunching, g, legend):
    number_ticks = 6
    
    ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I, axis=0), 'k-', g.z, np.amax(g.bunching, axis=0), 'grey', linewidth=1.5)
    # ax_bunching.plot(g.z, np.amax(g.bunching, axis=0), 'grey',linewidth=1.5) #only max
    ax_bunching.set_ylabel(r'Bunching')
    ax_bunching.set_ylim(ymin=0)
    # ax_bunching.set_ylim([0,0.8])
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.grid(True)


def subfig_evo_rad_pow_en(ax_rad_pow, g, legend, log=1):
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
            ax_rad_en.plot(g.z, g.pulse_energy*1e6, 'k--', linewidth=1.5)
            ax_rad_en.set_ylabel(r'E [$\mu$J]')
        else:
            ax_rad_en.plot(g.z, g.pulse_energy*1e3, 'k--', linewidth=1.5)
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

    ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W ' '\n' r'$E_{end}$= %.2e J' % (np.amax(g.p_int[:, -1]), np.mean(g.p_int[:, -1], axis=0) * g('xlamds') * g('zsep') * g.nSlices / speed_of_light), fontsize=12, horizontalalignment='right', verticalalignment='bottom', transform=ax_rad_pow.transAxes)


def subfig_evo_rad_pow(ax_rad_pow, g, legend, log=1):
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
    ax_rad_pow.text(0.98, 0.02, r'$P_{end}$= %.2e W' % (np.amax(g.p_int[:, -1])), fontsize=12, horizontalalignment='right', verticalalignment='bottom', transform=ax_rad_pow.transAxes)


def subfig_evo_rad_spec(ax_spectrum, g, legend, log=1):
    
    if 'spec' not in dir(g):
        g.calc_spec()
    
    ax_spectrum.plot(g.z, np.amax(g.spec, axis=0), 'r-', linewidth=1.5)
    ax_spectrum.text(0.5, 0.98, r"(on axis)", fontsize=10, horizontalalignment='center', verticalalignment='top', transform=ax_spectrum.transAxes)  # horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel(r'P$(\lambda)_{max}$ [a.u.]')
    plt.yticks(plt.yticks()[0][0:-1])
    
    if np.amax(np.amax(g.spec, axis=0)) > 0 and log:
        ax_spectrum.set_yscale('log')
    ax_spectrum.grid(True)
    
    spectrum_lamdwidth_fwhm = np.zeros_like(g.z)
    spectrum_lamdwidth_std = np.zeros_like(g.z)
    
    for zz in range(g.nZ):
        spectrum_lamdwidth_fwhm[zz] = None
        spectrum_lamdwidth_std[zz] = None
        if np.sum(g.spec[:,zz])!=0:
            pos, width, arr = fwhm3(g.spec[:, zz])
            if width != None:
                if arr[0] == arr[-1]:
                    dlambda = abs(g.freq_lamd[pos] - g.freq_lamd[pos-1])
                else:
                    dlambda = abs( (g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / (arr[0] - arr[-1]) )
                spectrum_lamdwidth_fwhm[zz] = dlambda * width / g.freq_lamd[pos]  
                # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[arr[0]] - g.freq_lamd[arr[-1]]) / g.freq_lamd[pos]  # the FWHM of spectral line (error when peakpos is at the edge of lamdscale)
            
            spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)

        # try:
            # peak = fwhm3(g.spec[:, zz])
            # spectrum_lamdwidth_fwhm[zz] = abs(g.freq_lamd[0] - g.freq_lamd[1]) * peak[1] / g.freq_lamd[peak[0]]  # the FWHM of spectral line (error when paekpos is at the edge of lamdscale)
        # except:
            # spectrum_lamdwidth_fwhm[zz] = 0

        # try:
            # spectrum_lamdwidth_std[zz] = std_moment(g.freq_lamd, g.spec[:, zz]) / n_moment(g.freq_lamd, g.spec[:, zz], 0, 1)
        # except:
            # spectrum_lamdwidth_std[zz] = 0

    ax_spec_bandw = ax_spectrum.twinx()
    ax_spec_bandw.plot(g.z, spectrum_lamdwidth_fwhm * 100, 'm--', label="fwhm")
    ax_spec_bandw.plot(g.z, 2*spectrum_lamdwidth_std * 100, 'm:', label="std")
    ax_spec_bandw.grid(False)
    plt.yticks(plt.yticks()[0][0:-1])
    ax_spec_bandw.set_ylim(ymin=0)
    
    if legend:
        ax_spec_bandw.legend()
    ax_spec_bandw.set_ylabel(r'$\Delta\lambda/\lambda, \%$'+'\n'+r'(-- fwhm, $\cdots2\sigma$)')
    if spectrum_lamdwidth_fwhm[-1] != None:
        ax_spec_bandw.text(0.98, 0.98, r"$(\Delta\lambda/\lambda)_{end}^{fwhm}$= %.2e"%(spectrum_lamdwidth_fwhm[-1]), fontsize=12, horizontalalignment='right', verticalalignment='top', transform=ax_spec_bandw.transAxes)

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
    


def subfig_rad_size(ax_size_t, g, legend):
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
    ax_size_t.set_ylabel(r'$\sim$size$_{transv}$ [$\mu$m]'+'\n'+u'(\u2014 $2\sigma$)')
    ax_size_t.grid(True)
    plt.yticks(plt.yticks()[0][0:-1])
    
    if g.nSlices > 1:
        ax_size_s = ax_size_t.twinx()
        size_long_fwhm = np.zeros_like(g.z)
        size_long_std = np.zeros_like(g.z)
        s = g.t * speed_of_light * 1.0e-15 * 1e6
        delta_s = (s[1] - s[0])
        for zz in range(g.nZ):
            #size_long_fwhm[zz] = fwhm(g.s,g.p_int[:, zz])
            if np.sum(g.p_int[:,zz])!=0:
                # try:
                _, width, _ = fwhm3(g.p_int[:, zz])
                if width != None:
                    size_long_fwhm[zz] = abs(delta_s) * width
                else:
                    size_long_fwhm[zz] = None
                # except:
                    # size_long_fwhm[zz] = 0
                
                # try:
                size_long_std[zz] = std_moment(s, g.p_int[:, zz])
                # except:
                    # size_long_std[zz] = 0
            else:
                size_long_fwhm[zz] = None
                size_long_std[zz] = None
            
        ax_size_s.plot(g.z, size_long_fwhm, color='navy', linestyle='--', linewidth=1.0, label="fwhm")
        ax_size_s.plot(g.z, 2*size_long_std, color='navy', linestyle=':', linewidth=1.0, label="std")
        ax_size_s.set_ylim(ymin=0)
        ax_size_s.set_ylabel(r'size$_{long}$ [$\mu$m]'+'\n'+r'(-- fwhm, $\cdots2\sigma$)')
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


def subfig_evo_rad_pow_sz(ax_power_evo, g, legend, norm=1):
    if g.nSlices > 1:
        z = g.z
        s = g.s
        power = g.p_int
        if norm == 1:
            max_power = np.max(power, 0)[np.newaxis, :]
            max_power[max_power == 0] = 1  # avoid division by zero
            power = power / max_power
            # power[isnan(power)]=0
        ax_power_evo.pcolormesh(z, s * 1e6, power)
        ax_power_evo.set_xlabel('z [m]')
        ax_power_evo.set_ylabel('s [$\mu$m]')
        ax_power_evo.axis('tight')
        ax_power_evo.grid(True)
    else:
        pass


def subfig_evo_rad_spec_sz(ax_spectrum_evo, g, legend, norm=1):
    if g.nSlices > 1:
        z = g.z
        l = g.freq_lamd
        spectrum = g.spec
        if norm == 1:
            max_spectrum = np.max(spectrum, 0)[np.newaxis, :]
            max_spectrum[max_spectrum == 0] = 1  # avoid division by zero
            spectrum = spectrum / max_spectrum
            # spectrum[isnan(spectrum)]=0
        ax_spectrum_evo.pcolormesh(z, l, spectrum)
        ax_spectrum_evo.set_xlabel('z [m]')
        ax_spectrum_evo.set_ylabel('$\lambda$ [nm]')
        ax_spectrum_evo.axis('tight')
        ax_spectrum_evo.grid(True)
    else:
        pass


def plot_gen_out_scanned_z(g, figsize=(10, 14), legend=True, fig_name=None, z=inf, savefig=False):

    if g('itdp') == True:
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
        if g.fileName() is '':
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

    ax_energy.plot(s, g.el_energy[:, zi] * m_e_GeV, 'b-', s, (g.el_energy[:, zi] + g.el_e_spread[:, zi]) * m_e_GeV, 'r--', s, (g.el_energy[:, zi] - g.el_e_spread[:, zi]) * m_e_GeV, 'r--')
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

def plot_dfl_all(dfl, **kwargs):
    
    plot_dfl(dfl, **kwargs)
    dfl.fft_z()
    plot_dfl(dfl, **kwargs)
    dfl.fft_xy()
    plot_dfl(dfl, **kwargs)
    dfl.fft_z()
    plot_dfl(dfl, **kwargs)
    dfl.fft_xy()

def plot_dfl(dfl, domains=None, z_lim=[], xy_lim=[], figsize=4, cmap=def_cmap, legend=True, phase=False, fig_name=None, auto_zoom=False, column_3d=True, savefig=False, showfig=True, return_proj=False, line_off_xy = True, log_scale=0, debug=1, cmin=0, vartype_dfl=np.complex64):
    '''
    Plots dfl radiation object in 3d.

    dfl is RadiationField() object
    z_lim sets the boundaries to CUT the dfl object in z to ranges of e.g. [2,5] um or nm depending on freq_domain=False of True
    xy_lim sets the boundaries to SCALE the dfl object in x and y to ranges of e.g. [2,5] um or urad depending on far_field=False of True
    figsize rescales the size of the figure
    legend not used yet
    phase can replace Z projection or spectrum with phase front distribution
    z dimentions correspondingly
    fig_name is the desired name of the output figure, would be used as suffix to the image filename if savefig==True
    auto_zoom automatically scales xyz the images to the (1%?) of the intensity limits
    column_3d plots top and side views of the radiation distribution
    savefig and showfig allow to save figure to image (savefig='png' (default) or savefig='eps', etc...) or to display it (slower)
    return_proj returns [xy_proj,yz_proj,xz_proj,x,y,z] array.
    vartype_dfl is the data type to store dfl in memory [either complex128 (two 64-bit floats) or complex64 (two 32-bit floats)], may save memory
    '''
    import matplotlib.colors as colors
    
    if showfig == False and savefig == False:
        return
    
    filePath = dfl.filePath

    text_present = 1
    _logger.info('plotting radiation field (dfl)')
    start_time = time.time()
    
    # print('dfl type is ',type(dfl))
    # if isinstance(dfl, RadiationField):
    # # if dfl.__class__ != RadiationField:
        # raise ValueError('wrong radiation object: should be RadiationField')
    
    dfl = deepcopy(dfl)
    
    if domains == None:
        domains = dfl.domains()
    else:
        dfldomain_check(domains)
        
    if 'k' in domains:
        far_field = True
    else:
        far_field = False
    if 'f' in domains:
        freq_domain = True
    else:
        freq_domain = False
    
    suffix = ''
   # if fig_name is None:
       # suffix = ''
   # else:
       # suffix = '_'+fig_name
        
    if dfl.Nz() != 1:
        # Make sure it is time-dependent
        ncar_z = dfl.Nz()
        leng_z = dfl.Lz()
        z = np.linspace(0, leng_z, ncar_z)
    else:
        column_3d = False
        phase = True
        freq_domain = False
        z_lim = []
    xlamds = dfl.xlamds

    # number of mesh points
    ncar_x = dfl.Nx()
    leng_x = dfl.Lx()  # transverse size of mesh [m]
    ncar_y = dfl.Ny()
    leng_y = dfl.Ly()
    E_pulse = dfl.E()

    if dfl.Nz() != 1:
        if freq_domain:
            if dfl.domain_z == 't':
                dfl.fft_z(debug=debug)
            
           # z = dfl.scale_z() * 1e9
           # dfl.fld = dfl.fld[::-1, :, :]
           # z = z[::-1]
           # unit_z = r'nm'
           # z_label = r'$\lambda$ [' + unit_z + ']'

            z = h_eV_s * speed_of_light / dfl.scale_z()
            unit_z = r'eV'
            z_label = r'$E_{{ph}}$ [{}]'.format(unit_z)
            
            z_labelv = r'[arb. units]'
            z_title = 'Spectrum'
            z_color = 'red'
            suffix += '_fd'
        else:
            if dfl.domain_z == 'f':
                dfl.fft_z(debug=debug)
            z = dfl.scale_z() * 1e6

            unit_z = r'$\mu$m'
            z_label = '$s$ [' + unit_z + ']'
            z_labelv = r'Power [W]'
            z_title = 'Z projection'
            z_color = 'blue'
    else:
        z = 0

    if z_lim != []:
        if len(z_lim) == 1:
            z_lim = [z_lim, z_lim]
        if z_lim[0] > z_lim[1]:
            z_lim[0] = -inf
            z_lim[1] = inf
        if z_lim[1] < np.amin(z) or z_lim[1] > np.amax(z):
            z_lim[1] = np.amax(z)
            # print('      set top lim to max')
        if z_lim[0] > np.amax(z) or z_lim[0] < np.amin(z):
            z_lim[0] = np.amin(z)
            # print('      set low lim to min')
        if debug > 1:
            print('      setting z-axis limits to ' + str(np.amin(z)) + ':' + str(z_lim[0]) + '-' + str(z_lim[1]) + ':' + str(np.amax(z)))  # tmp
        z_lim_1 = np.where(z <= z_lim[0])[0][-1]
        z_lim_2 = np.where(z >= z_lim[1])[0][0]

        if z_lim_1 == z_lim_2 and z_lim_1 == 0:
            z_lim_2 = z_lim_1 + 1
        elif z_lim_1 == z_lim_2 and z_lim_1 != 0:
            z_lim_1 = z_lim_2 - 1
        dfl.fld = dfl.fld[z_lim_1:z_lim_2, :, :]
        z = z[z_lim_1:z_lim_2]
        ncar_z = dfl.shape[0]
        suffix += '_zoom_%.2f-%.2f' % (np.amin(z), np.amax(z))

    if far_field:
        if dfl.domain_xy == 's':
            dfl.fft_xy(debug=debug)
        x = dfl.scale_x() * 1e6
        y = dfl.scale_y() * 1e6

        unit_xy = r'$\mu$rad'
        x_label = r'$\theta_x$ [' + unit_xy + ']'
        y_label = r'$\theta_y$ [' + unit_xy + ']'
        suffix += '_ff'
        x_title = 'X divergence'
        y_title = 'Y divergence'
        xy_title = 'Far field intensity'
        x_y_color = 'green'
        # if debug>1: print('        done in %.2f seconds' %(time.time()-calc_time))
    else:
        if dfl.domain_xy == 'k':
            dfl.fft_xy(debug=debug)
        x = dfl.scale_x() * 1e6
        y = dfl.scale_y() * 1e6
        
        unit_xy = r'$\mu$m'
        x_label = 'x [' + unit_xy + ']'
        y_label = 'y [' + unit_xy + ']'
        x_title = 'X projection'
        y_title = 'Y projection'
        xy_title = 'Intensity'
        x_y_color = 'blue'
    
    
    dfl.fld = dfl.fld.astype(np.complex64)
    xy_proj = dfl.int_xy()
    xy_proj_ph = np.angle(np.sum(dfl.fld,axis=0))  # tmp  # tmp
    yz_proj = dfl.int_zy()
    xz_proj = dfl.int_zx()
    z_proj = dfl.int_z()
    
    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])
    
    if log_scale:
        suffix += '_log'

    if fig_name is None:
        if dfl.fileName() is '':
            fig = plt.figure('Radiation distribution' + suffix)
        else:
            fig = plt.figure('Radiation distribution' + suffix + ' ' + dfl.fileName())
    else:
        fig = plt.figure(fig_name + suffix)
    del dfl

    fig.clf()
    fig.set_size_inches(((3 + 2 * column_3d) * figsize, 3 * figsize), forward=True)

    # cmap = plt.get_cmap(def_cmap)  # jet inferno viridis #change to convenient
    cmap_ph = plt.get_cmap('hsv')
    
    if line_off_xy:
        x_line = xy_proj[int((ncar_y - 1) / 2), :]
        y_line = xy_proj[:, int((ncar_x - 1) / 2)]
        x_title += ' lineoff'
        y_title += ' lineoff'
    else:
        x_line = np.sum(xy_proj, axis=0)
        y_line = np.sum(xy_proj, axis=1)
    
    if np.max(x_line) != 0 and np.max(y_line) != 0:
        x_line, y_line = x_line / np.max(x_line), y_line / np.max(y_line)
        
        
    if cmin not in [None, False, 0]:
        cmap = matplotlib.cm.get_cmap(cmap)
        cmap.set_under("w")
        xy_proj[xy_proj < xy_proj.max() * cmin] = -1e-10
        yz_proj[yz_proj < yz_proj.max() * cmin] = -1e-10
        xz_proj[xz_proj < xz_proj.max() * cmin] = -1e-10

    if log_scale:
        xy_proj[xy_proj <= 0] = None
        yz_proj[yz_proj <= 0] = None
        xz_proj[xz_proj <= 0] = None
        z_proj[z_proj <= 0] = None

    ax_int = fig.add_subplot(2, 2 + column_3d, 1)
    if log_scale:
        intplt = ax_int.pcolormesh(x, y, xy_proj, norm=colors.LogNorm(vmin=xy_proj.min(), vmax=xy_proj.max()), cmap=cmap)
    else:
        intplt = ax_int.pcolormesh(x, y, xy_proj, cmap=cmap, vmin=0)
    ax_int.set_title(xy_title, fontsize=15)
    ax_int.set_xlabel(r'' + x_label)
    ax_int.set_ylabel(y_label)
    if np.size(z) > 1 and text_present:
        ax_int.text(0.01, 0.01, r'$E_{p}$=%.2e J' % (E_pulse), horizontalalignment='left', verticalalignment='bottom', fontsize=12, color='white', transform=ax_int.transAxes)

    if phase == True:
        ax_ph = fig.add_subplot(2, 2 + column_3d, 4 + column_3d, sharex=ax_int, sharey=ax_int)
        ax_ph.pcolormesh(x, y, xy_proj_ph, cmap=cmap_ph)
        ax_ph.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
        ax_ph.set_title('Phase', fontsize=15)
    else:
        ax_z = fig.add_subplot(2, 2 + column_3d, 4 + column_3d)
        if log_scale:
            ax_z.semilogy(z, z_proj, linewidth=1.5, color=z_color)
        else:
            ax_z.plot(z, z_proj, linewidth=1.5, color=z_color)
        ax_z.set_title(z_title, fontsize=15)
        ax_z.set_xlabel(z_label)
        ax_z.set_ylabel(z_labelv)
        ax_z.set_ylim(ymin=0)

    ax_proj_x = fig.add_subplot(2, 2 + column_3d, 3 + column_3d, sharex=ax_int)
    ax_proj_x.set_title(x_title, fontsize=15)

    if sum(x_line) != 0:
        try:
            x_line_f, rms_x = gauss_fit(x, x_line)  # fit with Gaussian, and return fitted function and rms
        except RuntimeWarning:
            x_line_f = np.zeros_like(x_line)
            rms_x = 0
        try:
            fwhm_x = fwhm3(x_line)[1] * dx  # measure FWHM
        except ValueError:
            fwhm_x = 0
    else:
        x_line_f = np.zeros_like(x_line)
        rms_x = 0
        fwhm_x = 0
    
    
    
    
    if log_scale:
        ax_proj_x.semilogy(x, x_line, linewidth=2, color=x_y_color)
        ax_proj_x.semilogy(x, x_line_f, color='grey')
        ax_proj_x.set_ylim(ymin = np.amin(x_line), ymax=1)
    else:
        ax_proj_x.plot(x, x_line, linewidth=2, color=x_y_color)
        ax_proj_x.plot(x, x_line_f, color='grey')
        ax_proj_x.set_ylim(ymin=0, ymax=1)
    
    if text_present:
        try:
            ax_proj_x.text(0.95, 0.95, 'fwhm= \n' + str(round_sig(fwhm_x, 3)) + r' [' + unit_xy + ']\nrms= \n' + str(round_sig(rms_x, 3)) + r' [' + unit_xy + ']', horizontalalignment='right', verticalalignment='top', transform=ax_proj_x.transAxes, fontsize=12)
        except:
            pass
    
    ax_proj_x.set_xlabel(x_label)

    ax_proj_y = fig.add_subplot(2, 2 + column_3d, 2, sharey=ax_int)
    ax_proj_y.set_title(y_title, fontsize=15)

    if sum(y_line) != 0:
        try:
            y_line_f, rms_y = gauss_fit(y, y_line)  # fit with Gaussian, and return fitted function and rms
        except RuntimeWarning:
            y_line_f = np.zeros_like(y_line)
            rms_y = 0
        try:
            fwhm_y = fwhm3(y_line)[1] * dy  # measure FWHM
        except ValueError:
            fwhm_y = 0
    else:
        y_line_f = np.zeros_like(y_line)
        rms_y = 0
        fwhm_y = 0
    
    if log_scale:
        ax_proj_y.semilogx(y_line, y, linewidth=2, color=x_y_color)
        ax_proj_y.semilogx(y_line_f, y, color='grey')
        ax_proj_y.set_xlim(xmin=np.amin(y_line), xmax=1)
    else:
        ax_proj_y.plot(y_line, y, linewidth=2, color=x_y_color)
        ax_proj_y.plot(y_line_f, y, color='grey')
        ax_proj_y.set_xlim(xmin=0, xmax=1)
    
    if text_present:
        try:
            ax_proj_y.text(0.95, 0.95, 'fwhm= ' + str(round_sig(fwhm_y, 3)) + r' [' + unit_xy + ']\nrms= ' + str(round_sig(rms_y, 3)) + r' [' + unit_xy + ']', horizontalalignment='right', verticalalignment='top', transform=ax_proj_y.transAxes, fontsize=12)
        except:
            pass
    
    ax_proj_y.set_ylabel(y_label)
    
    # if log_scale:
        # ax_proj_x.set_yscale('log')
        # ax_proj_y.set_xscale('log')
        # if not phase:
            # ax_z.set_yscale('log')

    if column_3d:
        
        if log_scale:
            cut_off = 1e-6
            yz_proj[yz_proj < yz_proj.max() * cut_off] = 0
            xz_proj[xz_proj < xz_proj.max() * cut_off] = 0
            # cut-off = np.amin([yz_proj[yz_proj!=0].min(), xz_proj[xz_proj!=0].min()]) / 10
            # yz_proj += minmin
            # xz_proj += minmin
            min_xz_proj=xz_proj[xz_proj!=0].min()
            min_yz_proj=yz_proj[yz_proj!=0].min()
            
        
        # if np.amin(xz_proj) == 0:
            # min_xz_proj = 0
        # else:
            # min_xz_proj=xz_proj[xz_proj!=0].min()
        # if np.amin(yz_proj) == 0:
            # min_yz_proj = 0
        # else:
            # min_yz_proj=yz_proj[yz_proj!=0].min()
        
        if phase == True:
            ax_proj_xz = fig.add_subplot(2, 2 + column_3d, 6)
        else:
            ax_proj_xz = fig.add_subplot(2, 2 + column_3d, 6, sharex=ax_z)
        if log_scale:
            ax_proj_xz.pcolormesh(z, x, np.swapaxes(xz_proj, 1, 0), norm=colors.LogNorm(vmin=min_xz_proj, vmax=xz_proj.max()), cmap=cmap)
        else:
            ax_proj_xz.pcolormesh(z, x, np.swapaxes(xz_proj, 1, 0), cmap=cmap, vmin=0)
        ax_proj_xz.set_title('Top view', fontsize=15)
        ax_proj_xz.set_xlabel(z_label)
        ax_proj_xz.set_ylabel(x_label)
        

        ax_proj_yz = fig.add_subplot(2, 2 + column_3d, 3, sharey=ax_int, sharex=ax_proj_xz)
        if log_scale:
            ax_proj_yz.pcolormesh(z, y, np.swapaxes(yz_proj, 1, 0), norm=colors.LogNorm(vmin=min_yz_proj, vmax=yz_proj.max()), cmap=cmap)
        else:
            ax_proj_yz.pcolormesh(z, y, np.swapaxes(yz_proj, 1, 0), cmap=cmap, vmin=0)
        ax_proj_yz.set_title('Side view', fontsize=15)
        ax_proj_yz.set_xlabel(z_label)
        ax_proj_yz.set_ylabel(y_label)

    cbar = 0
    if cbar:
        fig.subplots_adjust(top=0.95, bottom=0.05, right=0.85, left=0.1)
        cbar_int = fig.add_axes([0.89, 0.15, 0.015, 0.7])
        cbar = plt.colorbar(intplt, cax=cbar_int)  # pad = -0.05 ,fraction=0.01)
        # cbar.set_label(r'[$ph/cm^2$]',size=10)
        cbar.set_label(r'a.u.', size=10)

    if auto_zoom != False:
        size_x = np.max(abs(x[np.nonzero(x_line > 0.005)][[0, -1]]))
        size_y = np.max(abs(x[np.nonzero(x_line > 0.005)][[0, -1]]))
        size_xy = np.max([size_x, size_y])
        # print(size_xy) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        # print(zlim_calc) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if z_lim == []:
            zlim_calc = z[np.nonzero(z_proj > np.max(z_proj) * 0.005)][[0, -1]]
            if column_3d == True:
                ax_proj_xz.set_xlim(zlim_calc)
                ax_proj_xz.set_ylim([-size_xy, size_xy])
            if phase == False:
                # _logger.debug(ind_str + 'scaling z to {:}'.format(zlim_calc))
                ax_z.set_xlim(zlim_calc)
        # elif phase == False and z_lim == []:
            # ax_z.set_xlim(zlim_calc)
            # _logger.debug(ind_str + 'scaling xy to {:}'.format(size_xy))
        # elif column_3d == True:
        
        # ax_int.axis('equal')
        ax_int.axis([-size_xy, size_xy, -size_xy, size_xy])
        suffix += '_zmd'
    else:
        if column_3d == True:
            ax_proj_xz.axis('tight')
            ax_proj_yz.axis('tight')
        elif column_3d == False and phase == False:
            ax_z.axis('tight')
        # ax_int.set_aspect('equal')
        ax_int.autoscale(tight=True)

    if len(xy_lim) == 2:
        ax_int.axis([-xy_lim[0], xy_lim[0], -xy_lim[1], xy_lim[1]])
        ax_proj_xz.set_ylim([-xy_lim[0], xy_lim[0]])
    elif len(xy_lim) == 1:
        ax_int.axis([-xy_lim[0], xy_lim[0], -xy_lim[0], xy_lim[0]])
        ax_proj_xz.set_ylim([-xy_lim[0], xy_lim[0]])

    fig.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.draw()

    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving *{:}.{:}'.format(suffix,savefig))
        fig.savefig(filePath + suffix + '.' + str(savefig), format=savefig)

    _logger.info(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))

    plt.draw()
    if showfig == True:
        _logger.debug(ind_str + 'showing dfl')
        rcParams["savefig.directory"] = os.path.dirname(filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)

    if return_proj:
        return [xy_proj, yz_proj, xz_proj, x, y, z]
    else:
        return


def plot_gen_stat(proj_dir, run_inp=[], stage_inp=[], param_inp=[], s_param_inp=['p_int', 'pulse_energy', 'r_size_weighted', 'spec', 'spec_phot_density', 'error'], z_param_inp=['p_int', 'phi_mid_disp', 'spec', 'spec_phot_density', 'bunching', 'wigner'], dfl_param_inp=['dfl_spec'], run_param_inp=['p_int', 'spec', 'spec_phot_density', 'pulse_energy'], s_inp=['max'], z_inp=[0,'end'], run_s_inp=['max'], run_z_inp=['end'], spec_pad=1, savefig=1, saveval=1, showfig=0, debug=1):
    '''
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
    '''
    import copy
    rc('text', usetex=False)
    dict_name = {'p_int': 'radiation power', 'pulse_energy': 'radiation pulse energy', 'el_e_spread': 'el.beam energy spread', 'el_energy': 'el.beam energy average', 'bunching': 'el.beam bunching', 'spec': 'radiation on-axis spectral density', 'spec_phot_density': 'radiation spectral photon density', 'dfl_spec': 'total radiation photon spectral density (dfl)', 'r_size': 'radiation transv size', 'r_size_weighted': 'radiation transv size (weighted)', 'xrms': 'el.beam x size', 'yrms': 'el.beam y size', 'error': 'genesis simulation error', 'p_mid': 'radiation power on-axis', 'phi_mid': 'radiation phase on-axis', 'increment': 'radiation power increment'}
    dict_unit = {'p_int': '[W]', 'pulse_energy': '[J]', 'el_e_spread': '(gamma)', 'el_energy': '(gamma)', 'bunching': '', 'spec': '[arb.units]', 'spec_phot_density': '(estimation) [ph/eV]', 'dfl_spec': '[ph/eV]', 'r_size': '[m]', 'xrms': '[m]', 'yrms': '[m]', 'error': ''}

    figsize = (14, 7)
    figsize = (8, 6)

    if debug > 0:
        print ('statistical postprocessing started')
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
                    outlist[irun].calc_spec(npad = spec_pad)
                    run_range_good.append(irun)
                except:
                    print('     could not read '+out_file)

        run_range = run_range_good
        
        if len(run_range_good) == 0:
            continue

        nSlices = np.array([int(outlist[run].nSlices) for run in run_range])
        nZ = np.array([int(outlist[run].nZ) for run in run_range])

        f_nSlices = np.argmax(np.bincount(nSlices)) #most abundant number of slices
        f_nZ = np.argmax(np.bincount(nZ)) #most abundant number of z records
        
        index = np.where(np.logical_and(nSlices==f_nSlices, nZ==f_nZ))[0]
        run_range_good = [run_range[i] for i in index]
        if len(list(set(run_range)-set(run_range_good))) > 0:
            print('run_range', run_range)
            print('run_range_good', run_range_good)
            print('discarding runs',  list(set(run_range)-set(run_range_good)))
        
        run_range = run_range_good
        
        # if len(run_range)!=0 and debug>0:
            # print('stage = ', stage)

        # check if all gout have the same number of slices nSlice and history records nZ
        # for irun in run_range[1:]:
            # if outlist[irun].nSlices != outlist[run_range[0]].nSlices or outlist[irun].nZ != outlist[run_range[0]].nZ:
                # raise ValueError('Non-uniform out objects')

        if run_range == [] or len(run_range) == 1:
            continue

        if debug > 0:
            print('    processing runs ' + str(run_range) + ' of stage ' + str(stage))

       # for irun in run_range:
           # out_file=proj_dir+'run_'+str(irun)+'/run.'+str(irun)+'.s'+str(stage)+'.gout'
           # outlist[irun] = read_out_file(out_file,read_level=1)
           # print(outlist[irun].sliceKeys)

        # if param_inp==[]:
           # if debug>1: print(outlist[run_range[0]].sliceKeys_used)
           # param_range=outlist[run_range[0]].sliceKeys_used
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
            # s_param_range=param_range
        # else:
        

        
        s_param_range = s_param_inp
        if debug > 0:
            print('    processing S parameters ' + str(s_param_range))
        if debug > 1:
            print('      s_inp ' + str(s_inp))

        for param in s_param_range:
            for s_ind in s_inp:
                s_value = []
                s_fig_name = 'stage_' + str(stage) + '__Z__' + dict_name.get(param, param).replace(' ', '_').replace('.', '_') + '__' + str(s_ind)
                for irun in run_range:
                    if not hasattr(outlist[irun], param):
                        continue
                    else:
                        if debug > 1:
                            print ('parameter = %s; s = %s; run = %s;' %(param, s_ind, irun))
                        param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                    if debug > 1:
                        print('shape param_matrix', shape(param_matrix))
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
                        print('plotting array shapes', shape(outlist[irun].z), shape(np.swapaxes(s_value, 0, 1)))
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
                        np.savetxt(saving_path + s_fig_name + '.txt', np.vstack([outlist[irun].z, np.mean(s_value, 0), s_value]).T, fmt="%E", newline='\n', comments='')
                    if not showfig:
                        plt.close('all')
        # if z_param_inp==[]:
            # z_param_range=param_range
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
                    out=outlist[irun]
                    W=wigner_out(out, z=z_ind, debug=0, pad=wig_pad)
                    w += W.wig
                W.wig= w / len(outlist)

                W.filePath = proj_dir + 'results' + os.path.sep + 'stage_' + str(stage) + '__WIG__' + str(z_ind) + '__m'
                wig_fig_name = 'stage_' + str(stage) + '__WIG__' + str(z_ind) + '__m'
                plot_wigner(W, z=z_ind, x_units='um', y_units='ev', fig_name=wig_fig_name, savefig=savefig, showfig=showfig, debug=0)
                if saveval != False:
                    if debug > 1:
                        print('      saving ' + wig_fig_name + '.txt')
                    np.savetxt(saving_path + wig_fig_name + '.txt', W.wig, fmt='%E ', newline='\n')
                    np.savetxt(saving_path + wig_fig_name + '_sc.txt', np.vstack([speed_of_light*h_eV_s*1e9/W.freq_lamd, W.s]).T, fmt='%E ', newline='\n', header=' E[eV], s[m]')
     

        for param in z_param_range:
            for z_ind in z_inp:
                z_value = []
                z_fig_name = 'stage_' + str(stage) + '__S__' + dict_name.get(param, param).replace(' ', '_').replace('.', '_') + '__' + str(z_ind) + '__m'
                for irun in run_range:
                    if not hasattr(outlist[irun], param):
                        break
                    else:
                        if debug > 1:
                            print ('parameter = %s; z = %s; run = %s;' %(param, z_ind, irun))
                        param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                    if debug > 1:
                        print('shape param_matrix', shape(param_matrix))
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
                            print('plotting array shapes freq', shape(freq_scale), shape(np.swapaxes(z_value, 0, 1)))
                        fig = plt.plot(freq_scale, np.swapaxes(z_value, 0, 1), '0.8')
                        fig = plt.plot(freq_scale, z_value[0], '0.5', linewidth=1)
                        fig = plt.plot(freq_scale, np.mean(z_value, 0), 'k', linewidth=2)
                        plt.xlim([np.min(freq_scale), np.max(freq_scale)])
                        plt.xlabel('$E_{photon}$ [eV]')
                    else:
                        s_scale = outlist[irun].s * 1e6
                        if debug > 1:
                            print('plotting array shapes', shape(s_scale), shape(np.swapaxes(z_value, 0, 1)))
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
                            np.savetxt(saving_path + z_fig_name + '.txt', np.vstack([outlist[irun].freq_ev, np.mean(z_value, 0), z_value]).T, fmt="%E", newline='\n', comments='')
                        else:
                            np.savetxt(saving_path + z_fig_name + '.txt', np.vstack([outlist[irun].s * 1e6, np.mean(z_value, 0), z_value]).T, fmt="%E", newline='\n', comments='')
                    if not showfig:
                        plt.close('all')
        # if run_param_inp==[]:
            # run_param_range=[]
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
                    run_fig_name = 'stage_' + str(stage) + '__RUN__' + dict_name.get(param, param).replace(' ', '_').replace('.', '_') + '__' + str(s_ind) + '__um__' + str(z_ind) + '__m'
                    for irun in run_range:
                        if not hasattr(outlist[irun], param):
                            break
                        else:
                            if debug > 1:
                                print ('parameter = %s; z = %s; s = %s; run = %s' %(param, z_ind, s_ind, irun))
                            param_matrix = copy.deepcopy(getattr(outlist[irun], param))
                            if debug > 1:
                                print('shape param_matrix', shape(param_matrix))
                            if debug > 1:
                                print('length', len(param_matrix), len(outlist[irun].z))

                        if len(param_matrix) != len(outlist[irun].z):  # case if the array is 1D (no s/z matrix presented)
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

                        plt.ylabel(dict_name.get(param, param) + '  ' + dict_unit.get(param, '') + ' (' + str(s_ind) + ' um, ' + str(z_ind) + ' m)')
                        if savefig != False:
                            if debug > 1:
                                print('      saving ' + run_fig_name + '.' + savefig)
                            plt.draw()
                            plt.savefig(saving_path + run_fig_name + '.' + savefig, format=savefig)
                        if saveval != False:
                            if debug > 1:
                                print('      saving ' + run_fig_name + '.txt')
                            np.savetxt(saving_path + run_fig_name + '.txt', np.vstack([run_range, run_value_arr]).T, fmt="%E", newline='\n', comments='')
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
                        np.savetxt(saving_path + dfl_fig_name + '.txt', np.vstack([freq_scale, np.mean(dfl_value, 0), dfl_value]).T, fmt="%E", newline='\n', comments='')
                if not showfig:
                    plt.close('all')
    if showfig:
        plt.draw()
        plt.show()
    else:
        plt.close('all')

    if debug > 0:
        print('done in %.2f seconds' % (time.time() - start_time))


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

    corr_fig_name = 'corr_' + param_1 + '_s' + str(stage_1) + '_at' + str(z_1) + '_' + str(s_1) + '__' + param_2 + '_s' + str(stage_2) + '_at' + str(z_2) + '_' + str(s_2)

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
    if savefig != False:
        print('      saving ' + corr_fig_name + '.' + savefig)
        plt.savefig(saving_path + corr_fig_name + '.' + savefig, format=savefig)
    if saveval != False:
        print('      saving ' + corr_fig_name + '.txt')
        np.savetxt(saving_path + corr_fig_name + '.txt', np.vstack([var_1, var_2]).T, fmt="%E", newline='\n', comments=param_1 + '_s' + str(stage_1) + '_at' + str(z_1) + '_' + str(s_1) + ' ' + param_2 + '_s' + str(stage_2) + '_at' + str(z_2) + '_' + str(s_2))

    if showfig:
        plt.show()
    else:
        plt.close('all')

    return fig


def plot_dpa_bucket_out(out, dpa=None, slice_pos='max_I', repeat=1, GeV=1, figsize=4, cmap=def_cmap, scatter=True, energy_mean=None, legend=True, fig_name=None, savefig=False, showfig=True, bins=[50,50], debug=1):
    
    if dpa == None:
        dpa=read_dpa_file_out(out)
        
    
    if out.nSlices > 1:
        if type(slice_pos) == str:
            if slice_pos == 'max_I':
                slice_num = np.argmax(out.I)
            elif slice_pos == 'max_P':
                slice_num = np.argmax(out.power)
            elif slice_pos == 'max_B':
                slice_num = np.argmax(out.bunching[:,-1])
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
    suffix = '_%.2fum_%2.fm' % (slice_pos_act*1e6,np.amax(out.z))
    if scatter: suffix = '_scatter' + suffix
    return plot_dpa_bucket(dpa=dpa, slice_num=slice_num, repeat=repeat, GeV=GeV, figsize=figsize, cmap=cmap, scatter=scatter, energy_mean=energy_mean, legend=legend, fig_name=fig_name, savefig=savefig, showfig=showfig, suffix=suffix, bins=bins, debug=debug)


def plot_dpa_bucket(dpa, slice_num=None, repeat=1, GeV=1, figsize=4, cmap=def_cmap, scatter=False, energy_mean=None, legend=True, fig_name=None, savefig=False, showfig=True, suffix='', bins=(50,50), debug=1, return_mode_gamma=0):
    part_colors = ['darkred', 'orange', 'g', 'b', 'm', 'c', 'y']
    # cmap='BuPu'
    y_bins = bins[0]
    z_bins = bins[1]
    
    if showfig == False and savefig == False:
        return

    _logger.info('plotting dpa bucket')
    start_time = time.time()
    
    if dpa.__class__ != GenesisParticlesDump:
        raise ValueError('wrong particle object: should be GenesisParticlesDump')

    if shape(dpa.ph)[0] == 1:
        slice_num = 0
    if slice_num is None:
        slice_num = int(shape(dpa.ph)[0]/2)
        _logger.debug(ind_str + 'no slice number provided, using middle of the distribution - slice number {}'.format(slice_num))
    else:
        assert (slice_num <= shape(dpa.ph)[0]), 'slice_num larger than the dpa shape'

    if fig_name == None:
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

    nbins = shape(dpa.ph)[1]
    phase = deepcopy(dpa.ph[slice_num, :, :])
    energy = deepcopy(dpa.e[slice_num, :, :])
    _logger.debug(ind_str + 'nbins =  {}'.format(nbins))



    if GeV:
        energy *= m_e_MeV
        if energy_mean == None: 
            energy_mean = round(np.mean(energy), 0)
    else:
        if energy_mean == None: 
            energy_mean = round(np.mean(energy), 1)
    energy -= energy_mean

    phase_flat=phase.flatten()
    energy_flat=energy.flatten()
    for irep in range(repeat-1):
        phase_flat=np.append(phase_flat,phase.flatten() + 2 * np.pi * (irep+1))
        energy_flat=np.append(energy_flat,energy.flatten())
    
    # phase_hist = np.ravel(phase)
    # for irep in range(repeat-1):
        # phase_hist = np.concatenate((phase_hist, np.ravel(phase) + 2 * np.pi * (irep+1)))

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
        

    if scatter == True:
        for irep in range(repeat):
            for ibin in range(nbins):
                ax_main.scatter(phase[ibin, :] + 2 * np.pi * (irep), energy[ibin, :], color=part_colors[ibin], marker='.')

        # ax_z_hist.set_xlim([edges[0], edges[-1]])

    elif scatter == False:
        ax_main.hist2d(phase_flat, energy_flat, bins=[z_bins * repeat, y_bins], cmin=0, cmap=cmap)
        
    
                
    ax_main.set_xlabel('$\phi$ [rad]')
    if GeV:
        ax_main.set_ylabel('E [MeV] + ' + str(energy_mean / 1000) + ' [GeV]')
    else:
        ax_main.set_ylabel('$\gamma$ + ' + str(energy_mean))

    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving to {}'.format(dpa.fileName() + suffix + '.' + savefig))
        plt.savefig(dpa.filePath + suffix + '.' + savefig, format=savefig)

    if showfig:
        rcParams["savefig.directory"] = os.path.dirname(dpa.filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


def plot_edist(edist, figsize=4, fig_name=None, savefig=False, showfig=True, scatter=False, plot_x_y=True, plot_xy_s=True, bins=(50, 50, 50, 50), flip_t=False, x_units='um', y_units='ev', cmin=0, y_offset=None, cmap=def_cmap, debug=1):

    if showfig == False and savefig == False:
        return
    _logger.info('plotting edist file')
    start_time = time.time()
    # suffix=''
    # if edist.__class__ != GenesisElectronDist:
        # raise ValueError('wrong distribution object: should be GenesisElectronDist')

    if np.size(bins) == 1:
        bins = (bins, bins, bins, bins)  # x,y,t,e

    if fig_name == None:
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
        # s = -edist.t * speed_of_light * 1e6
    # else:
        # s = edist.t * speed_of_light * 1e6

    hist, edges = np.histogram(s, bins=bins[2])  # calculate current histogram
    edges = edges[0:-1]  # remove the last bin edge to save equal number of points
    hist_int = np.trapz(hist, edges) / mult  # normalize
    hist = np.rint(hist.astype(float) / (hist_int / float(edist.charge())))

    ax_curr = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 1)
    #ax_curr.hist(s, bins,color='b')
    ax_curr.plot(edges, hist/1000, color='k',linewidth=2)
    ax_curr.set_xlabel(s_label)
    ax_curr.set_ylabel('I [kA]')

    ax_se = fig.add_subplot(2, 1 + plot_x_y + plot_xy_s, 2 + plot_x_y + plot_xy_s, sharex=ax_curr)
    if y_units in ['ev', 'eV']:
        energy = edist.g * m_e_MeV
    else:  # elif beam_E_plot=='gamma':
        energy = edist.g
        
    if y_offset == None:
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
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving ' + edist.filePath + '.' + savefig)
        plt.savefig(edist.filePath + '.' + savefig, format=savefig)

    if showfig:
        rcParams["savefig.directory"] = os.path.dirname(edist.filePath)
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)

    _logger.info(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
    # return fig


def plot_beam(beam, figsize=3, showfig=True, savefig=False, fig=None, plot_xy=None, debug=0):

    _logger.info('plotting beam')
    
    if showfig == False and savefig == False:
        return
    
    # if beam.__class__ != GenesisBeam:
        # raise ValueError('wrong beam object: should be GenesisBeam')
    
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
    
    g0 = np.mean(beam.g).astype(int) #mean
    g_dev = beam.g - g0 #deviation from mean

    fig.set_size_inches((4 * figsize, (3 + plot_xy) * figsize), forward=True)
    ax = fig.add_subplot(2 + plot_xy, 2, 1)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.I, 'r', lw=3)
    
    ax.set_ylim(ymin=0)
    
    if hasattr(beam,'eloss'):
        if (beam.eloss!=0).any():
            ax = ax.twinx()
            p2, = plt.plot(1.e6 * np.array(beam.s), 1.e-3 * np.array(beam.eloss), 'g', lw=3)
            ax.legend([p1, p2], [r'$I [A]$', r'Wake $[KV/m]$'], fontsize=fontsize, loc='best')
        else:
            ax.legend([r'$I [A]$'], fontsize=fontsize, loc='best')
    else:
        ax.legend([r'$I [A]$'], fontsize=fontsize, loc='best')
    plt.plot(1.e6 * beam.s[idx], beam.I[idx], 'bs')
    # ax.set_xlim([np.amin(beam.s),np.amax(beam.x)])
    ax = fig.add_subplot(2 + plot_xy, 2, 2, sharex=ax)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
    #p1,= plt.plot(1.e6 * np.array(beam.s),1.e-3 * np.array(beam.eloss),'r',lw=3)
    p1, = plt.plot(1.e6 * np.array(beam.s), g_dev, 'r', lw=3)
    plt.plot(1.e6 * beam.s[idx], g_dev[idx], 'bs')
    ax = ax.twinx()
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.dg, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.dg[idx], 'bs')
    ax.legend([p1, p2], [r'$\gamma$ + '+str(g0), r'$\delta \gamma$'], loc='best')

    ax = fig.add_subplot(2 + plot_xy, 2, 3, sharex=ax)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.emit_xn * 1e6, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.emit_yn * 1e6, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.emit_xn[idx] * 1e6, 'bs')
    ax.set_ylim(ymin=0)

    ax.legend([p1, p2], [r'$\varepsilon_x [\mu m]$', r'$\varepsilon_y [\mu m]$'], fontsize=fontsize, loc='best')
    # ax3.legend([p3,p4],[r'$\varepsilon_x$',r'$\varepsilon_y$'])

    ax = fig.add_subplot(2 + plot_xy, 2, 4, sharex=ax)
    plt.grid(True)
    ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
    p1, = plt.plot(1.e6 * np.array(beam.s), beam.beta_x, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.s), beam.beta_y, 'g', lw=3)
    plt.plot(1.e6 * beam.s[idx], beam.beta_x[idx], 'bs')
    ax.set_ylim(ymin=0)
    ax.legend([p1, p2], [r'$\beta_x [m]$', r'$\beta_y [m]$'], fontsize=fontsize, loc='best')

    if plot_xy:

        ax = fig.add_subplot(3, 2, 5, sharex=ax)
        plt.grid(True)
        ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
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
        ax.set_xlabel(r'$s [\mu m]$',fontsize=fontsize)
        p1, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.xp), 'r', lw=3)
        p2, = plt.plot(1.e6 * np.array(beam.s), 1.e6 * np.array(beam.yp), 'g', lw=3)

        ax.legend([p1, p2], [r'$x_p [\mu rad]$', r'$y_p [\mu rad]$'], fontsize=fontsize, loc='best')

    ax.set_xlim([1.e6 * np.amin(beam.s), 1e6 * np.amax(beam.s)])

    fig.subplots_adjust(hspace=0.2, wspace=0.3)

    # if savefig != False:
        # if hasattr(beam,'filePath'):
            # if savefig == True:
                # filetype = 'png'
            # else:
            #
            # path = beam.filePath
            # name = beam.fileName()
        # else:
            # if if savefig == True:
            
            # if debug > 1:
                # print('      saving ' + beam.fileName() + '.' + savefig)
            # plt.savefig(beam.filePath + '.' + savefig, format=savefig)
        
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

def plot_wigner(wig_or_out, z=np.inf, x_units='um', y_units='ev', x_lim=(None,None), y_lim=(None,None), downsample=1, autoscale=None, figsize=4, cmap='seismic', abs_value=0, fig_name=None, savefig=False, showfig=True, plot_proj=1, plot_text=1, plot_moments=0, debug=1):
    '''
    plots wigner distribution (WD) with marginals
    wig_or_out -  may be WignerDistribution() or GenesisOutput() object
    z - (if isinstance(wig_or_out, GenesisOutput)) location at which WD will be calculated
    x_units - (um or fs) - units to display power scale
    y_units - (nm or eV) - units to display spectrum scale
    x_lim, y_lim - scaling limits in given units, (min,max) or [min,max], e.g: (None,6)
    abs_value - if True, absolute value of WD is displayed (usually, it has both positive and negative values)
    cmap - colormar if abs_value==False (http://matplotlib.org/users/colormaps.html)
    '''
    if showfig == False and savefig == False:
        return
    
    _logger.info('plotting Wigner distribution')
        
    if isinstance(wig_or_out, GenesisOutput):
        W=wigner_out(wig_or_out,z)
    elif isinstance(wig_or_out, WignerDistribution):
        W=wig_or_out
    else:
        raise ValueError('Unknown object for Wigner plot')
    
    
    if fig_name is None:
        if W.fileName() is '':
            fig_text = 'Wigner distribution'
        else:
            fig_text = 'Wigner distribution ' + W.fileName()
    else:
        fig_text = fig_name
    if W.z!=None:
        fig_text += ' ' + str(W.z) + 'm'
        
    if autoscale:
        fig_text += ' autsc'
        
    fig = plt.figure(fig_text)
    plt.clf()
    fig.set_size_inches((4.5*figsize, 3.25*figsize), forward=True)
        
    power = W.power()
    spec = W.spectrum()
    wigner = W.wig
    wigner_lim = np.amax(abs(W.wig))
    inst_freq = W.inst_freq()
    group_delay = W.group_delay()
    
    if x_units=='fs':
        power_scale = W.s / speed_of_light * 1e15
        p_label_txt = 't [fs]'
        group_delay = group_delay / speed_of_light * 1e15
    else:
        power_scale = W.s*1e6
        p_label_txt = 's [$\mu$m]'
        group_delay = group_delay*1e6
    
    if y_units in ['ev', 'eV']:
        spec_scale = W.phen
        f_label_txt = '$E_{photon}$ [eV]'
        inst_freq = inst_freq
    else:
        spec_scale = W.freq_lamd
        f_label_txt = '$\lambda& [nm]'
        inst_freq = h_eV_s * speed_of_light * 1e9 / inst_freq
    
    if plot_proj:
        # definitions for the axes
        left, width = 0.18, 0.57
        bottom, height = 0.14, 0.55
        left_h = left + width + 0.02 - 0.02
        bottom_h = bottom + height + 0.02 - 0.02
        

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.15, height]
        
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        axScatter = plt.axes(rect_scatter, sharex=axHistx, sharey=axHisty)
    else:
        axScatter = plt.axes()
    
    if abs_value:
        wigplot = axScatter.pcolormesh(power_scale[::downsample], spec_scale[::downsample], abs(wigner[::downsample,::downsample]))
        if plot_text:
            axScatter.text(0.02, 0.98, r'$W_{{max}}$= {:.2e}'.format(np.amax(wigner)), horizontalalignment='left', verticalalignment='top', transform=axScatter.transAxes, color='w')
    else:
        # cmap='RdBu_r'
        # axScatter.imshow(wigner, cmap=cmap, vmax=wigner_lim, vmin=-wigner_lim)
        wigplot = axScatter.pcolormesh(power_scale[::downsample], spec_scale[::downsample], wigner[::downsample,::downsample], cmap=cmap, vmax=wigner_lim, vmin=-wigner_lim)
        if plot_text:
            axScatter.text(0.02, 0.98, r'$W_{{max}}$= {:.2e}'.format(np.amax(wigner)), horizontalalignment='left', verticalalignment='top', transform=axScatter.transAxes)#fontsize=12,
            
    if plot_moments:
        axScatter.plot(power_scale[::downsample], inst_freq[::downsample], "-k")
        axScatter.plot(group_delay[::downsample], spec_scale[::downsample], "-g")

        
    if autoscale == 1:
        autoscale = 1e-2
    
    if autoscale != None:
        max_power = np.amax(power)
        max_spectrum = np.amax(spec)
        idx_p = np.where(power > max_power*autoscale)[0]
        # print(max_spectrum*autoscale)
        # print(spectrum)
        idx_s = np.where(spec > max_spectrum*autoscale)[0]
        
        x_lim = [ power_scale[idx_p[0]], power_scale[idx_p[-1]] ]
        y_lim = [ spec_scale[idx_s[0]], spec_scale[idx_s[-1]] ]
        
        x_lim = np.array(x_lim)
        y_lim = np.array(y_lim)
        x_lim.sort()
        y_lim.sort()
    else:
        x_lim = (np.amin(power_scale), np.amax(power_scale))
        y_lim = (np.amin(spec_scale), np.amax(spec_scale))
        
    if plot_proj:
        axHistx.plot(power_scale,power)
        if plot_text:
            axHistx.text(0.02, 0.95, r'E= {:.2e} J'.format(W.energy()), horizontalalignment='left', verticalalignment='top', transform=axHistx.transAxes)#fontsize=12,
        axHistx.set_ylabel('Power [W]')
        
        if spec.max() <= 0:
            axHisty.plot(spec, spec_scale)
        else:
            axHisty.plot(spec/spec.max(), spec_scale)
        axHisty.set_xlabel('Spectrum [a.u.]')
        
        axScatter.axis('tight')
        axScatter.set_xlabel(p_label_txt)
        axScatter.set_ylabel(f_label_txt)

        axHistx.set_ylim(ymin=0)
        axHisty.set_xlim(xmin=0)

        for tl in axHistx.get_xticklabels():
            tl.set_visible(False)

        for tl in axHisty.get_yticklabels():
            tl.set_visible(False)
        
        axHistx.yaxis.major.locator.set_params(nbins=4)
        axHisty.xaxis.major.locator.set_params(nbins=2)
        
        axHistx.set_xlim(x_lim[0], x_lim[1])
        axHisty.set_ylim(y_lim[0], y_lim[1])
    else:
        axScatter.axis('tight')
        axScatter.set_xlabel(p_label_txt)
        axScatter.set_ylabel(f_label_txt)
    

    
    # axScatter.set_xlim(x_lim[0], x_lim[1])
    # axScatter.set_ylim(y_lim[0], y_lim[1])

    
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if W.z is None:
            fig.savefig(W.filePath + '_wig.' + str(savefig), format=savefig)
        else:
            fig.savefig(W.filePath + '_wig_' + str(W.z) + 'm.' + str(savefig), format=savefig)

    plt.draw()
    
    if showfig == True:
        dir = os.path.dirname(W.filePath)
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)
        

'''
tmp for HXRSS
'''
def read_plot_dump_proj(exp_dir, stage, run_ids, plot_phase=1, showfig=True, savefig=0, debug=1):

    if showfig == 0 and savefig == 0:
        return None

    t_l_int_arr = []
    t_l_pha_arr = []
    f_l_int_arr = []
    for run_id in run_ids:

        array = np.loadtxt(exp_dir + 'run_' + str(run_id) + '/run.' + str(run_id) + '.s' + str(stage) + '.dfl.t.txt', skiprows=1)
        array = np.rollaxis(array, 1)
        t_l_scale, t_l_int_a, t_l_pha_a = array[0], array[1], array[2]

        array = np.loadtxt(exp_dir + 'run_' + str(run_id) + '/run.' + str(run_id) + '.s' + str(stage) + '.dfl.f.txt', skiprows=1)
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
    t_domain = plt.figure(fig_name,figsize=(15,7))
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

    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if debug > 1:
            print('      saving ' + fig_name + '.' + savefig)
        plt.savefig(exp_dir + 'results/' + fig_name + '.' + savefig, format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')

    fig_name = 'stage_' + str(stage) + '__FILT__spectrum'
    f_domain = plt.figure(fig_name,figsize=(15,7))
    ax1 = f_domain.add_subplot(2, 1, 1)
    ax1.plot(h_eV_s * speed_of_light / f_l_scale, f_l_int_arr, '0.5')
    ax1.plot(h_eV_s * speed_of_light / f_l_scale, f_l_int_mean, 'k', linewidth=1.5)
    ax1.grid(True)

    plt.ylabel(r'$P(\lambda)$ [a.u.]')
    ax2 = f_domain.add_subplot(2, 1, 2, sharex=ax1)
    # plt.xlabel(r'$S [\mu m]$')

    ax2.plot(h_eV_s * speed_of_light / f_l_scale, f_l_ftlt_abs**2, 'r')
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

    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if debug > 1:
            print('      saving ' + fig_name + '.' + savefig)
        plt.savefig(exp_dir + 'results/' + fig_name + '.' + savefig, format=savefig)

    if showfig:
        plt.show()
    else:
        plt.close('all')


def plot_dfl_waistscan(sc_res, fig_name=None, showfig=True, savefig=0, debug=1):

    if showfig == False and savefig == False:
        return

    if fig_name is None:
        if sc_res.fileName() is '':
            fig = plt.figure('Waist scan')
        else:
            fig = plt.figure(sc_res.fileName() + ' waist scan')
    else:
        fig = plt.figure(fig_name)

    plt.clf()
    ax_int = fig.add_subplot(1, 1, 1)
    ax_int.plot(sc_res.z_pos, sc_res.phdens_max, 'k', label='max', linewidth=2)
    ax_int.plot(sc_res.z_pos, sc_res.phdens_onaxis, 'grey', label='on-axis')
    ax_int.set_xlabel('z [m]')
    ax_int.set_ylabel('photon density [arb.units]')
    ax_int.legend(loc='lower left')
    ax_size = ax_int.twinx()
    ax_size.plot(sc_res.z_pos, sc_res.fwhm_x * 1e6, 'g--', label='fwhm_x')
    ax_size.plot(sc_res.z_pos, sc_res.fwhm_y * 1e6, 'b--', label='fwhm_y')
    ax_size.plot(sc_res.z_pos, sc_res.std_x * 1e6, 'g:', label='std_x')
    ax_size.plot(sc_res.z_pos, sc_res.std_y * 1e6, 'b:', label='std_y')

    ax_size.set_ylabel('size [um]')
    ax_size.legend(loc='lower right')

    plt.draw()

    if savefig != False:
        if savefig == True:
            savefig = 'png'
        if debug > 0:
            print('      saving *.' + savefig)
        if debug > 1:
            print('        to ' + sc_res.filePath + '_{:.2f}m-{:.2f}m-waistscan.'.format(sc_res.z_pos[0], sc_res.z_pos[-1]) + str(savefig))
        fig.savefig(sc_res.filePath + '_{:.2f}m-{:.2f}m-waistscan.'.format(sc_res.z_pos[0], sc_res.z_pos[-1]) + str(savefig), format=savefig)
    # if debug>0: print('      done in %.2f seconds' % (time.time() - start_time))
    if showfig:
        if debug > 0:
            print('      showing fig')
        plt.show()
    else:
        plt.close('all')


def plot_trf(trf, mode='tr', autoscale=0, showfig=True, savefig=None, fig_name=None):
    '''
    plots TransferFunction() object,
    mode: 
        'tr' - transmission
        'ref' - reflection
    autoscale = scale down to several FWHMma in frequency and several bumps in time
    showfig - display on screen or not
    savefig - path to save png (if any)
    '''
    n_width = 8
    
    l = len(trf.k)
    L = 2*pi / (trf.k[1] - trf.k[0])
    trf_s_td = np.linspace(0, -L, l) * 1e6
    trf_s_fd = trf.ev()
    #trf_s_fd = trf.k
    
    if autoscale:
        trf_s_fd_xlim = np.array([trf.mid_k - n_width*trf.dk, trf.mid_k + n_width*trf.dk])
        trf_s_fd_xlim = h_eV_s*speed_of_light / (2*pi / trf_s_fd_xlim)
        trf_s_fd_xlim = np.sort(trf_s_fd_xlim)
    
    if mode == 'tr':
        trf_fd = deepcopy(trf.tr)
    elif mode == 'ref':
        trf_fd = deepcopy(trf.ref)
    else:
        raise ValueError('mode argument should be "tr" or "ref"')
    
    trf_fd_tmp = trf_fd / (abs(trf_s_td[-1]) / l)
    trf_td = np.fft.ifft(np.fft.fftshift(trf_fd_tmp))
    trf_td = abs(trf_td)**2
    del trf_fd_tmp
    
    if hasattr(trf,'cryst'):
        title = trf.cryst.lattice.element_name + ' ' + str(trf.cryst.ref_idx) + ' ' + mode
    else:
        title = ''
    
    if fig_name is None:
        trf_fig = plt.figure('Filter '+title)
    else:
        trf_fig = plt.figure(fig_name)
    
    trf_fig.set_size_inches((9, 11), forward=True)
    if title != '':
        trf_fig.suptitle(title)
    
    ax_fd_abs = trf_fig.add_subplot(3, 1, 1)
    ax_fd_abs.clear()
    ax_fd_ang = trf_fig.add_subplot(3, 1, 2, sharex=ax_fd_abs)
    ax_fd_ang.clear()
    ax_td = trf_fig.add_subplot(3, 1, 3)
    ax_td.clear()
    
    trf_fig.subplots_adjust(hspace=0)
    trf_fig.subplots_adjust(top=0.95, bottom=0.2, right=0.85, left=0.15)
    
    ax_fd_abs.plot(trf_s_fd, np.abs(trf_fd)**2,'k')
    ax_fd_ang.plot(trf_s_fd, np.angle(trf_fd),'g')
    
    ax_td.semilogy(trf_s_td, trf_td)
    
    ax_fd_abs.set_ylabel(r'|amplitude|$^2$')
    ax_fd_ang.set_ylabel('phase')
    ax_fd_ang.set_xlabel('ph.energy')
    
    ax_td.set_ylabel('impulse responce (power)')
    ax_td.set_xlabel('s [um]')
    
    ax_fd_abs.axis('tight')
    ax_fd_abs.set_ylim([0, 1])
    ax_fd_ang.set_ylim([-np.pi, np.pi])
    if autoscale:
        ax_fd_abs.set_xlim(trf_s_fd_xlim)
    if autoscale:
        ax_td.set_xlim(-n_width*pi / trf.dk * 1e6, 0)
        idx = np.argwhere(trf_s_td>-n_width*pi / trf.dk * 1e6)[-1]
        ax_td.set_ylim(np.amin(trf_td[1:idx]), np.amax(trf_td[1:idx]))
    
    ax_fd_abs.grid(True)
    ax_fd_ang.grid(True)
    ax_td.grid(True)
    
    for label in ax_fd_abs.get_xticklabels():
        label.set_visible(False)
    
    #ax_td.axis('tight')
    
    pos1 = ax_td.get_position()  # get the original position
    pos2 = [pos1.x0 + 0, pos1.y0 - 0.1,  pos1.width / 1.0, pos1.height / 0.9]
    ax_td.set_position(pos2)
    
    plt.draw()
    if savefig != None and savefig.__class__ == str:
        trf_fig.savefig(savefig, format='png')
    #    if savefig == True:
    #        savefig = 'png'
    #    fig.savefig(g.filePath + '_z_' + str(z) + 'm.' + str(savefig), format=savefig)
    
    if showfig:
        plt.show()
    else:
        plt.close('all')
        

def plot_stokes_values(S, fig=None, s_lin=0, norm=0, showfig=True, gw=1, direction='z'):
    
    # if type(S) != StokesParameters:
        # raise ValueError('Not a StokesParameters object')
    if direction == 'z':
       sc = S.sc_z * 1e6
    elif direction == 'x':
       sc = S.sc_x * 1e6
    elif direction == 'y':
       sc = S.sc_y * 1e6
    # sc = S.sc * 1e6
    
    
    if np.size(sc) > 1:
        if fig == None:
            plt.figure('Stokes S')
        elif type(fig) == matplotlib.figure.Figure:
            plt.figure(fig.number)
        else:
            plt.figure(fig)
        plt.clf()
        
        if gw:
            mult = 1e-9
            plt.ylabel('$S_0$ [GW]')
        else:
            mult = 1
            plt.ylabel('$S_0$ [W]')
        plt.xlabel('s [$\mu$m]')
        
        if s_lin:
            # plt.step(sc, np.sqrt(S.s1**2+S.s2**2), linewidth=2, where='mid',color=[0.5,0.5,0.5], linestyle='--')
            plt.step(sc, np.sqrt(S.s1**2+S.s2**2)*mult, linewidth=2, where='mid',color='m', linestyle='--')
 
        plt.step(sc, S.s1*mult, linewidth=2, where='mid',color='g')
        plt.step(sc, S.s2*mult, linewidth=2, where='mid',color='r')
        plt.step(sc, S.s3*mult, linewidth=2, where='mid',color='c')
        plt.step(sc, S.s0*mult, linewidth=2, where='mid',color='b')
        # plt.step(sc, S.s1, linewidth=2, where='mid',color='m')
        # plt.step(sc, S.s2, linewidth=2, where='mid',color='r')
        # plt.step(sc, S.s3, linewidth=2, where='mid',color='c')
        # plt.step(sc, S.s0, linewidth=2, where='mid',color='k')
        
        if s_lin:
            plt.legend(['$\sqrt{S_1^2+S_2^2}$','$S_1$','$S_2$','$S_3$','$S_0$'], loc='lower center', ncol=5, mode="expand", borderaxespad=0.5, frameon=1).get_frame().set_alpha(0.4)
        else:
            plt.legend(['$S_1$','$S_2$','$S_3$','$S_0$'], fontsize=13, ncol=4, loc='upper left', frameon=1).get_frame().set_alpha(0.4)
           # plt.legend(['$S_1$','$S_2$','$S_3$','$S_0$'], loc='lower center', ncol=5, mode="expand", borderaxespad=0.5, frameon=1).get_frame().set_alpha(0.4)
        
        if showfig:
            plt.show()
        else:
            plt.close('all')
        
        
def plot_stokes_angles(S, fig=None, showfig=True, direction='z', scatter=True):
    
   # if type(S) != StokesParameters:
       # raise ValueError('Not a StokesParameters object')
    if direction == 'z':
        sc = S.sc_z * 1e6
    elif direction == 'x':
        sc = S.sc_x * 1e6
    elif direction == 'y':
        sc = S.sc_y * 1e6
    # sc = S.sc * 1e6
    
    if np.size(sc) > 1:
        if fig == None:
            plt.figure('Stokes angles')
        else:
            plt.figure(fig.number)
        plt.clf()

       # plt.step(sc, S.chi(), sc, S.psi(),linewidth=2)
        if scatter:
            psize = S.P_pol()
            psize /= np.amax(psize)
            plt.scatter(sc, S.chi(),psize,linewidth=2,color='g')
            plt.scatter(sc, S.psi(),psize,linewidth=2,color='b')
        else:
            plt.step(sc, S.chi(), linewidth=2, where='mid', color='g')
            plt.step(sc, S.psi(), linewidth=2, where='mid', color='b')
        plt.legend(['$\chi$','$\psi$'])#,loc='best')
        plt.xlabel('s [$\mu$m]')
        plt.ylabel('[rad]')
        plt.ylim([-np.pi/2,np.pi/2])
        plt.xlim([np.amin(sc),np.amax(sc)])
        
        if showfig:
            plt.show()
        else:
            plt.close('all')

def plot_stokes_3d(stk_params, x_plane='max_slice', y_plane='max_slice', z_plane='max_slice', interpolation=None,
                   cmap_lin='brightwheel', cmap_circ='seismic', figsize=4, fig_name='Visualization Stokes parameters',
                   cbars=True, savefig=False, showfig=True, text_present=True, debug=1, **kwargs):
    '''
    Plot 6 images with normalized Stokes parameters on them

    :param stk_params: 3d ocelot.optics.wave.StokesParameters() type object
    :param x_plane: this variable responds on which value on x-axis the 3d stk_params will intersect.
                    It can take 3 different recognition:
                    'max_slice': the intersection of 3d stk_params will contain the max value s0 in stk_params
                    'proj': at the third subplot will be shown the projection of 3d stk_params in x direction
                    <number> in [m]: the position of intersection on x-axis
    :param y_plane: this variable responds on which value on y-axis the 3d stk_params will intersect.
                    It can take 3 different recognition:
                    'max_slice': the intersection of 3d stk_params will contain the max value s0 in stk_params
                    'proj': at the third subplot will be shown the projection of 3d stk_params in y direction
                    <number> in [m]: the position of intersection on y-axis
    :param z_plane: this variable responds on which value on z-axis the 3d stk_params will intersect.
                    It can take 3 different recognition:
                    'max_slice': the intersection of 3d stk_params will contain the max value s0 in stk_params
                    'proj': at the third subplot will be shown the projection of 3d stk_params in z direction
                    <number> in [m]: the position of intersection on z-axis
    :param interpolation: str type variable wich responds for interpolation before plotting linear polarized part
    :param cmap_lin: numpy array with shape (nwidth, nheight, 4) that contains the 4 rgba values in hue (width)
                    and lightness (height).
                    Can be obtained by a call to get_cmap2d(name).
                    or:
                    name where name is one of the following strings:
                    'brightwheel', 'darkwheel', 'hardwheel', 'newwheel',
                    'smoothwheel', 'wheel'
    :param cmap_circ:--------------------
    :param figsize: size of the figure
    :param fig_name: name of the figure
    :param cbars: bool type variable which responds for showing of colorbars
    :param savefig: bool type variable which responds for saving of the figure
    :param showfig: bool type variable which responds for showing of the figure
    :param text_present: bool type variable which responds for showing text on subplots
    :param debug:
    :param kwargs:
    '''

    if showfig == False and savefig == False:
        return
    _logger.info('plotting stokes parameters')
    start_time = time.time()

    # Plotting colorbars
    if cbars:
        _logger.info('plotting colorbars for stokes parameters')
        package_dir = os.path.dirname(__file__)
        cbar_lin_path = os.path.join(package_dir, 'colorbars', 'colorbar_linear_polarization.npy')
        _logger.info('loading ' + cbar_lin_path)
        cbar_lin = np.load(cbar_lin_path)

        cbfig = plt.figure('cbars_stokes_parameters_visualization')
        cbfig.set_size_inches((2 * figsize, 2 * figsize), forward=True)
        cbax1 = cbfig.add_axes([0.15, 0.15, 0.7, 0.7])
        cbar_lin_im = imshow2d(cbar_lin, ax=cbax1, origin='lower', cmap2d=cmap_lin ,interpolation=interpolation,
                               extent=[-1, 1, -1, 1], aspect='auto')
        cbax1.set_title('Colorbar for linear polarization')
        polar_cbax = cbfig.add_axes([0.2, 0.2, 0.6, 0.6], polar=True, frameon=False)
        polar_cbax.set_yticklabels([])
        polar_cbax.tick_params(labelsize=12)
        # polar_cbax.grid(True)


    # Plotting data
    fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches((5 * figsize, 3 * figsize), forward=True)

    z, y, x = stk_params.s0.shape
    ax1 = fig.add_subplot(2, 3, 1)
    linear_plt = plot_stokes_sbfg_lin(ax1, stk_params, slice=z_plane, plane='z', cmap2d=cmap_lin,
                                      plot_title=None, x_label='x', y_label='y', text_present=text_present,
                                      interpolation=interpolation, result=1, **kwargs)

    ax2 = fig.add_subplot(2, 3, 2)
    plot_stokes_sbfg_lin(ax2, stk_params, slice=x_plane, plane='x', cmap2d=cmap_lin, plot_title='Linear polarization',
                         x_label='z', y_label='y', text_present=text_present, interpolation=interpolation, **kwargs)

    ax3 = fig.add_subplot(2, 3, 3)
    plot_stokes_sbfg_lin(ax3, stk_params, slice=y_plane, plane='y', cmap2d=cmap_lin, plot_title=None,
                         x_label='z', y_label='x', text_present=text_present, interpolation=interpolation, **kwargs)

    ax4 = fig.add_subplot(2, 3, 4, sharex=ax1, sharey=ax1)
    circular_plt = plot_stokes_sbfg_circ(ax4, stk_params, slice=z_plane, plane='z', cmap=cmap_circ, plot_title=None,
                                         x_label='x', y_label='y', text_present=text_present, result=1,
                                         interpolation=interpolation, **kwargs)

    if cbars:
        cbax2 = fig.add_axes([0.93, 0.11, 0.02, 0.323])  # This is the position for the colorbar [x, y, width, height]
        cbar_circ_im = plt.colorbar(circular_plt, cax=cbax2)
        cbax2.set_ylabel('Normalized S3 (S3/S0)')
        cbax2.tick_params(axis='both', which='major', labelsize=10)

    ax5 = fig.add_subplot(2, 3, 5, sharex=ax2, sharey=ax2)
    plot_stokes_sbfg_circ(ax5, stk_params, slice=x_plane, plane='x', cmap=cmap_circ, plot_title='Circular polarization',
                          x_label='z', y_label='y', text_present=text_present, interpolation=interpolation, **kwargs)

    ax6 = fig.add_subplot(2, 3, 6, sharex=ax3, sharey=ax3)
    plot_stokes_sbfg_circ(ax6, stk_params, slice=y_plane, plane='y', cmap=cmap_circ, plot_title=None, x_label='z',
                          y_label='x', text_present=text_present, interpolation=interpolation, **kwargs)

    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    _logger.info(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving figure')
        fig.savefig()

    if showfig:
        _logger.debug(ind_str + 'showing Stokes Parameters')
        plt.show()
    else:
        plt.close('all')

    # if cbar:
    #     _logger.info('plotting colorbars for stokes parameters')
    #     package_dir = os.path.dirname(__file__)
    #     cbar_lin_path = os.path.join(package_dir, 'colorbars', 'colorbar_linear_polarization.npy')
    #     _logger.info('loading' + cbar_lin_path)
    #     cbar_lin = np.load(cbar_lin_path)
    #
    #     cbfig = plt.figure('cbars_stokes_parameters_visualization')
    #     cbfig.set_size_inches((5 * figsize, 3 * figsize), forward=True)
    #     cbax1 = cbfig.add_subplot(111)
    #     cbar_lin_im = imshow2d(cbar_lin, ax=cbax1, origin='lower', cmap2d=cmap_lin, interpolation=interpolation,
    #                            extent=[-1, 1, -1, 1], aspect='equal')
    #     cbax1.set_title('colorbar of the visualization of linear polarization')
    #     # cbax2 = cbfig.add_axes([0.85, 0.1, 0.02, 0.8])  # This is the position for the colorbar [x, y, width, height]
    #     # cbax2.set_label('colorbar of the visualization of circular polarization')
    #     # cbar_circ_im = plt.colorbar(circular_plt, cax=cbax2)
    #     # cbar_circ_im.set_label('Normalized S3 (S3/S0)')
    #     # cbax2.tick_params(axis='both', which='major', labelsize=10)
    #     plt.draw()
    #     plt.show()



def plot_stokes_sbfg_lin(ax, stk_params, slice, plane, cmap2d='brightwheel', plot_title=None, x_label='', y_label='',
                         result=0, text_present=True, interpolation=None, **kwargs):
    '''
    Plot normalized intensity and angle of the linear polarization of the light

    :param ax: matplotlib.pyplot.AxesSubplot on which the data will be plotted
    :param stk_params: 3d ocelot.optics.wave.StokesParameters() type object
    :param plane: the direction in which the projection/intersection of {stk_params} will be done
    :param slice: this variable responds on which value on {plane} direction the 3d stk_params will intersect.
                    It can take 3 different recognition:
                    'max_slice': the intersection of 3d stk_params will contain the max value s0 in stk_params
                    'proj': at the third subplot will be shown the projection of 3d stk_params in {plane} direction
                    <number> in [m]: the position of intersection on {plane} direction
    :param cmap2d: numpy array with shape (nwidth, nheight, 4) that contains the 4 rgba values in hue (width)
                    and lightness (height).
                    Can be obtained by a call to get_cmap2d(name).
                    or:
                    name where name is one of the following strings:
                    'brightwheel', 'darkwheel', 'hardwheel', 'newwheel',
                    'smoothwheel', 'wheel'
    :param plot_title: title of the plot
    :param x_label: label of the x axis
    :param y_label: label of the y axis
    :param result: a bool type variable; if bool == True the function will return linear_plt of AxesImage type
    :param text_present: bool type variable which responds for showing text on subplots
    :param interpolation: str type variable wich responds for interpolation before plotting linear polarized part
    :param kwargs:
    :return:
    '''
    # Getting intersections of stk_params for ploting data
    z_max, y_max, x_max = np.unravel_index(stk_params.s0.argmax(), stk_params.s0.shape)  # getting max element position

    if plane in ['x', 2]:
        swap_axes = True
        extent = [stk_params.sc_z[0]*1e6, stk_params.sc_z[-1]*1e6, stk_params.sc_y[0]*1e6, stk_params.sc_y[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(x_max, plane=plane)[:, :, 0]
            slice_pos = stk_params.sc_x[x_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[:, :, 0]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_x, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[:, :, 0]
            slice_pos = stk_params.sc_x[slice_pos]
    elif plane in ['y', 1]:
        swap_axes = True
        extent = [stk_params.sc_z[0]*1e6, stk_params.sc_z[-1]*1e6, stk_params.sc_x[0]*1e6, stk_params.sc_x[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(y_max, plane=plane)[:, 0, :]
            slice_pos = stk_params.sc_y[y_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[:, 0, :]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_y, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[:, 0, :]
            slice_pos = stk_params.sc_y[slice_pos]
    elif plane in ['z', 0]:
        swap_axes = False
        extent = [stk_params.sc_x[0]*1e6, stk_params.sc_x[-1]*1e6, stk_params.sc_y[0]*1e6, stk_params.sc_y[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(z_max, plane=plane)[0, :, :]
            slice_pos = stk_params.sc_z[z_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[0, :, :]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_z, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[0, :, :]
            slice_pos = stk_params.sc_z[slice_pos]
    else:
        _logger.error(ind_str + 'argument "plane" should be in ["x","y","z",0,1,2]')
        raise ValueError('argument "plane" should be in ["x","y","z",0,1,2]')

    # Normalization
    max_s0 = np.amax(stk_params.s0)
    if swap_axes:
        lin_pol_plane = np.swapaxes((stk_params_plane.P_pol_l() / max_s0), 0, 1)
        psi_plane = np.swapaxes((2 * stk_params_plane.psi() / np.pi), 0, 1)
    else:
        lin_pol_plane = stk_params_plane.P_pol_l() / max_s0
        psi_plane = 2 * stk_params_plane.psi() / np.pi

    m, n = psi_plane.shape
    linear_plt = imshow2d(np.array([psi_plane, lin_pol_plane]), ax=ax, cmap2d=cmap2d, extent=extent,
                          interpolation=interpolation, aspect='auto', lightvmin=0, lightvmax=1, origin='lower', **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title, fontsize=15)
    ax.set_xlabel(x_label + ' [$\mu$m]')
    ax.set_ylabel(y_label + ' [$\mu$m]')
    if text_present:
        dic = {'proj': 'projection', 'max_slice': 'slice at {:.3f} $\mu$m (max int)'}
        ax.text(0.97, 0.97, dic.get(slice, 'slice at {:.3f} $\mu$m').format(slice_pos*1e6), horizontalalignment='right',
                  verticalalignment='top', transform=ax.transAxes, fontsize=10)
    if result:
        return linear_plt


def plot_stokes_sbfg_circ(ax, stk_params, slice, plane, cmap='seismic', plot_title=None, x_label='', y_label='',
                          result=0, text_present=True, interpolation=None, **kwargs):
    '''
    Plot normalized Stokes parameter S3

    :param ax: matplotlib.pyplot.AxesSubplot on which the data will be plotted
    :param stk_params: 3d ocelot.optics.wave.StokesParameters() type object
    :param plane: the direction in which the projection/intersection of {stk_params} will be done
    :param slice: this variable responds on which value on {plane} direction the 3d stk_params will intersect.
                    It can take 3 different recognition:
                    'max_slice': the intersection of 3d stk_params will contain the max value s0 in stk_params
                    'proj': at the third subplot will be shown the projection of 3d stk_params in {plane} direction
                    <number> in [m]: the position of intersection on {plane} direction
    :param cmap: colormap which will be used for plotting data
    :param plot_title: title of the plot
    :param x_label: label of the x axis
    :param y_label: label of the y axis
    :param result: a bool type variable; if bool == True the function will return linear_plt of AxesImage type
    :param text_present: bool type variable which responds for showing text on subplots
    :param interpolation: str type variable wich responds for interpolation before plotting linear polarized part
    :param kwargs:
    :return:
    '''

    # Getting intersections of stk_params for ploting data
    z_max, y_max, x_max = np.unravel_index(stk_params.s0.argmax(), stk_params.s0.shape)  # getting max element position

    if plane in ['x', 2]:
        swap_axes = True
        extent = [stk_params.sc_z[0]*1e6, stk_params.sc_z[-1]*1e6, stk_params.sc_y[0]*1e6, stk_params.sc_y[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(x_max, plane=plane)[:, :, 0]
            slice_pos = stk_params.sc_x[x_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[:, :, 0]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_x, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[:, :, 0]
            slice_pos = stk_params.sc_x[slice_pos]
    elif plane in ['y', 1]:
        swap_axes = True
        extent = [stk_params.sc_z[0]*1e6, stk_params.sc_z[-1]*1e6, stk_params.sc_x[0]*1e6, stk_params.sc_x[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(y_max, plane=plane)[:, 0, :]
            slice_pos = stk_params.sc_y[y_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[:, 0, :]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_y, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[:, 0, :]
            slice_pos = stk_params.sc_y[slice_pos]
    elif plane in ['z', 0]:
        swap_axes = False
        extent = [stk_params.sc_x[0]*1e6, stk_params.sc_x[-1]*1e6, stk_params.sc_y[0]*1e6, stk_params.sc_y[-1]*1e6]
        if slice == 'max_slice':
            stk_params_plane = stk_params.slice_2d_idx(z_max, plane=plane)[0, :, :]
            slice_pos = stk_params.sc_z[z_max]
        elif slice == 'proj':
            stk_params_plane = stk_params.proj(plane=plane, mode='mean')[0, :, :]
            slice_pos = 0
        else:
            slice_pos = find_nearest_idx(stk_params.sc_z, slice)
            stk_params_plane = stk_params.slice_2d_idx(slice_pos, plane=plane)[0, :, :]
            slice_pos = stk_params.sc_z[slice_pos]
    else:
        _logger.error(ind_str + 'argument "plane" should be in ["x","y","z",0,1,2]')
        raise ValueError('argument "plane" should be in ["x","y","z",0,1,2]')

    # Normalization
    max_s0 = np.amax(stk_params.s0)
    if swap_axes:
        s3_plane = np.swapaxes((stk_params_plane.s3 / max_s0), 0, 1)
    else:
        s3_plane = stk_params_plane.s3 / max_s0

    m, n = s3_plane.shape
    circular_plt = ax.imshow(s3_plane, cmap=cmap, vmin=-1, vmax=1, interpolation=interpolation, aspect='auto',
                             extent=extent, origin='lower', **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title, fontsize=15)
    ax.set_xlabel(x_label + ' [$\mu$m]')
    ax.set_ylabel(y_label + ' [$\mu$m]')
    if text_present:
        dic = {'proj': 'projection', 'max_slice': 'slice at {:.3f} $\mu$m (max int)'}
        ax.text(0.97, 0.97, dic.get(slice, 'slice at {:.3f} $\mu$m').format(slice_pos*1e6), horizontalalignment='right',
                  verticalalignment='top', transform=ax.transAxes, fontsize=10)
    if result:
        return circular_plt


'''
    scheduled for removal
'''


def show_output(g, show_field=False, show_slice=0):

    print ('plotting slice', show_slice)

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

    print ('wavelength ', xlamds)

    if show_field:
        #from mpi4py import MPI

        #comm = MPI.COMM_WORLD
        #slices = readRadiationFile_mpi(comm=comm, fileName=file+'.dfl', npoints=npoints)
        slices = readRadiationFile(fileName=g.path + '.dfl', npoints=npoints)
        print ('slices:', slices.shape)

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
            P[i] = abs(s * s.conjugate()) * (dgrid**2 / npoints)**2

        t = 1.0e+15 * float(g('zsep')) * float(g('xlamds')) * np.arange(0, len(P)) / c
        plt.plot(t, P)
        plt.title('Pulse/axis')

        fig.add_subplot(133)
        spec = np.abs(np.fft.fft(slices[:, int(npoints / 2), int(npoints / 2)]))**2
        freq_ev = h * fftfreq(len(spec), d=zsep * xlamds / c)
        plt.plot(freq_ev, spec)
        plt.title('Spectrum/axis')


def show_plots(displays, fig):
    '''
    putting arbitrarily many plots on single figure
    '''
    n1 = (len(displays) - 1) / 2 + 1
    n2 = (len(displays) - 1) / n1 + 1
    # print n1, n2
    fmt = str(n1) + str(n2)
    print (fmt)

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


def round_sig(x, sig=2):
    from math import log10, floor
    return round(x, sig - int(floor(log10(x))) - 1)


def gauss_fit(X, Y):
    import numpy as np
    import scipy.optimize as opt

    def gauss(x, p):  # p[0]==mean, p[1]==stdev p[2]==peak
        return p[2] / (p[1] * np.sqrt(2 * np.pi)) * np.exp(-(x - p[0])**2 / (2 * p[1]**2))

    p0 = [0, np.max(X) / 2, np.max(Y)]
    errfunc = lambda p, x, y: gauss(x, p) - y
    p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))
    fit_mu, fit_stdev, ampl = p1
    Y1 = gauss(X, p1)
    RMS = fit_stdev
    return (Y1, RMS)




    # ax_size_l = ax_size_t.twinx() #longitudinal size
    # ax_size_l.plot(g.z, rad_longit_size*2, color='indigo', linestyle='dashed',linewidth=1.5)
    # ax_size_l.set_ylabel('longitudinal [$\mu$m]')
