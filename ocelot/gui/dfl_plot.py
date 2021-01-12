"""
user interface for viewing radiation field
"""

import sys
import os
import csv
import time
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors  # for wigner log scale
import numpy as np
import logging

from ocelot.gui.settings_plot import *

from ocelot.adaptors.genesis import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and
from ocelot.common.math_op import *  # import of mathematical functions like gauss_fit
from ocelot.utils.xfel_utils import *
from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.optics.wave import *

from ocelot.gui.colormaps2d.colormap2d import *

# in order to run decorators properly
import functools

_logger = logging.getLogger(__name__)
__author__ = "Svitozar Serkez, Andrei Trebushinin, Mykola Veremchuk"


@if_plottable
def plot_dfl_all(dfl, **kwargs):
    """
    plots given RadiationField() object in 4 domain combinations
    """
    plot_dfl(dfl, **kwargs)
    dfl.fft_z()
    plot_dfl(dfl, **kwargs)
    dfl.fft_xy()
    plot_dfl(dfl, **kwargs)
    dfl.fft_z()
    plot_dfl(dfl, **kwargs)
    dfl.fft_xy()


@if_plottable
def plot_dfl(dfl, domains=None, z_lim=[], xy_lim=[], figsize=4, cmap=def_cmap, legend=True, phase=False, fig_name=None,
             auto_zoom=False, column_3d=True, savefig=False, showfig=True, return_proj=False, line_off_xy=True,
             slice_xy=False, log_scale=0, cmap_cutoff=0, vartype_dfl=None, **kwargs):
    """
    Plots dfl radiation object in 3d using matplotlib.

    :param dfl: RadiationField() object
    :param domains: longitudinal domain + transverse domain ('t' or 'f' + 's' or 'k') (example: 'tk' - time/inversespace domain)
    :param z_lim: sets the boundaries to CUT the dfl object in z to ranges of e.g. [2,5] um or nm depending on freq_domain=False of True
    :param xy_lim: sets the boundaries to SCALE the dfl object in x and y to ranges of e.g. [2,5] um or urad depending on far_field=False of True
    :param figsize: rescales the size of the figure
    :param cmap: color map which will be used for plotting (http://matplotlib.org/users/colormaps.html)
    :param legend: not used yet
    :param phase: bool type variable, can replace Z projection or spectrum with phase front distribution z dimensions correspondingly
    :param fig_name: the desired name of the output figure, would be used as suffix to the image filename if savefig==True
    :param auto_zoom: bool type variable, automatically scales xyz the images to the (1%?) of the intensity limits
    :param column_3d: bool type variable, plots top and side views of the radiation distribution
    :param savefig: bool type variable, allow to save figure to image (savefig='png' (default) or savefig='eps', etc...)
    :param showfig: bool type variable, allow to display figure (slower)
    :param return_proj: bool type variable, returns [xy_proj,yz_proj,xz_proj,x,y,z] array.
    :param line_off_xy: bool type variable, if True, the transverse size of radiation are calculated at x=0 and y=0 position, otherwise marginal distribution are used
    :param slice_xy: bool type variable, if True, slices will be plotted; if False, projections will be plotted
    :param log_scale: bool type variable, if True, log scale will be used for potting
    :param cmap_cutoff: 0 <= cmap_cutoff <= 1; all pixels that have intensity lower than cmap_cutoff will be seted to white color
    :param vartype_dfl: the data type to store dfl in memory [either complex128 (two 64-bit floats) or complex64 (two 32-bit floats)], may save memory
    :param kwargs: 
    :return:
    """
    import matplotlib.colors as colors

    if showfig == False and savefig == False:
        return

    filePath = dfl.filePath

    text_present = 1
    _logger.info('plotting radiation field (dfl)')
    start_time = time.time()

    if dfl.Nx() == 1 or dfl.Ny() == 1:
        _logger.warning(ind_str + 'plot_dfl() works only with RadiationFields, with dfl.Nx(), dfl.Ny() > 1')
    # print('dfl type is ',type(dfl))
    # if isinstance(dfl, RadiationField):
    # # if dfl.__class__ != RadiationField:
    #     raise ValueError('wrong radiation object: should be RadiationField')

    if vartype_dfl is not None:
        dfl_copy = RadiationField()
        dfl_copy.copy_param(dfl, version=2)
        dfl_copy.fld = dfl.fld.astype(vartype_dfl)
    else:
        dfl_copy = deepcopy(dfl)

    if domains is None:
        domains = dfl_copy.domains()
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
    #     suffix = ''
    # else:
    #     suffix = '_'+fig_name

    if dfl_copy.Nz() != 1:
        # Make sure it is time-dependent
        ncar_z = dfl_copy.Nz()
        leng_z = dfl_copy.Lz()
        z = np.linspace(0, leng_z, ncar_z)
    else:
        column_3d = False
        phase = True
        freq_domain = False
        z_lim = []
    xlamds = dfl_copy.xlamds

    # number of mesh points
    ncar_x = dfl_copy.Nx()
    leng_x = dfl_copy.Lx()  # transverse size of mesh [m]
    ncar_y = dfl_copy.Ny()
    leng_y = dfl_copy.Ly()
    E_pulse = dfl_copy.E()

    if dfl_copy.Nz() != 1:
        if freq_domain:
            if dfl_copy.domain_z == 't':
                dfl_copy.fft_z()

            # z = dfl_copy.scale_z() * 1e9
            # dfl_copy.fld = dfl_copy.fld[::-1, :, :]
            # z = z[::-1]
            # unit_z = r'nm'
            # z_label = r'$\lambda$ [' + unit_z + ']'

            z = h_eV_s * speed_of_light / dfl_copy.scale_z()
            unit_z = r'eV'
            z_label = r'$E_{{ph}}$ [{}]'.format(unit_z)

            z_labelv = r'[arb. units]'
            z_title = 'Spectrum'
            z_color = 'red'
            suffix += '_fd'
        else:
            if dfl_copy.domain_z == 'f':
                dfl_copy.fft_z()
            z = dfl_copy.scale_z() * 1e6

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
        _logger.debug(ind_str + 'setting z-axis limits to ' + str(np.amin(z)) + ':' + str(z_lim[0]) + '-' + str(
            z_lim[1]) + ':' + str(np.amax(z)))  # tmp
        z_lim_1 = np.where(z <= z_lim[0])[0][-1]
        z_lim_2 = np.where(z >= z_lim[1])[0][0]

        if z_lim_1 == z_lim_2 and z_lim_1 == 0:
            z_lim_2 = z_lim_1 + 1
        elif z_lim_1 == z_lim_2 and z_lim_1 != 0:
            z_lim_1 = z_lim_2 - 1
        dfl_copy.fld = dfl_copy.fld[z_lim_1:z_lim_2, :, :]
        z = z[z_lim_1:z_lim_2]
        ncar_z = dfl_copy.shape[0]
        suffix += '_zoom_%.2f-%.2f' % (np.amin(z), np.amax(z))

    if far_field:
        if dfl_copy.domain_xy == 's':
            dfl_copy.fft_xy()
        x = dfl_copy.scale_x() * 1e6
        y = dfl_copy.scale_y() * 1e6

        unit_xy = r'$\mu$rad'
        x_label = r'$\theta_x$ [' + unit_xy + ']'
        y_label = r'$\theta_y$ [' + unit_xy + ']'
        suffix += '_ff'
        x_title = 'X divergence'
        y_title = 'Y divergence'
        xy_title = 'Far field intensity'
        x_y_color = 'green'
    else:
        if dfl_copy.domain_xy == 'k':
            dfl_copy.fft_xy()
        x = dfl_copy.scale_x() * 1e6
        y = dfl_copy.scale_y() * 1e6

        unit_xy = r'$\mu$m'
        x_label = 'x [' + unit_xy + ']'
        y_label = 'y [' + unit_xy + ']'
        x_title = 'X projection'
        y_title = 'Y projection'
        xy_title = 'Intensity'
        x_y_color = 'blue'

    dfl_copy.fld = dfl_copy.fld.astype(np.complex64)
    xy_proj = dfl_copy.int_xy()
    xy_proj_ph = np.angle(np.sum(dfl_copy.fld, axis=0))  # tmp  # tmp
    if slice_xy:
        yz_proj = dfl_copy.intensity()[:, :, dfl_copy.Nx() // 2]
        xz_proj = dfl_copy.intensity()[:, dfl_copy.Ny() // 2, :]
        xz_title = 'Top slice y=0'
        yz_title = 'Side slice x=0'
        z_proj = dfl_copy.intensity()[:, dfl_copy.Ny() // 2, dfl_copy.Nx() // 2]
        z_title += ' (on-axis)'
    else:
        yz_proj = dfl_copy.int_zy()
        xz_proj = dfl_copy.int_zx()
        xz_title = 'Top projection'
        yz_title = 'Side projection'
        z_proj = dfl_copy.int_z()

    dx = abs(x[1] - x[0])
    dy = abs(y[1] - y[0])

    if log_scale:
        suffix += '_log'

    if fig_name is None:
        if dfl_copy.fileName() == '':
            fig = plt.figure('Radiation distribution' + suffix)
        else:
            fig = plt.figure('Radiation distribution' + suffix + ' ' + dfl_copy.fileName())
    else:
        fig = plt.figure(fig_name + suffix)
    del dfl_copy

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

    if cmap_cutoff not in [None, False, 0]:
        cmap = matplotlib.cm.get_cmap(cmap)
        cmap.set_under("w")
        xy_proj[xy_proj < xy_proj.max() * cmap_cutoff] = -1e-10
        yz_proj[yz_proj < yz_proj.max() * cmap_cutoff] = -1e-10
        xz_proj[xz_proj < xz_proj.max() * cmap_cutoff] = -1e-10

    if log_scale:
        xy_proj[xy_proj <= 0] = None
        yz_proj[yz_proj <= 0] = None
        xz_proj[xz_proj <= 0] = None
        z_proj[z_proj <= 0] = None

    ax_int = fig.add_subplot(2, 2 + column_3d, 1)
    if log_scale:
        intplt = ax_int.pcolormesh(x, y, xy_proj, norm=colors.LogNorm(vmin=np.nanmin(xy_proj), vmax=np.nanmax(xy_proj)),
                                   cmap=cmap)
    else:
        intplt = ax_int.pcolormesh(x, y, xy_proj, cmap=cmap, vmin=0)
    ax_int.set_title(xy_title, fontsize=15)
    ax_int.set_xlabel(r'' + x_label)
    ax_int.set_ylabel(y_label)
    if np.size(z) > 1 and kwargs.get('showtext', True):
        ax_int.text(0.01, 0.01, r'$E_{p}$=%.2e J' % (E_pulse), horizontalalignment='left', verticalalignment='bottom',
                    fontsize=12, color='white', transform=ax_int.transAxes)

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
        ax_proj_x.set_ylim(ymin=np.amin(x_line), ymax=1)
    else:
        ax_proj_x.plot(x, x_line, linewidth=2, color=x_y_color)
        ax_proj_x.plot(x, x_line_f, color='grey')
        ax_proj_x.set_ylim(ymin=0, ymax=1)

    if kwargs.get('showtext', True):
        try:
            ax_proj_x.text(0.95, 0.95,
                           'fwhm={:.3g} '.format(fwhm_x) + r' [{:}]'.format(unit_xy) + '\nrms={:.3g}'.format(
                               rms_x) + r' [{:}]'.format(unit_xy), horizontalalignment='right', verticalalignment='top',
                           transform=ax_proj_x.transAxes, fontsize=12)
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
        ax_proj_y.set_xlim(xmin=np.nanmin(y_line), xmax=1)
    else:
        ax_proj_y.plot(y_line, y, linewidth=2, color=x_y_color)
        ax_proj_y.plot(y_line_f, y, color='grey')
        ax_proj_y.set_xlim(xmin=0, xmax=1)

    if kwargs.get('showtext', True):
        try:
            ax_proj_y.text(0.95, 0.95,
                           'fwhm={:.3g} '.format(fwhm_y) + r' [{:}]'.format(unit_xy) + '\nrms={:.3g}'.format(
                               rms_y) + r' [{:}]'.format(unit_xy), horizontalalignment='right', verticalalignment='top',
                           transform=ax_proj_y.transAxes, fontsize=12)
        except:
            pass

    ax_proj_y.set_ylabel(y_label)

    # if log_scale:
    #     ax_proj_x.set_yscale('log')
    #     ax_proj_y.set_xscale('log')
    #     if not phase:
    #         ax_z.set_yscale('log')

    if column_3d:

        if log_scale:
            cut_off = 1e-6
            yz_proj[yz_proj < np.nanmax(yz_proj) * cut_off] = 0
            xz_proj[xz_proj < np.nanmax(xz_proj) * cut_off] = 0
            # cut-off = np.amin([yz_proj[yz_proj!=0].min(), xz_proj[xz_proj!=0].min()]) / 10
            # yz_proj += minmin
            # xz_proj += minmin
            min_xz_proj = np.nanmin(xz_proj[xz_proj != 0])
            min_yz_proj = np.nanmin(yz_proj[yz_proj != 0])

        # if np.amin(xz_proj) == 0:
        #     min_xz_proj = 0
        # else:
        #     min_xz_proj=xz_proj[xz_proj!=0].min()
        # if np.amin(yz_proj) == 0:
        #     min_yz_proj = 0
        # else:
        #     min_yz_proj=yz_proj[yz_proj!=0].min()

        if phase == True:
            ax_proj_xz = fig.add_subplot(2, 2 + column_3d, 6)
        else:
            ax_proj_xz = fig.add_subplot(2, 2 + column_3d, 6, sharex=ax_z)
        if log_scale:
            ax_proj_xz.pcolormesh(z, x, np.swapaxes(xz_proj, 1, 0),
                                  norm=colors.LogNorm(vmin=min_xz_proj, vmax=np.nanmax(xz_proj)), cmap=cmap)
        else:
            ax_proj_xz.pcolormesh(z, x, np.swapaxes(xz_proj, 1, 0), cmap=cmap, vmin=0)
        ax_proj_xz.set_title(xz_title, fontsize=15)
        ax_proj_xz.set_xlabel(z_label)
        ax_proj_xz.set_ylabel(x_label)

        ax_proj_yz = fig.add_subplot(2, 2 + column_3d, 3, sharey=ax_int, sharex=ax_proj_xz)
        if log_scale:
            ax_proj_yz.pcolormesh(z, y, np.swapaxes(yz_proj, 1, 0),
                                  norm=colors.LogNorm(vmin=min_yz_proj, vmax=np.nanmax(yz_proj)), cmap=cmap)
        else:
            ax_proj_yz.pcolormesh(z, y, np.swapaxes(yz_proj, 1, 0), cmap=cmap, vmin=0)
        ax_proj_yz.set_title(yz_title, fontsize=15)
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
        #     ax_z.set_xlim(zlim_calc)
        #     _logger.debug(ind_str + 'scaling xy to {:}'.format(size_xy))
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
        _logger.debug(ind_str + 'saving *{:}.{:}'.format(suffix, savefig))
        fig.savefig(filePath + suffix + '.' + str(savefig), format=savefig)
    _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))

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


@if_plottable
def plot_wigner(wig_or_out, z=np.inf, x_units='um', y_units='ev', x_lim=(None, None), y_lim=(None, None), downsample=1,
                autoscale=None, figsize=3, cmap='seismic', fig_name=None, savefig=False, showfig=True,
                plot_proj=1, plot_text=1, plot_moments=0, plot_cbar=0, log_scale=0, **kwargs):
    """
    Plots wigner distribution (WD) with marginals

    :param wig_or_out: may be WignerDistribution() or GenesisOutput() object
    :param z: (if isinstance(wig_or_out, GenesisOutput)) location at which WD will be calculated
    :param x_units: [m or fs] units to display power scale
    :param y_units: [nm or eV] units to display spectrum scale
    :param x_lim: scaling limits for x in given units, (min,max) or [min,max], e.g: (None,6)
    :param x_lim: scaling limits for y in given units, (min,max) or [min,max], e.g: (None,6)
    :param downsample: speeds up plotting by displaying only 1/downsample**2 points
    :param autoscale: find x_lim and x_lim values automatically. Only (values > max_value * autoscale) will be displayed
    :param figsize: rescales the size of the figure
    :param cmap: colormar (http://matplotlib.org/users/colormaps.html)
    :param fig_name: the desired name of the output figure, would be used as suffix to the image filename if savefig==True
    :param savefig: bool type variable, allow to save figure to image (savefig='png' (default) or savefig='eps', etc...)
    :param showfig: bool type variable, allow to display figure (slower)
    :param plot_proj: plot marginal distributions
    :param plot_text: show text
    :param plot_moments: plot moments as lines on top of Wigner distribution
    :param plot_cbar: plots colorbar
    :param log_scale: plots wigner distribution in logarithmic scale
    :param kwargs:
    :return: None
    """
    if showfig == False and savefig == False:
        return

    _logger.info('plotting Wigner distribution')

    if not hasattr(wig_or_out, 'wig') and hasattr(wig_or_out, 'calc_radsize'):
        W = wigner_out(wig_or_out, z)
    elif hasattr(wig_or_out, 'wig'):
        W = wig_or_out
    else:
        raise ValueError('Unknown object for Wigner plot')

    if fig_name is None:
        if W.fileName() == '':
            fig_text = 'Wigner distribution'
        else:
            fig_text = 'Wigner distribution ' + W.fileName()
    else:
        fig_text = fig_name
    if W.z != None:
        fig_text += ' ' + str(W.z) + 'm'

    if autoscale:
        fig_text += ' autsc'

    fig = plt.figure(fig_text)
    plt.clf()
    fig.set_size_inches((4.5 * figsize, 3.25 * figsize), forward=True)

    power = W.power()
    spec = W.spectrum()
    wigner = W.wig
    wigner_lim = np.amax(abs(W.wig))
    if plot_moments:
        inst_freq = W.inst_freq()
        group_delay = W.group_delay()

    if x_units == 'fs':
        power_scale = -W.s / speed_of_light * 1e15
        p_label_txt = 't [fs]'
        if plot_moments:
            group_delay = group_delay / speed_of_light * 1e15
    else:
        power_scale = W.s * 1e6
        p_label_txt = 's [$\mu$m]'
        if plot_moments:
            group_delay = group_delay * 1e6

    if y_units in ['ev', 'eV']:
        spec_scale = W.phen
        f_label_txt = '$E_{photon}$ [eV]'
        if plot_moments:
            inst_freq = inst_freq
    else:
        spec_scale = W.freq_lamd
        f_label_txt = '$\lambda$ [nm]'
        if plot_moments:
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

    # cmap='RdBu_r'
    # axScatter.imshow(wigner, cmap=cmap, vmax=wigner_lim, vmin=-wigner_lim)

    if log_scale != 0:
        if log_scale == 1:
            log_scale = 0.01
        wigplot = axScatter.pcolormesh(power_scale[::downsample], spec_scale[::downsample],
                                       wigner[::downsample, ::downsample], cmap=cmap,
                                       norm=colors.SymLogNorm(linthresh=wigner_lim * log_scale, linscale=2,
                                                              vmin=-wigner_lim, vmax=wigner_lim),
                                       vmax=wigner_lim, vmin=-wigner_lim)
    else:
        wigplot = axScatter.pcolormesh(power_scale[::downsample], spec_scale[::downsample],
                                       wigner[::downsample, ::downsample], cmap=cmap, vmax=wigner_lim, vmin=-wigner_lim)

    if plot_cbar:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cbaxes = inset_axes(axScatter, width="50%", height="3%", loc=1)
        fig.colorbar(wigplot, cax=cbaxes, orientation='horizontal')

    if plot_text:
        if hasattr(wig_or_out, 'is_spectrogram'):
            if wig_or_out.is_spectrogram:
                axScatter.text(0.98, 0.98, 'Spectrogram', horizontalalignment='right', verticalalignment='top',
                               transform=axScatter.transAxes, color='red')
        axScatter.text(0.02, 0.98, r'$W_{{max}}$= {:.2e}'.format(np.amax(wigner)), horizontalalignment='left',
                       verticalalignment='top', transform=axScatter.transAxes)  # fontsize=12,
        if hasattr(W, 'on_axis'):
            if W.on_axis == True:
                axScatter.text(0.5, 0.98, r"(on axis)", fontsize=10, horizontalalignment='center',
                               verticalalignment='top', transform=axScatter.transAxes)
            else:
                axScatter.text(0.5, 0.98, r"(assuming full spatial coherence)", fontsize=10,
                               horizontalalignment='center', verticalalignment='top', transform=axScatter.transAxes)

    if plot_moments:
        weight_power = power / np.max(power)
        weight_power[weight_power < np.nanmax(weight_power) / 1e2] = 0
        idx_power_fine = np.where(weight_power > np.nanmax(weight_power) / 1e2)
        weight_spec = spec / np.max(spec)
        weight_spec[weight_spec < np.nanmax(weight_spec) / 1e2] = 0
        idx_spec_fine = np.where(weight_spec > np.nanmax(weight_spec) / 1e2)

        plt.scatter(power_scale[idx_power_fine], inst_freq[idx_power_fine], s=weight_power[idx_power_fine], c='black',
                    linewidths=2)
        plt.scatter(group_delay[idx_spec_fine], spec_scale[idx_spec_fine], s=weight_spec[idx_spec_fine], c='green',
                    linewidths=2)
        # axScatter.plot(power_scale[::downsample], inst_freq[::downsample], "-k")
        # axScatter.plot(group_delay[::downsample], spec_scale[::downsample], "-g")

    if autoscale == 1:
        autoscale = 1e-2

    if autoscale not in [0, None]:
        max_power = np.amax(power)
        max_spectrum = np.amax(spec)
        idx_p = np.where(power > max_power * autoscale)[0]
        idx_s = np.where(spec > max_spectrum * autoscale)[0]

        x_lim_appl = [power_scale[idx_p[0]], power_scale[idx_p[-1]]]
        x_lim_appl = np.array(x_lim_appl)
        x_lim_appl.sort()

        y_lim_appl = [spec_scale[idx_s[0]], spec_scale[idx_s[-1]]]
        y_lim_appl = np.array(y_lim_appl)
        y_lim_appl.sort()

    else:
        x_lim_appl = (np.amin(power_scale), np.amax(power_scale))
        y_lim_appl = (np.amin(spec_scale), np.amax(spec_scale))

    if x_units == 'fs':
        x_lim_appl = np.flipud(x_lim_appl)

    if x_lim[0] is not None:
        x_lim_appl[0] = x_lim[0]
    if x_lim[1] is not None:
        x_lim_appl[1] = x_lim[1]
    if y_lim[0] is not None:
        y_lim_appl[0] = y_lim[0]
    if y_lim[1] is not None:
        y_lim_appl[1] = y_lim[1]

    if plot_proj:
        axHistx.plot(power_scale, power)
        if plot_text:
            axHistx.text(0.02, 0.95, r'E= {:.2e} J'.format(W.energy()), horizontalalignment='left',
                         verticalalignment='top', transform=axHistx.transAxes)  # fontsize=12,
        axHistx.set_ylabel('Power [W]')

        if spec.max() <= 0:
            axHisty.plot(spec, spec_scale)
        else:
            axHisty.plot(spec / spec.max(), spec_scale)

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

        axHistx.set_xlim(x_lim_appl[0], x_lim_appl[1])
        axHisty.set_ylim(y_lim_appl[0], y_lim_appl[1])

        if log_scale != 0:
            axHistx.set_ylim(np.nanmin(power), np.nanmax(power))
            axHisty.set_xlim(np.nanmin(spec), np.nanmax(spec))
            axHisty.set_xscale('log')
            axHistx.set_yscale('log')

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
            save_path = W.filePath + '_wig.' + str(savefig)
            # fig.savefig(W.filePath + '_wig.' + str(savefig), format=savefig)
        else:
            save_path = W.filePath + '_wig_' + str(W.z) + 'm.' + str(savefig)
            # fig.savefig(W.filePath + '_wig_' + str(W.z) + 'm.' + str(savefig), format=savefig)
        _logger.debug(ind_str + 'saving to {}'.format(save_path))
        fig.savefig(save_path, format=savefig)

    plt.draw()

    if showfig == True:
        dir = os.path.dirname(W.filePath)
        rcParams["savefig.directory"] = dir
        plt.show()
    else:
        # plt.close('all')
        plt.close(fig)


@if_plottable
def plot_dfl_waistscan(sc_res, fig_name=None, figsize=4, showfig=True, savefig=False):
    _logger.info('plot dfl waist scan')
    if showfig == False and savefig == False:
        return

    if fig_name is None:
        if sc_res.fileName() == '':
            fig = plt.figure('Waist scan')
        else:
            fig = plt.figure(sc_res.fileName() + ' waist scan')
    else:
        fig = plt.figure(fig_name)

    plt.clf()
    fig.set_size_inches((3 * figsize, 2 * figsize), forward=True)
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
        _logger.debug(ind_str + 'saving *.' + savefig)
        _logger.debug(ind_str + 'to ' + sc_res.filePath + '_{:.2f}m-{:.2f}m-waistscan.'.format(sc_res.z_pos[0],
                                                                                               sc_res.z_pos[-1]) + str(
            savefig))
        fig.savefig(
            sc_res.filePath + '_{:.2f}m-{:.2f}m-waistscan.'.format(sc_res.z_pos[0], sc_res.z_pos[-1]) + str(savefig),
            format=savefig)
    if showfig:
        _logger.debug(ind_str + 'showing fig')
        plt.show()
    else:
        plt.close('all')


@if_plottable
def plot_trf(trf, mode='tr', autoscale=0, showfig=True, savefig=None, fig_name=None):
    """
    plots TransferFunction() object,
    mode: 
        'tr' - transmission
        'ref' - reflection
    autoscale = scale down to several FWHMma in frequency and several bumps in time
    showfig - display on screen or not
    savefig - path to save png (if any)
    """
    n_width = 8

    l = len(trf.k)
    L = 2 * pi / (trf.k[1] - trf.k[0])
    trf_s_td = np.linspace(0, -L, l) * 1e6
    trf_s_fd = trf.ev()
    # trf_s_fd = trf.k

    if autoscale:
        trf_s_fd_xlim = np.array([trf.mid_k - n_width * trf.dk, trf.mid_k + n_width * trf.dk])
        trf_s_fd_xlim = h_eV_s * speed_of_light / (2 * pi / trf_s_fd_xlim)
        trf_s_fd_xlim = np.sort(trf_s_fd_xlim)

    if mode == 'tr':
        trf_fd = deepcopy(trf.tr)
    elif mode == 'ref':
        trf_fd = deepcopy(trf.ref)
    else:
        raise ValueError('mode argument should be "tr" or "ref"')

    trf_fd_tmp = trf_fd / (abs(trf_s_td[-1]) / l)
    trf_td = np.fft.ifft(np.fft.fftshift(trf_fd_tmp))
    trf_td = abs(trf_td) ** 2
    del trf_fd_tmp

    if hasattr(trf, 'cryst'):
        title = trf.cryst.lattice.element_name + ' ' + str(trf.cryst.ref_idx) + ' ' + mode
    else:
        title = ''

    if fig_name is None:
        trf_fig = plt.figure('Filter ' + title)
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

    ax_fd_abs.plot(trf_s_fd, np.abs(trf_fd) ** 2, 'k')
    ax_fd_ang.plot(trf_s_fd, np.angle(trf_fd), 'g')

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
        ax_td.set_xlim(-n_width * pi / trf.dk * 1e6, 0)
        idx = np.argwhere(trf_s_td > -n_width * pi / trf.dk * 1e6)[-1]
        ax_td.set_ylim(np.amin(trf_td[1:idx]), np.amax(trf_td[1:idx]))

    ax_fd_abs.grid(True)
    ax_fd_ang.grid(True)
    ax_td.grid(True)

    for label in ax_fd_abs.get_xticklabels():
        label.set_visible(False)

    # ax_td.axis('tight')

    pos1 = ax_td.get_position()  # get the original position
    pos2 = [pos1.x0 + 0, pos1.y0 - 0.1, pos1.width / 1.0, pos1.height / 0.9]
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


@if_plottable
def plot_stokes_values(S, fig=None, d_pol=0, norm=0, showfig=True, gw=1, direction='z', plot_func='step', **kwargs):
    # if type(S) != StokesParameters:
    #     raise ValueError('Not a StokesParameters object')
    if direction == 'z':
        sc = S.sc_z * 1e6
        Scp = S[:, 0, 0]  ##### tbd: calculate middle?
    elif direction == 'x':
        sc = S.sc_x * 1e6
        Scp = S[0, 0, :]
    elif direction == 'y':
        sc = S.sc_y * 1e6
        Scp = S[0, :, 0]

    if np.size(sc) <= 1:
        _logger.warning('plot_stokes_values needs more than a single point to plot (np.size(sc) <= 1)')
        return

    if d_pol != 0:
        gw = 0
        norm = 1

    if fig == None:
        plt.figure('Stokes S')
        plt.clf()
    elif type(fig) == matplotlib.figure.Figure:
        plt.figure(fig.number)
    else:
        plt.figure(fig)
    plt.clf()

    if gw:
        mult = 1e-9
        plt.ylabel('$S_0$ [GW]')
    elif norm:
        mult = 1 / np.amax(Scp.s0)
    else:
        mult = 1
        plt.ylabel('$S_0$ [W]')
    plt.xlabel('s [$\mu$m]')

    kwargs = {'linewidth': 2}

    if plot_func == 'step':
        plot_function = plt.step
        kwargs['where'] = 'mid'
    elif plot_func == 'line':
        plot_function = plt.plot
    else:
        raise ValueError

    if d_pol == 'lin':
        # plt.step(sc, np.sqrt(S.s1**2+S.s2**2), linewidth=2, where='mid',color=[0.5,0.5,0.5], linestyle='--')
        plot_function(sc, Scp.deg_pol_l(), linestyle='-', color='#1f77b4', **kwargs)
    elif d_pol == 1:
        plot_function(sc, Scp.deg_pol(), linestyle='-', color='#1f77b4', **kwargs)
    else:
        pass

    plot_function(sc, Scp.s1 * mult, color='g', **kwargs)
    plot_function(sc, Scp.s2 * mult, color='r', **kwargs)
    plot_function(sc, Scp.s3 * mult, color='c', **kwargs)
    plot_function(sc, Scp.s0 * mult, color='b', **kwargs)
    # plt.step(sc, S.s1, linewidth=2, where='mid',color='m')
    # plt.step(sc, S.s2, linewidth=2, where='mid',color='r')
    # plt.step(sc, S.s3, linewidth=2, where='mid',color='c')
    # plt.step(sc, S.s0, linewidth=2, where='mid',color='k')

    if d_pol == 'lin':
        plt.legend(['$D_{lin}$', '$S_1$', '$S_2$', '$S_3$', '$S_0$'], loc='lower center', ncol=5,
                   mode="expand", borderaxespad=0.5, frameon=1).get_frame().set_alpha(0.4)
    elif d_pol == 1:
        plt.legend(['$D_{pol}$', '$S_1$', '$S_2$', '$S_3$', '$S_0$'], loc='lower center', ncol=5,
                   mode="expand", borderaxespad=0.5, frameon=1).get_frame().set_alpha(0.4)
    else:
        plt.legend(['$S_1$', '$S_2$', '$S_3$', '$S_0$'], fontsize=13, ncol=4, loc='upper left',
                   frameon=1).get_frame().set_alpha(0.4)
        # plt.legend(['$S_1$','$S_2$','$S_3$','$S_0$'], loc='lower center', ncol=5, mode="expand", borderaxespad=0.5, frameon=1).get_frame().set_alpha(0.4)
    plt.xlim([np.amin(sc), np.amax(sc)])
    if norm:
        plt.ylim([-1, 1])
    plt.draw()
    if showfig:
        plt.show()
    else:
        plt.close('all')


@if_plottable
def plot_stokes_angles(S, fig=None, showfig=True, direction='z', plot_func='scatter'):
    # if type(S) != StokesParameters:
    #     raise ValueError('Not a StokesParameters object')
    if direction == 'z':
        sc = S.sc_z * 1e6
        Scp = S[:, 0, 0]
    elif direction == 'x':
        sc = S.sc_x * 1e6
        Scp = S[0, 0, :]
    elif direction == 'y':
        sc = S.sc_y * 1e6
        Scp = S[0, :, 0]
    # sc = S.sc * 1e6

    if np.size(sc) <= 1:
        _logger.warning('plot_stokes_angles needs more than a single point to plot (np.size(sc) <= 1)')
        return

    if fig == None:
        plt.figure('Stokes angles')
        plt.clf()
    else:
        plt.figure(fig.number)
    plt.clf()

    kwargs = {'linewidth': 2}
    if plot_func == 'scatter':
        psize = Scp.deg_pol()
        kwargs['s'] = psize
        plot_function = plt.scatter
    elif plot_func == 'step':
        plot_function = plt.step
        kwargs['where'] = 'mid'
    elif plot_func == 'line':
        plot_function = plt.plot
    else:
        raise ValueError

    # plt.step(sc, S.chi(), sc, S.psi(),linewidth=2)

    plot_function(sc, Scp.psi(), color='b', **kwargs)
    # if plot_func == 'scatter':
    # kwargs['s'] = psize
    plot_function(sc, Scp.chi(), color='g', **kwargs)

    # if scatter:
    #     psize = Scp.P_pol()
    #     psize /= np.amax(psize)
    #     plt.scatter(sc, Scp.chi(),psize,linewidth=2,color='g')
    #     plt.scatter(sc, Scp.psi(),psize,linewidth=2,color='b')
    # else:
    #     plt.step(sc, Scp.chi(), linewidth=2, where='mid', color='g')
    #     plt.step(sc, Scp.psi(), linewidth=2, where='mid', color='b')
    plt.legend(['$\chi$', '$\psi$'])  # ,loc='best')
    plt.xlabel('s [$\mu$m]')
    plt.ylabel('[rad]')
    plt.ylim([-np.pi / 2, np.pi / 2])
    plt.xlim([np.amin(sc), np.amax(sc)])
    plt.draw()
    if showfig:
        plt.show()
    else:
        plt.close('all')


@if_plottable
def plot_stokes_3d(stk_params, x_plane='max_slice', y_plane='max_slice', z_plane='max_slice', interpolation=None,
                   cmap_lin='brightwheel', cmap_circ='bwr', figsize=4, fig_name='Visualization Stokes parameters',
                   normalization='s0_max', cbars=True, savefig=False, showfig=True, text_present=True, **kwargs):
    """
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
    :param kwargs:
    """

    if showfig == False and savefig == False:
        return
    _logger.info('plotting stokes parameters')
    start_time = time.time()

    cbax1_dir = kwargs.pop('cbax1_dir', 1)

    ny_plots = 2
    # Plotting data
    fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches((5 * figsize, 3 * figsize), forward=True)

    z, y, x = stk_params.s0.shape
    ax1 = fig.add_subplot(ny_plots, 3, 1)
    linear_plt = plot_stokes_sbfg_lin(ax1, stk_params, slice=z_plane, plane='z', cmap2d=cmap_lin, plot_title=None,
                                      x_label='x', y_label='y', text_present=text_present, interpolation=interpolation,
                                      normalization=normalization, result=1, **kwargs)

    ax2 = fig.add_subplot(ny_plots, 3, 2, sharey=ax1)
    plot_stokes_sbfg_lin(ax2, stk_params, slice=x_plane, plane='x', cmap2d=cmap_lin, plot_title='Linear polarization',
                         x_label='z', y_label='y', text_present=text_present, interpolation=interpolation,
                         normalization=normalization, **kwargs)

    ax3 = fig.add_subplot(ny_plots, 3, 3, sharex=ax2)
    plot_stokes_sbfg_lin(ax3, stk_params, slice=y_plane, plane='y', cmap2d=cmap_lin, plot_title=None, x_label='z',
                         y_label='x', text_present=text_present, interpolation=interpolation,
                         normalization=normalization, **kwargs)

    ax4 = fig.add_subplot(ny_plots, 3, 4, sharex=ax1, sharey=ax1)
    circular_plt = plot_stokes_sbfg_circ(ax4, stk_params, slice=z_plane, plane='z', cmap=cmap_circ, plot_title=None,
                                         x_label='x', y_label='y', text_present=text_present, result=1,
                                         interpolation=interpolation, normalization=normalization, **kwargs)

    ax5 = fig.add_subplot(ny_plots, 3, 5, sharex=ax2, sharey=ax2)
    plot_stokes_sbfg_circ(ax5, stk_params, slice=x_plane, plane='x', cmap=cmap_circ, plot_title='Circular polarization',
                          x_label='z', y_label='y', text_present=text_present, interpolation=interpolation,
                          normalization=normalization, **kwargs)

    ax6 = fig.add_subplot(ny_plots, 3, 6, sharex=ax3, sharey=ax3)
    plot_stokes_sbfg_circ(ax6, stk_params, slice=y_plane, plane='y', cmap=cmap_circ, plot_title=None, x_label='z',
                          y_label='x', text_present=text_present, interpolation=interpolation,
                          normalization=normalization, **kwargs)

    if cbars:
        cbax1 = fig.add_axes([0.91, 0.56, 0.04, 0.321])

        if cbax1_dir == 1:
            ph = np.ones((100, 100)) * np.linspace(1, -1, 100)[:, np.newaxis]
            I = np.ones((100, 100)) * np.linspace(0, 1, 100)[np.newaxis, :]
            imshow2d(np.array([ph, I]), ax=cbax1, cmap2d=cmap_lin, huevmin=-1, huevmax=1, lightvmin=0, lightvmax=1,
                     extent=[0, 1, np.pi / 2, -np.pi / 2], aspect='auto')
            plt.yticks(np.linspace(np.pi / 2, -np.pi / 2, 3),
                       ['$-\pi/2$', '0', '$\pi/2$'])  # ['0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$']
            cbax1.yaxis.set_label_position("right")
            cbax1.yaxis.tick_right()
            cbax1.set_ylabel('$\psi$')
            if normalization == 's0':
                cbax1.set_xlabel('$ \sqrt{S_1^2+S_2^2} / S_0$')
            elif normalization == 's0_max':
                cbax1.set_xlabel('$ \sqrt{S_1^2+S_2^2} / max(S_0)$')
            else:
                cbax1.set_xlabel('$ \sqrt{S_1^2+S_2^2}$')
            cbax1.tick_params(axis='both', which='major', labelsize=10)

        else:
            ph = np.ones((100, 100)) * np.linspace(1, -1, 100)[np.newaxis, :]
            I = np.ones((100, 100)) * np.linspace(0, 1, 100)[:, np.newaxis]
            imshow2d(np.array([ph, I]), ax=cbax1, cmap2d=cmap_lin, huevmin=-1, huevmax=1, lightvmin=0, lightvmax=1,
                     extent=[np.pi / 2, -np.pi / 2, 1, 0], aspect='auto')
            plt.xticks(np.linspace(-np.pi / 2, np.pi / 2, 3),
                       ['$-\pi/2$', '0', '$\pi/2$'])  # ['0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$']
            cbax1.yaxis.set_label_position("right")
            cbax1.yaxis.tick_right()
            cbax1.set_xlabel('$\psi$')
            if normalization == 's0':
                cbax1.set_ylabel('$ \sqrt{S_1^2+S_2^2} / S_0$')
            elif normalization == 's0_max':
                cbax1.set_ylabel('$ \sqrt{S_1^2+S_2^2} / max(S_0)$')
            else:
                cbax1.set_ylabel('$ \sqrt{S_1^2+S_2^2}$')
            cbax1.tick_params(axis='both', which='major', labelsize=10)

        cbax2 = fig.add_axes([0.91, 0.11, 0.04, 0.321])  # This is the position for the colorbar [x, y, width, height]
        cbar_circ_im = plt.colorbar(circular_plt, cax=cbax2)
        if normalization == 's0':
            cbax2.set_ylabel('$S_3 / S_0$')
        elif normalization == 's0_max':
            cbax2.set_ylabel('$S_3 / max(S_0)$')
        else:
            cbax2.set_ylabel('S3')
        cbax2.tick_params(axis='both', which='major', labelsize=10)

    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    _logger.info(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
    plt.draw()
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving figure')
        fig.savefig(savefig)

    if showfig:
        _logger.debug(ind_str + 'showing Stokes Parameters')
        plt.show()
    else:
        plt.close('all')


@if_plottable
def plot_stokes_sbfg_lin(ax, stk_params, slice, plane, cmap2d='brightwheel', plot_title=None, x_label='', y_label='',
                         result=0, text_present=True, interpolation=None, normalization='s0_max', **kwargs):
    """
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
    """
    # Getting intersections of stk_params for ploting data
    z_max, y_max, x_max = np.unravel_index(stk_params.s0.argmax(), stk_params.s0.shape)  # getting max element position

    if plane in ['x', 2]:
        swap_axes = True
        extent = [stk_params.sc_z[0] * 1e6, stk_params.sc_z[-1] * 1e6, stk_params.sc_y[0] * 1e6,
                  stk_params.sc_y[-1] * 1e6]
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
        extent = [stk_params.sc_z[0] * 1e6, stk_params.sc_z[-1] * 1e6, stk_params.sc_x[0] * 1e6,
                  stk_params.sc_x[-1] * 1e6]
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
        extent = [stk_params.sc_x[0] * 1e6, stk_params.sc_x[-1] * 1e6, stk_params.sc_y[0] * 1e6,
                  stk_params.sc_y[-1] * 1e6]
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
    if normalization is None:
        norm = 1
    elif normalization == 's0':
        norm = stk_params_plane.s0
    elif normalization == 's0_max':
        norm = np.amax(stk_params.s0)
    else:
        raise ValueError('"normalization" should be in [None, "s0", "s0_max"]')

    if swap_axes:
        lin_pol_plane = np.swapaxes((stk_params_plane.P_pol_l() / norm), 0, 1)
        psi_plane = np.swapaxes((2 * stk_params_plane.psi() / np.pi), 0, 1)
    else:
        lin_pol_plane = stk_params_plane.P_pol_l() / norm
        psi_plane = 2 * stk_params_plane.psi() / np.pi

    m, n = psi_plane.shape
    linear_plt = imshow2d(np.array([psi_plane, lin_pol_plane]), ax=ax, cmap2d=cmap2d, extent=extent,
                          interpolation=interpolation, aspect='auto', huevmin=-1, huevmax=1, lightvmin=0, lightvmax=1,
                          origin='lower', **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title, fontsize=15)
    ax.set_xlabel(x_label + ' [$\mu$m]')
    ax.set_ylabel(y_label + ' [$\mu$m]')
    if text_present:
        # print(plane, slice_pos*1e6)
        dic = {'proj': 'projection', 'max_slice': 'slice at {:}={:.3f} $\mu$m'}
        ax.text(0.97, 0.97, dic.get(slice, 'slice at {:}={:.3f} $\mu$m (max int)').format(plane, slice_pos * 1e6),
                horizontalalignment='right',
                verticalalignment='top', transform=ax.transAxes, fontsize=10)
    if result:
        return linear_plt


@if_plottable
def plot_stokes_sbfg_circ(ax, stk_params, slice, plane, cmap='seismic', plot_title=None, x_label='', y_label='',
                          result=0, text_present=True, interpolation=None, normalization='s0_max', **kwargs):
    """
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
    """

    # Getting intersections of stk_params for ploting data
    z_max, y_max, x_max = np.unravel_index(stk_params.s0.argmax(), stk_params.s0.shape)  # getting max element position

    if plane in ['x', 2]:
        swap_axes = True
        extent = [stk_params.sc_z[0] * 1e6, stk_params.sc_z[-1] * 1e6, stk_params.sc_y[0] * 1e6,
                  stk_params.sc_y[-1] * 1e6]
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
        extent = [stk_params.sc_z[0] * 1e6, stk_params.sc_z[-1] * 1e6, stk_params.sc_x[0] * 1e6,
                  stk_params.sc_x[-1] * 1e6]
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
        extent = [stk_params.sc_x[0] * 1e6, stk_params.sc_x[-1] * 1e6, stk_params.sc_y[0] * 1e6,
                  stk_params.sc_y[-1] * 1e6]
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
    if normalization is None:
        norm = 1
    elif normalization == 's0':
        norm = stk_params_plane.s0
    elif normalization == 's0_max':
        norm = np.amax(stk_params.s0)
    else:
        raise ValueError('"normalization" should be in [None, "s0", "s0_max"]')

    if swap_axes:
        s3_plane = np.swapaxes((stk_params_plane.s3 / norm), 0, 1)
    else:
        s3_plane = stk_params_plane.s3 / norm

    m, n = s3_plane.shape
    circular_plt = ax.imshow(s3_plane, cmap=cmap, vmin=-1, vmax=1, interpolation=interpolation, aspect='auto',
                             extent=extent, origin='lower', **kwargs)
    if plot_title is not None:
        ax.set_title(plot_title, fontsize=15)
    ax.set_xlabel(x_label + ' [$\mu$m]')
    ax.set_ylabel(y_label + ' [$\mu$m]')
    if text_present:
        dic = {'proj': 'projection', 'max_slice': 'slice at {:}={:.3f} $\mu$m (max int)'}
        ax.text(0.97, 0.97, dic.get(slice, 'slice at {:}={:.3f} $\mu$m').format(plane, slice_pos * 1e6),
                horizontalalignment='right',
                verticalalignment='top', transform=ax.transAxes, fontsize=10)
    if result:
        return circular_plt


def plot_hprofile(*args, **kwargs):
    plot_1d_hprofile(*args, **kwargs)


def plot_1d_hprofile(height_profile, figsize=4, fig_name='Height profile', savefig=False, showfig=True,
                     **kwargs):
    """
    This function plotting the height map and PSD of HeightProfile object

    :param height_profile: HeightProfile object from ocelot
    :param figsize: size of figure
    :param fig_name: name of figure
    :param savefig: bool type flag, responding for saving figure
    :param showfig: bool type flag, responding for showing figure
    :param kwargs:
    """
    if (showfig == False) and (savefig == False):
        return

    _logger.info('plotting height_profile')
    _logger.warning(ind_str + 'in beta')
    start_time = time.time()

    # Plotting data
    fig = plt.figure(fig_name)
    # fig.canvas.set_window_title(fig_name)
    fig.clf()
    fig.set_size_inches((4 * figsize, 1.5 * figsize), forward=True)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(height_profile.s * 1000, height_profile.h, **kwargs)
    ax1.set_title('height errors', fontsize=15)
    ax1.set_xlabel('s [mm]')
    ax1.set_ylabel('height [m]')

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.loglog(*height_profile.psd(), marker='*', **kwargs)
    ax2.set_title('power spectral density', fontsize=15)
    ax2.set_xlabel('k [1/m]')
    ax2.set_ylabel('PSD [m^3]')
    ax2.text(0.97, 0.97, 'surface height errors RMS = {0:.02e} [m]'.format(height_profile.hrms()),
             horizontalalignment='right',
             verticalalignment='top', transform=ax2.transAxes, fontsize=12)

    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.draw()
    _logger.debug(ind_str + 'done in {:.2f} seconds'.format(time.time() - start_time))
    if savefig != False:
        if savefig == True:
            savefig = 'png'
        _logger.debug(ind_str + 'saving figure')
        fig.savefig(fig_name + '.' + savefig)

    if showfig:
        _logger.debug('showing HeightProfile')
        plt.show()
    else:
        plt.close('all')
