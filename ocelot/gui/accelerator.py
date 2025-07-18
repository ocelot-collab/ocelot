"""
user interface for viewing/editing electron optics layouts
"""

from ocelot.cpbd.physics_proc import *
from scipy import stats
from ocelot.cpbd.beam import global_slice_analysis
import sys
import os
import csv
import numbers
import matplotlib
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PatchCollection
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ocelot.cpbd.optics import *
import numpy as np
from ocelot.cpbd.elements import *
from copy import deepcopy


import matplotlib.font_manager as font_manager


def plot_lattice(lat, axis, alpha=1.0, params={'kmax': 2.0, 'ang_max': 0.5e-2}, s_start=0.0):
    axis.grid(True)

    pos = s_start
    offs = np.array([0, 0.0])

    ang_max = params['ang_max']  # max dipole strength in lattice
    min_dipole_height = 0.1  # dipole of any strength will show at least this strength

    min_solenoid_height = 0.1
    sol_max = 0.1

    kmax = params['kmax']
    min_quad_height = 0.1

    rendered_seq = []
    rendered_len = 0.0
    total_len = 0.0

    for e in lat.sequence:

        # if e.type in ['bend','sbend', 'rbend', 'quadrupole', 'undulator', 'drift', 'monitor','hcor','vcor', 'cavity','edge', 'solenoid']:
        if e.__class__ in [Bend, SBend, RBend, Quadrupole, Undulator, Drift, Monitor, Hcor, Vcor, Cavity,
                           Solenoid]:
            e.s = total_len
            rendered_seq.append(e)
            rendered_len += e.l
        total_len += e.l

    for e in rendered_seq:
        dx = e.l

        if e.__class__ in [Bend, SBend, RBend, Hcor, Vcor]:
            axis.add_patch(mpatches.Rectangle(offs + np.array([pos, 0.0]), dx,
                                              np.sign(-e.angle) * min_dipole_height - e.angle / ang_max * (
                1 - min_dipole_height),
                color='#0099FF', alpha=alpha))

        if e.__class__ in [Solenoid]:
            axis.add_patch(mpatches.Rectangle(offs + np.array([pos, 0.0]), dx,
                                              np.sign(-e.k) * min_solenoid_height - e.k / sol_max * (
                1 - min_solenoid_height),
                color='#FF99FF', alpha=alpha))

        if e.__class__ in [Quadrupole]:

            if e.k1 >= 0:
                axis.add_patch(mpatches.Ellipse(offs + np.array([pos, 0.0]), dx, min_quad_height + abs(e.k1 / kmax) * 2,
                                                color='green', alpha=alpha))
            else:
                Path = mpath.Path
                h = abs(e.k1 / kmax) + min_quad_height / 2

                verts = np.array([
                    (dx, h),
                    (-dx, h),
                    (-dx / 4, 0),
                    (-dx, -h),
                    (dx, -h),
                    (dx / 4, 0),
                    (dx, h)
                ])

                codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.LINETO, Path.LINETO, Path.CURVE3, Path.CURVE3]

                path = mpath.Path(verts + offs + np.array([pos, 0.0]), codes)
                patch = mpatches.PathPatch(path, color='green', alpha=alpha)
                axis.add_patch(patch)

        if e.__class__ in [Undulator]:
            nper = 16
            dxs = dx / nper / 2.0

            height = 0.7
            gap = 0.1
            if e.Kx ** 2 + e.Ky ** 2 < 1.e-6:
                gap = 0.2

            for iseg in np.arange(0, nper, 2):
                axis.add_patch(
                    mpatches.Rectangle(offs + np.array([pos + iseg * dxs, gap]), dxs, height, color='red', alpha=alpha))
                axis.add_patch(
                    mpatches.Rectangle(offs + np.array([pos + iseg * dxs, -height - gap]), dxs, height, color='blue',
                                       alpha=alpha))
                axis.add_patch(
                    mpatches.Rectangle(offs + np.array([pos + (iseg + 1) * dxs, gap]), dxs, height, color='blue',
                                       alpha=alpha))
                axis.add_patch(mpatches.Rectangle(offs + np.array([pos + (iseg + 1) * dxs, -height - gap]), dxs, height,
                                                  color='red', alpha=alpha))

        if e.__class__ in [Cavity]:
            nper = 16
            dxs = dx / nper / 2.0

            height = 0.7
            gap = 0.1
            axis.add_patch(mpatches.Ellipse(offs + np.array([pos + dx / 2, 0.0]),
                                            dx / 3, 0.5,
                                            color='#FF0033', alpha=alpha))

            axis.add_patch(mpatches.Ellipse(offs + np.array([pos + dx / 6, 0.0]),
                                            dx / 3, 0.5,
                                            color='#FF0033', alpha=alpha))

            axis.add_patch(mpatches.Ellipse(offs + np.array([pos + 5 * dx / 6, 0.0]),
                                            dx / 3, 0.5,
                                            color='#FF0033', alpha=alpha))

        if e.__class__ in [Monitor]:
            Path = mpath.Path
            h = 0.005
            offs_mon = np.array([0.0, 0.03])

            verts = np.array([
                (0, h),
                (0, -h),
                (h, 0),
                (-h, 0)
            ])

            codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]

            verts += offs_mon

            path = mpath.Path(verts + offs + np.array([pos, 0.0]), codes)
            patch = mpatches.PathPatch(path, color='black', lw=2, alpha=alpha * 0.5)
            axis.add_patch(patch)
            axis.add_patch(
                mpatches.Circle(offs + offs_mon + np.array([pos, 0.0]), h / 2, color='black', alpha=alpha * 0.5))

        pos += dx

    lw = 0.005
    axis.add_patch(mpatches.Rectangle(offs - (0, lw / 2), 1.0, lw, color='black', alpha=alpha * 0.4))

    # axis.set_ylim(-1,1)
    # axis.set_xticks([])
    # axis.set_yticks([])


dict_plot = {Quadrupole: {"scale": 0.7, "color": "r", "edgecolor": "r", "label": "quad"},
             Sextupole: {"scale": 0.5, "color": "g", "edgecolor": "g", "label": "sext"},
             Octupole: {"scale": 0.5, "color": "g", "edgecolor": "g", "label": "oct"},
             Cavity: {"scale": 0.7, "color": "orange", "edgecolor": "lightgreen", "label": "cav"},
             TWCavity: {"scale": 0.7, "color": "orange", "edgecolor": "lightgreen", "label": "twcav"},
             Bend: {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k", "label": "bend"},
             RBend: {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k", "label": "bend"},
             SBend: {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k", "label": "bend"},
             Matrix: {"scale": 0.7, "color": "pink", "edgecolor": "k", "label": "mat"},
             Multipole: {"scale": 0.7, "color": "g", "edgecolor": "k", "label": "mult"},
             Undulator: {"scale": 0.7, "color": "pink", "edgecolor": "k", "label": "und"},
             Monitor: {"scale": 0.5, "color": "orange", "edgecolor": "orange", "label": "mon"},
             Hcor: {"scale": 0.7, "color": "c", "edgecolor": "c", "label": "cor"},
             Vcor: {"scale": 0.7, "color": "c", "edgecolor": "c", "label": "cor"},
             Drift: {"scale": 0., "color": "k", "edgecolor": "k", "label": ""},
             Marker: {"scale": 0., "color": "k", "edgecolor": "k", "label": "mark"},
             Solenoid: {"scale": 0.7, "color": "g", "edgecolor": "g", "label": "sol"},
             TDCavity: {"scale": 0.7, "color": "magenta", "edgecolor": "g", "label": "tds"},
             UnknownElement: {"scale": 0.7, "color": "g", "edgecolor": "g", "label": "unk"},
             XYQuadrupole: {"scale": 0.7, "color": "r", "edgecolor": "r", "label": "xyquad"},
             Aperture: {"scale": 0.7, "color": "g", "edgecolor": "g", "label": "ap"},
             }


def plot_elems(fig, ax, lat, s_point=0, y_lim=None, y_scale=1, legend=True, font_size=18, excld_legend=None):
    legend_font_size = font_size

    if excld_legend is None:
        excld_legend = []

    dict_copy = deepcopy(dict_plot)
    alpha = 1
    ax.set_ylim((-1, 1.5))
    ax.tick_params(axis='both', labelsize=font_size)
    if y_lim is not None:
        ax.set_ylim(y_lim)
    points_with_annotation = []
    L = 0.
    q = []
    b = []
    c = []
    s = []
    u = []
    rf = []
    m = []
    sol = []
    for elem in lat.sequence:
        if elem.__class__ == Quadrupole:
            q.append(elem.k1)
        elif elem.__class__ in [Bend, RBend, SBend]:
            b.append(elem.angle)
        elif elem.__class__ in [Hcor, Vcor]:
            c.append(elem.angle)
        elif elem.__class__ == Sextupole:
            s.append(elem.k2)
        elif elem.__class__ == Solenoid:
            sol.append(elem.k)
        elif elem.__class__ == Undulator:
            u.append(elem.Kx + elem.Ky)
        elif elem.__class__ in [Cavity, TWCavity, TDCavity]:
            rf.append(elem.v)
        elif elem.__class__ == Multipole:
            m.append(sum(np.abs(elem.kn)))

    q_max = np.max(np.abs(q)) if len(q) != 0 else 0
    b_max = np.max(np.abs(b)) if len(b) != 0 else 0
    s_max = np.max(np.abs(s)) if len(s) != 0 else 0
    c_max = np.max(np.abs(c)) if len(c) != 0 else 0
    u_max = np.max(np.abs(u)) if len(u) != 0 else 0
    sol_max = np.max(np.abs(sol)) if len(sol) != 0 else 0
    rf_max = np.max(np.abs(rf)) if len(rf) != 0 else 0
    m_max = np.max(m) if len(m) != 0 else 0
    ncols = np.sign(len(q)) + np.sign(len(b)) + np.sign(len(s)) + np.sign(len(c)) + np.sign(len(u)) + np.sign(
        len(rf)) + np.sign(len(m)) + np.sign(len(sol))

    labels_dict = {}
    for elem in dict_copy.keys():
        labels_dict[elem] = dict_copy[elem]["label"]
    for elem in lat.sequence:
        if elem.__class__ in excld_legend:
            elem = Drift(l=elem.l)

        if elem.__class__ == Marker:
            L += elem.l
            continue
        l = elem.l
        if l == 0:
            l = 0.03
        # type = elem.type
        if elem.__class__ in dict_copy:
            scale = dict_copy[elem.__class__]["scale"]
            color = dict_copy[elem.__class__]["color"]
            label = dict_copy[elem.__class__]["label"]
            ecolor = dict_copy[elem.__class__]["edgecolor"]
        else:
            scale = dict_copy[UnknownElement]["scale"]
            color = dict_copy[UnknownElement]["color"]
            label = dict_copy[UnknownElement]["label"]
            ecolor = dict_copy[UnknownElement]["edgecolor"]
        ampl = 1
        s_coord = np.array(
            [L + elem.l / 2. - l / 2., L + elem.l / 2. - l / 2., L + elem.l / 2. + l / 2., L + elem.l / 2. + l / 2.,
             L + elem.l / 2. - l / 2.]) + s_point

        rect = np.array([-1, 1, 1, -1, -1])

        if elem.__class__ == Quadrupole:
            ampl = elem.k1 / q_max if q_max != 0 else 1
            point, = ax.fill(s_coord, (rect + 1) * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Solenoid:
            ampl = elem.k / sol_max if sol_max != 0 else 1
            point, = ax.fill(s_coord, (rect + 1) * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ in [Bend, RBend, SBend]:
            ampl = elem.angle / b_max if b_max != 0 else 1
            point, = ax.fill(s_coord, (rect + 1) * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[Bend]["label"] = ""
            dict_copy[RBend]["label"] = ""
            dict_copy[SBend]["label"] = ""

        elif elem.__class__ in [Hcor, Vcor]:
            ampl = elem.angle / c_max if c_max != 0 else 0.5
            if elem.angle == 0:
                ampl = 0.5
                point, = ax.fill(s_coord, rect * ampl * scale * y_scale, "lightcyan", edgecolor="k",
                                 alpha=0.5, label=label)
            else:
                point, = ax.fill(s_coord, (rect + 1) * ampl * scale * y_scale, color, edgecolor=ecolor,
                                 alpha=alpha, label=label)
            dict_copy[Hcor]["label"] = ""
            dict_copy[Vcor]["label"] = ""

        elif elem.__class__ == Sextupole:
            ampl = elem.k2 / s_max if s_max != 0 else 1
            point, = ax.fill(s_coord, (rect + 1) * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ in [Cavity, TWCavity, TDCavity]:
            ampl = 1
            point, = ax.fill(s_coord, rect * ampl * scale * y_scale, color,
                             alpha=alpha, edgecolor=ecolor,
                             label=label)
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Undulator:
            ampl = elem.Kx / u_max if u_max != 0 else 0.5
            point, = ax.fill(s_coord, rect * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Multipole:
            ampl = sum(elem.kn) / m_max if u_max != 0 else 0.5
            point, = ax.fill(s_coord, rect * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha, label=label)
            dict_copy[elem.__class__]["label"] = ""

        else:
            point, = ax.fill(s_coord, rect * ampl * scale * y_scale, color, edgecolor=ecolor,
                             alpha=alpha)

        annotation = ax.annotate(elem.__class__.__name__ + ": " + elem.id,
                                 xy=(L + l / 2. + s_point, 0),  # xycoords='data',
                                 # xytext=(i + 1, i), textcoords='data',
                                 horizontalalignment="left",
                                 arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=+0.2"),
                                 bbox=dict(boxstyle="round", facecolor="w", edgecolor="0.5", alpha=0.9),
                                 fontsize=legend_font_size
                                 )
        # by default, disable the annotation visibility
        annotation.set_visible(False)
        L += elem.l
        points_with_annotation.append([point, annotation])
        ax.set_xlabel("s [m]", fontsize=font_size)

    def on_move(event):
        visibility_changed = False
        for point, annotation in points_with_annotation:
            should_be_visible = (point.contains(event)[0] is True)
            if should_be_visible != annotation.get_visible():
                visibility_changed = True
                annotation.set_visible(should_be_visible)

        if visibility_changed:
            plt.draw()

    on_move_id = fig.canvas.mpl_connect('motion_notify_event', on_move)
    if legend:
        ax.legend(loc='upper center', ncol=ncols, shadow=False, prop=font_manager.FontProperties(size=legend_font_size))


def plot_disp(ax, tws, top_plot, font_size):
    S = [p.s for p in tws]  # map(lambda p:p.s, tws)
    d_Ftop = []
    Fmin = []
    Fmax = []
    for elem in top_plot:
        Ftop = [getattr(p, elem) for p in tws]

        Fmin.append(min(Ftop))
        Fmax.append(max(Ftop))
        greek = ""
        if "beta" in elem or "alpha" in elem or "mu" in elem:
            greek = "\\"
        if "mu" in elem:
            elem = elem.replace("mu", "mu_")
        top_label = r"$" + greek + elem + "$"
        ax.plot(S, Ftop, lw=2, label=top_label)
        d_Ftop.append(max(Ftop) - min(Ftop))
    d_F = max(d_Ftop)
    if d_F == 0:
        d_Dx = 1
        ax.set_ylim((min(Fmin) - d_Dx * 0.1, max(Fmax) + d_Dx * 0.1))
    if top_plot[0] == "E":
        top_ylabel = r"$" + "/".join(top_plot) + "$" + ", [GeV]"
    elif top_plot[0] in ["mux", 'muy']:
        top_ylabel = r"$" + "/".join(top_plot) + "$" + ", [rad]"
    else:
        top_ylabel = r"$" + "/".join(top_plot) + "$" + ", [m]"

    # yticks = ax.get_yticks()
    # yticks = yticks[2::1]
    # ax.set_yticks(yticks)

    ax.yaxis.set_major_locator(plt.MaxNLocator(3))

    ax.set_ylabel(top_ylabel, fontsize=font_size)
    ax.tick_params(axis='both', labelsize=font_size)

    leg2 = ax.legend(loc='upper right', shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg2.get_frame().set_alpha(0.2)


def plot_betas(ax, S, beta_x, beta_y, font_size):
    ax.set_ylabel(r"$\beta_{x,y}$ [m]", fontsize=font_size)
    ax.plot(S, beta_x, 'b', lw=2, label=r"$\beta_{x}$")
    ax.plot(S, beta_y, 'r', lw=2, label=r"$\beta_{y}$")
    ax.tick_params(axis='both', labelsize=font_size)
    leg = ax.legend(loc='upper left', shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.2)


def plot_opt_func(lat, tws, top_plot=["Dx"], legend=True, fig_name=None, grid=True, font_size=12, excld_legend=None,
                  figsize=None, plt_figure=None):
    """
    function for plotting: lattice (bottom section), vertical and horizontal beta-functions (middle section),
    other parameters (top section) such as "Dx", "Dy", "E", "mux", "muy", "alpha_x", "alpha_y", "gamma_x", "gamma_y"

    :param lat: MagneticLattice,
    :param tws: list if Twiss objects,
    :param top_plot:  ["Dx"] - parameters which displayed in top section. Can be any attribute of Twiss class, e.g. top_plot=["Dx", "Dy", "alpha_x"]
    :param legend: True - displaying legend of element types in bottom section,
    :param fig_name: None - name of figure,
    :param grid: True - grid
    :param font_size: 16 - font size for any element of plot
    :param excld_legend: None, exclude type of element from the legend, e.g. excld_legend=[Hcor, Vcor]
    :param figsize: None or e.g. (8, 6)
    :param plt_figure: None, exported plt.figure()
    :return:
    """
    if plt_figure is None:
        if fig_name is None:
            fig = plt.figure(figsize=figsize)
        else:
            fig = plt.figure(fig_name, figsize=figsize)
    else:
        fig = plt_figure

    plt.rc('axes', grid=grid)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85

    rect1 = [left, 0.65, width, 0.3]
    rect2 = [left, 0.19, width, 0.46]
    rect3 = [left, 0.07, width, 0.12]

    ax_top = fig.add_axes(rect1)
    ax_b = fig.add_axes(rect2, sharex=ax_top)  # left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_top)
    for ax in ax_b, ax_el, ax_top:
        if ax != ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    ax_b.grid(grid)
    ax_top.grid(grid)
    ax_el.set_yticks([])
    ax_el.grid(grid)

    fig.subplots_adjust(hspace=0)
    beta_x = [p.beta_x for p in tws]
    beta_y = [p.beta_y for p in tws]
    S = [p.s for p in tws]

    plt.xlim(S[0], S[-1])

    plot_disp(ax_top, tws, top_plot, font_size)

    plot_betas(ax_b, S, beta_x, beta_y, font_size)
    plot_elems(fig, ax_el, lat, s_point=S[0], legend=legend, y_scale=0.8, font_size=font_size,
               excld_legend=excld_legend)


def plot_opt_func_reduced(lat, tws, top_plot=["Dx"], legend=True, fig_name=None, grid=False, font_size=18):
    """
    function for plotting: lattice (bottom section), vertical and horizontal beta-functions (middle section),
    other parameters (top section) such as "Dx", "Dy", "E", "mux", "muy", "alpha_x", "alpha_y", "gamma_x", "gamma_y"
    lat - MagneticLattice,
    tws - list if Twiss objects,
    top_plot=["Dx"] - parameters which displayed in top section. Example top_plot=["Dx", "Dy", "alpha_x"]
    legend=True - displaying legend of element types in bottom section,
    fig_name=None - name of figure,
    grid=True - grid
    font_size=18 - font size.
    """
    if fig_name is None:
        fig = plt.figure()
    else:
        fig = plt.figure(fig_name)

    plt.rc('axes', grid=grid)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85

    rect1 = [left, 0.63, width, 0.3]
    rect2 = [left, 0.20, width, 0.43]
    rect3 = [left, 0.10, width, 0.10]

    ax_top = fig.add_axes(rect1)
    ax_b = fig.add_axes(rect2, sharex=ax_top)  # left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_top)
    for ax in ax_b, ax_el, ax_top:
        if ax != ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    ax_b.grid(grid)
    ax_top.grid(grid)
    ax_el.set_yticks([])
    ax_el.grid(grid)

    fig.subplots_adjust(hspace=0)
    beta_x = [p.beta_x for p in tws]  # list(map(lambda p:p.beta_x, tws))
    beta_y = [p.beta_y for p in tws]  # list(map(lambda p:p.beta_y, tws))
    D_x = [p.Dx for p in tws]  # list(map(lambda p:p.beta_x, tws))
    D_y = [p.Dy for p in tws]  # list(map(lambda p:p.beta_y, tws))
    S = [p.s for p in tws]  # list(map(lambda p:p.s, tws))

    plt.xlim(S[0], S[-1])

    ax_top.set_ylabel(r"$D_{x,y}$ [m]")
    ax_top.plot(S, D_x, 'b', lw=2, label=r"$D_{x}$")
    ax_top.plot(S, D_y, 'r', lw=2, label=r"$D_{y}$")
    leg = ax_top.legend(loc='upper left', shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.2)

    plot_betas(ax_b, S, beta_x, beta_y, font_size)
    plot_elems(fig, ax_el, lat, s_point=S[0], legend=legend, y_scale=0.8)


def plot_xy(ax, S, X, Y, font_size):
    ax.set_ylabel(r"$X, Y$, m")
    ax.plot(S, X, 'r', lw=2, label=r"$X$")
    ax.plot(S, Y, 'b', lw=2, label=r"$Y$")
    leg = ax.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.5)


def plot_API(lat, legend=True, fig_name=1, grid=True, font_size=12,
             excld_legend=None, figsize=None, add_extra_subplot=False):
    """
    Creates a figure with a main plot area (ax_xy) and beamline (ax_el).
    Optionally adds an extra subplot (ax_extra) above ax_xy.

    :param lat: MagneticLattice
    :param legend: Show legend on beamline
    :param fig_name: ID for the matplotlib figure
    :param grid: Enable grid on plots
    :param font_size: Font size for axis labels
    :param excld_legend: Element types to exclude from beamline legend
    :param figsize: Custom figure size
    :param add_extra_subplot: If True, includes an additional subplot above ax_xy
    :return: fig, axes (ax_xy, or (ax_extra, ax_xy))
    """

    fig = plt.figure(fig_name, figsize=figsize)
    plt.rc('axes', grid=grid)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)

    left, width = 0.1, 0.85

    if add_extra_subplot:
        rect1 = [left, 0.59, width, 0.34]  # ax_extra
        rect2 = [left, 0.2, width, 0.34]  # ax_xy
        rect3 = [left, 0.07, width, 0.10]  # ax_el
        ax_extra = fig.add_axes(rect1)
        ax_xy = fig.add_axes(rect2, sharex=ax_extra)
        ax_el = fig.add_axes(rect3, sharex=ax_extra)

        for ax in [ax_extra, ax_xy]:
            for label in ax.get_xticklabels():
                label.set_visible(False)

        axes = (ax_extra, ax_xy)

    else:
        rect2 = [left, 0.19, width, 0.69]  # ax_xy
        rect3 = [left, 0.07, width, 0.10]  # ax_el
        ax_xy = fig.add_axes(rect2)
        ax_el = fig.add_axes(rect3, sharex=ax_xy)

        for label in ax_xy.get_xticklabels():
            label.set_visible(False)

        axes = ax_xy

    # General settings
    for ax in ([ax_extra, ax_xy] if add_extra_subplot else [ax_xy]):
        ax.grid(grid)
        ax.tick_params(axis='both', labelsize=font_size)

    ax_el.set_yticks([])
    ax_el.grid(grid)
    ax_el.tick_params(axis='both', labelsize=font_size)

    fig.subplots_adjust(hspace=0)

    # Plot beamline (always present)
    plot_elems(fig, ax_el, lat, legend=legend, y_scale=0.8, font_size=font_size, excld_legend=excld_legend)
    return fig, axes


def compare_betas(lat, tws1, tws2, prefix1="beam1", prefix2="beam2", legend=True, fig_name=None, grid=True,
                  font_size=18):
    """
    function for plotting: lattice (bottom section), vertical and horizontal beta-functions (middle section),
    other parameters (top section) such as "Dx", "Dy", "E", "mux", "muy", "alpha_x", "alpha_y", "gamma_x", "gamma_y"
    lat - MagneticLattice,
    tws - list if Twiss objects,
    top_plot=["Dx"] - parameters which displayed in top section. Example top_plot=["Dx", "Dy", "alpha_x"]
    legend=True - displaying legend of element types in bottom section,
    fig_name=None - name of figure,
    grid=True - grid
    font_size=18 - font size.
    """
    # fig, ax_xy = plot_API(lat, legend=legend, fig_name=fig_name, grid=grid)
    if fig_name is None:
        fig = plt.figure()
    else:
        fig = plt.figure(fig_name)

    plt.rc('axes', grid=grid)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85

    rect1 = [left, 0.55, width, 0.37]
    rect2 = [left, 0.18, width, 0.37]
    rect3 = [left, 0.05, width, 0.13]

    ax_top = fig.add_axes(rect1)
    ax_b = fig.add_axes(rect2, sharex=ax_top)  # left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_top)
    for ax in ax_b, ax_el, ax_top:
        if ax != ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    ax_b.grid(grid)
    ax_top.grid(grid)
    ax_el.set_yticks([])
    ax_el.grid(grid)

    fig.subplots_adjust(hspace=0)
    beta_x1 = [p.beta_x for p in tws1]
    beta_y1 = [p.beta_y for p in tws1]
    beta_x2 = [p.beta_x for p in tws2]
    beta_y2 = [p.beta_y for p in tws2]
    s1 = [p.s for p in tws1]
    s2 = [p.s for p in tws2]
    # plt.plot(S, beta_x)

    plt.xlim(s1[0], s1[-1])

    # plot_disp(ax_top, tws, top_plot, font_size)
    ax_top.plot(s1, beta_y1, 'r', lw=2, label=prefix1 + r"  $\beta_{y}$")
    ax_top.plot(s2, beta_y2, 'b--', lw=2, label=prefix2 + r"  $\beta_{y}$")

    leg = ax_top.legend(shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.2)

    ax_b.set_ylabel(r"$\beta_{x,y}$, m")
    ax_b.plot(s1, beta_x1, 'r', lw=2, label=prefix1 + r"  $\beta_{x}$")

    ax_b.plot(s2, beta_x2, 'b--', lw=2, label=prefix2 + r"  $\beta_{x}$")
    leg = ax_b.legend(shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.2)

    plot_elems(fig, ax_el, lat, s_point=s1[0], legend=legend, y_scale=0.8)


def resonance(Qx, Qy, order=5):
    ORD = order
    qx1, qy1 = 0, 0
    qx2, qy2 = 2, 2
    X = []
    Y = []
    Order = []
    params = []
    # Qx = 0.22534
    # Qy = 0.301
    for order in range(1, ORD + 1):
        n = np.arange(-order, order + 1)
        m = np.array([order - abs(n), - order + abs(n)]).flatten()
        n = np.tile(n, 2)
        ziped = []
        for n, m in zip(n, m):
            if (n, m) in ziped:
                continue
            ziped.append((n, m))
        # ziped =  zip(n,m)
        # print ziped
        for p in range(-50, 50):
            for n, m in ziped:
                # print p, n,m
                if m != 0:
                    x = [qx1, qx2]
                    y = [(p - n * qx1) / float(m), (p - n * qx2) / float(m)]
                else:
                    x = [p / float(n), p / float(n)]
                    y = [qy1, qy2]
                params.append([n, m, p])
                X.append(x)
                Y.append(y)
                Order.append(order)
    return X, Y, Order, params


def plot_resonance_diag(ax, Qx, Qy, order):
    X, Y, Order, params = resonance(Qx, Qy, order)
    indsort = np.argsort(Order)

    X = np.array(X)

    for i, order in enumerate(Order):
        if order == 1:
            color = "k"
            lw = 3
        elif order == 2:
            color = "r"
            lw = 2
        elif order == 3:
            color = "b"
            lw = 2
        elif order == 4:
            color = "g"
            lw = 0.6
        else:  # order == 5:
            color = "c"
            lw = 0.3

        plt.plot(np.array(X[i]) + Qx, np.array(Y[i]) + Qy, color, lw=lw, picker=True)
        plt.xlim(Qx, Qx + 1)
        plt.ylim(Qy, Qy + 1)


def show_da(out_da, x_array, y_array, title=""):
    nx = len(x_array)
    ny = len(y_array)
    out_da = out_da.reshape(ny, nx)
    xmin, xmax, ymin, ymax = np.min(x_array), np.max(x_array), np.min(y_array), np.max(y_array)
    extent = xmin, xmax, ymin, ymax

    plt.figure(figsize=(10, 7))
    fig1 = plt.contour(out_da, linewidths=2, extent=extent)  # , colors = 'r')

    plt.grid(True)
    plt.title(title)
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    cb = plt.colorbar()
    cb.set_label('Nturns')

    plt.show()


def show_mu(contour_da, mux, muy, x_array, y_array, zones=None):
    from matplotlib import pyplot as plt

    nx = len(x_array)
    ny = len(y_array)
    t = np.linspace(0, 3.14, num=100)
    contour_da = contour_da.reshape(ny, nx)
    mux = mux.reshape(ny, nx)
    muy = muy.reshape(ny, nx)
    xmin, xmax, ymin, ymax = min(x_array), max(x_array), min(y_array), max(y_array)
    plt.figure(1, figsize=(10, 7))  # axisbg='darkslategray'
    extent = xmin, xmax, ymin, ymax

    my_cmap = plt.cm.Paired

    fig1 = plt.contour(contour_da, 1, extent=extent, linewidths=2, colors='k')  # , colors = 'r')
    fig1 = plt.contourf(mux, 40, cmap=my_cmap, extent=extent)  # , colors = 'r')
    cb = plt.colorbar()
    fig1 = plt.contourf(mux, 10, levels=[-1, -.0001], colors='w', extent=extent)
    if zones is not None:
        x_zone = zones[0]
        y_zone = zones[1]
        plt.plot(x_zone * np.cos(t), y_zone * np.sin(t), "g", lw=2)
        plt.plot(2 * x_zone * np.cos(t), 2 * y_zone * np.sin(t), "b", lw=2)
        plt.plot(3 * x_zone * np.cos(t), 3 * y_zone * np.sin(t), "r", lw=2)
        plt.plot(4 * x_zone * np.cos(t), 4 * y_zone * np.sin(t), "y", lw=2)
    plt.grid(True)
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    cb.set_label('Qx')
    plt.figure(2, figsize=(10, 7))

    fig1 = plt.contour(contour_da, 1, extent=extent, linewidths=2, colors='k')  # , colors = 'r')
    fig1 = plt.contourf(muy, 40, cmap=my_cmap, extent=extent)  # , colors = 'r')
    if zones is not None:
        x_zone = zones[0]
        y_zone = zones[1]
        plt.plot(x_zone * np.cos(t), y_zone * np.sin(t), "g", lw=2)
        plt.plot(2 * x_zone * np.cos(t), 2 * y_zone * np.sin(t), "b", lw=2)
        plt.plot(3 * x_zone * np.cos(t), 3 * y_zone * np.sin(t), "r", lw=2)
        plt.plot(4 * x_zone * np.cos(t), 4 * y_zone * np.sin(t), "y", lw=2)

    cb = plt.colorbar()
    fig1 = plt.contourf(muy, 10, levels=[-1, -.0001], colors='w', extent=extent)
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    plt.grid(True)
    cb.set_label('Qy')
    plt.show()


def show_density(
    x: np.ndarray,
    y: np.ndarray,
    q: np.ndarray | None = None,
    *,
    ax: plt.Axes | None = None,
    nbins_x: int = 250,
    nbins_y: int = 250,
    interpolation: str = "bilinear",
    xlabel: str | None = None,
    ylabel: str | None = None,
    nfig: int | None = 50,
    title: str | None = None,
    figsize: tuple[int, int] | None = None,
    grid: bool = True,
    show_xtick_label: bool = True,
    limits: list[list[float]] | None = None,
    cmap: str = "my_rainbow",
    density_mode: str = "hist",  # "hist" for imshow; "scatter" for per-particle colouring
    use_log: bool = True,
    downsample: int | None = None,  # new: for scatter mode, use N random particles
    seed: int | None = None
):
    """Visualise 2-D charge density or particle distribution.

    Parameters
    ----------
    x, y : np.ndarray
        Coordinates of macro-particles (same length).
    q : np.ndarray | None, optional
        Charge (weight) of each particle. If *None*, all particles get equal
        charge so that total charge = 1 (arbitrary units).
    ax : matplotlib.axes.Axes | None, optional
        Axis to plot into. If *None*, a new figure/axis is created.
    nbins_x, nbins_y : int
        Number of histogram bins in *x* and *y* directions.
    interpolation : str
        Interpolation method for *imshow* (only for ``mode='hist'``).
    xlabel, ylabel : str | None
        Axis labels.
    nfig : int | None
        Figure number (only if *ax* is *None*).
    title : str | None
        Figure title.
    figsize : tuple[int, int] | None
        Figure size in inches (only if *ax* is *None*).
    grid : bool
        Turn grid on/off.
    show_xtick_label : bool
        Show or hide x tick labels (useful for subplots).
    limits : [[xmin, xmax], [ymin, ymax]] | None
        Plot limits. If *None*, determined from data with ~5 % margin in *y*.
    cmap : str
        Matplotlib colormap name. Special value ``"my_rainbow"`` uses
        ``matplotlib.colormaps['rainbow']`` with *under* colour = white.
    density_mode : {"hist", "scatter"}
        * "hist"  – show 2-D charge density via ``imshow`` (default).
        * "scatter" – colour individual particles by charge density in their
          histogram cell (like the example we discussed).
    use_log : bool
        Use logarithmic colour normalisation.
    downsample : int | None
        If set and mode is 'scatter', use only a random subset of particles.
    seed : int | None
        Optional seed for reproducibility (used when downsampling).

    Returns
    -------
    im | PathCollection
        Handle to the created artist (image or scatter).
    """

    if q is None:
        q = np.full_like(x, 1.0 / len(x))
    else:
        q = np.asarray(q)
        if q.shape != x.shape:
            raise ValueError("q must have same shape as x and y")

    if ax is None:
        plt.figure(nfig, figsize=figsize)
        ax = plt.subplot(111)

    if title is not None:
        ax.set_title(title)

    # colour map setup
    if cmap == "my_rainbow":
        my_rainbow = matplotlib.colormaps['rainbow'].copy()
        my_rainbow.set_under('w')
    else:
        my_rainbow = matplotlib.colormaps[cmap]

    # limits
    if limits is None:
        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        dy = y_max - y_min
        limits = [[x_min, x_max], [y_min - 0.05 * dy, y_max + 0.05 * dy]]

    # histogram of total charge per bin
    H, xedges, yedges = np.histogram2d(x, y, bins=(nbins_x, nbins_y),
                                       range=limits, weights=q)

    if density_mode == "hist":
        H = H.T  # transpose so y is vertical
        vmin = np.min(H[np.nonzero(H)]) if use_log else None
        norm = LogNorm(vmin=vmin, vmax=H.max()) if use_log else None
        im = ax.imshow(H, interpolation=interpolation, origin='lower',
                       aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       cmap=my_rainbow, norm=norm)
    elif density_mode == "scatter":
        # downsample if requested
        if downsample is not None and downsample < len(x):
            rng = np.random.default_rng(seed)
            idx = rng.choice(len(x), size=downsample, replace=False)
            x, y, q = x[idx], y[idx], q[idx]

        # find bin index for each particle
        ix = np.clip(np.digitize(x, xedges) - 1, 0, nbins_x - 1)
        iy = np.clip(np.digitize(y, yedges) - 1, 0, nbins_y - 1)
        density = H[ix, iy]
        vmin = np.min(density[density > 0]) if use_log else None
        norm = LogNorm(vmin=vmin, vmax=density.max()) if use_log else None
        im = ax.scatter(x, y, c=density, s=0.3, cmap=my_rainbow,
                         norm=norm, alpha=0.5)
    else:
        raise ValueError("mode must be 'hist' or 'scatter'")

    # labels etc.
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.grid(grid)
    plt.setp(ax.get_xticklabels(), visible=show_xtick_label)

    # colour-bar switched off
    if ax.figure.get_axes()[0] is ax and 0:
        cbar = ax.figure.colorbar(im, ax=ax)
        if density_mode == "hist":
            cbar.set_label(r"Charge per bin (arb. units)")
        else:
            cbar.set_label(r"Local charge density (arb. units)")

    return im


def show_e_beam(p_array, nparts_in_slice=5000, smooth_param=0.05, nbins_x=200, nbins_y=200,
                interpolation="bilinear", inverse_tau=False,
                show_moments=False, figname=40,
                title=None, figsize=None, grid=True,
                filename=None, headtail=True,
                filter_base=2, filter_iter=2,
                tau_units="mm", p_units=None, cmap="my_rainbow", **kwargs
                ):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param nbins_x: number of bins for 2D hist. in horz. plane
    :param nbins_y: number of bins for 2D hist. in vertical plane
    :param interpolation: "bilinear", and acceptable values are 'none’, ‘nearest’, ‘bilinear’, ‘bicubic’, ‘spline16’,
                        ‘spline36’, ‘hanning’, ‘hamming’, ‘hermite’, ‘kaiser’, ‘quadric’, ‘catrom’, ‘gaussian’, ‘bessel’
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param show_moments: False, show moments (X_mean_slice and Y_mean_slice) in the density distribution
    :param figname: number of the figure or name
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param grid: True, show grid
    :param filename: None or str,  filename to save picture in the file
    :param headtail: True, shows where is the beam head is.
    :param filter_base: support of rectangle filter is 2*p+1
    :param filter_iter: the number of the filter iterations
    :param tau_units: unints of longitudinal coordinate tau - 'm', 'mm', 'um', "ps", "fs
    :param cmap: color map
    :return:
    """
    allowed = {'downsample', 'density_mode'}
    filtered = {k: v for k, v in kwargs.items() if k in allowed}
    if tau_units == "m":
        tau_factor = 1  # m
        tau_label = r"$s\,[\mathrm{m}]$"
    elif tau_units == "um":
        tau_factor = 1e6  # um
        tau_label = r"$s\,[\mathrm{\mu{}m}]$"
    elif tau_units == "ps":
        tau_factor = 1e12/speed_of_light  # pm
        tau_label = r"$t\,[\mathrm{ps}]$"
    elif tau_units == "fs":
        tau_factor = 1e15/speed_of_light  # fm
        tau_label = r"$t\,[\mathrm{fs}]$"
    else:
        tau_factor = 1e3  # mm
        tau_label = r"$s\,[\mathrm{mm}]$"

    if p_units == "MeV":
        pref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
        #Enew = p_array.p()[0] * pref + p_array.E
        p_factor = pref*1e3
        p_label = r'$E\,[MeV]$'
    else:
        p_factor = 1e2
        p_label = r'$\delta_E\,[\%]$'


    p_array_copy = deepcopy(p_array)
    if inverse_tau:
        p_array_copy.tau()[:] *= -1
    slice_params = global_slice_analysis(p_array_copy, nparts_in_slice, smooth_param, filter_base, filter_iter)
    nfig = kwargs.get('nfig', 40)
    if nfig != 40:
        figname = nfig
    fig = plt.figure(figname, figsize=figsize)
    if title is not None:
        fig.suptitle(title)
    ax_sp = plt.subplot(325)
    plt.title("Energy spread")
    plt.plot(slice_params.s * tau_factor, slice_params.se * 1e-3, "b")
    # plt.legend()
    plt.xlabel(tau_label)
    plt.ylabel(r"$\sigma_E\,[\mathrm{keV}]$")
    plt.grid(grid)

    ax_em = plt.subplot(323, sharex=ax_sp)
    plt.title("Emittances")
    emitxn_mm_mrad = np.round(slice_params.emitxn * 1e6, 2)
    emityn_mm_mrad = np.round(slice_params.emityn * 1e6, 2)
    plt.plot(slice_params.s * tau_factor, slice_params.exn*1e6, "r",
             label=fr"$\varepsilon_x^{{\mathrm{{proj}}}} = {emitxn_mm_mrad}\,\mathrm{{mm\cdot{{}}mrad}}$")
    plt.plot(slice_params.s * tau_factor, slice_params.eyn*1e6, "b",
             label=rf"$\varepsilon_y^{{\mathrm{{proj}}}} = {emityn_mm_mrad}\,\mathrm{{mm\cdot{{}}mrad}}$")
    plt.legend()
    plt.setp(ax_em.get_xticklabels(), visible=False)
    plt.ylabel(r"$\varepsilon_{x,y}\,[\mathrm{mm\cdot{}mrad}]$")
    plt.grid(grid)

    ax_c = plt.subplot(321, sharex=ax_sp)
    plt.title("Current")

    plt.plot(slice_params.s * tau_factor, slice_params.I, "b")
    imax = np.max(slice_params.I)
    imax_label = rf"$I_{{\mathrm{{max}}}}= {imax:.0f}\,\mathrm{{A}}$"
    leg = ax_c.legend([imax_label], handlelength=0, handletextpad=0, fancybox=True, loc="best")
    for item in leg.legend_handles:
        item.set_visible(False)


    plt.setp(ax_c.get_xticklabels(), visible=False)
    plt.ylabel("$I\,[\mathrm{A}]$")
    plt.grid(grid)

    ax_ys = plt.subplot(326, sharex=ax_sp)

    show_density(p_array_copy.tau() * tau_factor, p_array_copy.y() * 1e3, ax=ax_ys, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, xlabel=tau_label, ylabel='$y\,[\mathrm{mm}]$', nfig=50,
                 title="Side view", figsize=None, grid=grid, cmap=cmap, **filtered)
    if show_moments:
        plt.plot(slice_params.s * tau_factor, slice_params.my * 1e3, "k", lw=2)
        plt.plot(slice_params.s * tau_factor, (slice_params.my + slice_params.sig_y) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * tau_factor, (slice_params.my - slice_params.sig_y) * 1e3, "w", lw=1)

    ax_xs = plt.subplot(324, sharex=ax_sp)

    show_density(p_array_copy.tau() * tau_factor, p_array_copy.x() * 1e3, ax=ax_xs, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, ylabel=r'$x\,[\mathrm{mm}]$',
                 title="Top view", grid=grid, show_xtick_label=False, cmap=cmap, **filtered)
    if show_moments:
        plt.plot(slice_params.s * tau_factor, slice_params.mx * 1e3, "k", lw=2)
        plt.plot(slice_params.s * tau_factor, (slice_params.mx + slice_params.sig_x) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * tau_factor, (slice_params.mx - slice_params.sig_x) * 1e3, "w", lw=1)

    ax_ps = plt.subplot(322, sharex=ax_sp)

    show_density(p_array_copy.tau() * tau_factor, p_array_copy.p() * p_factor, ax=ax_ps, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, ylabel=p_label,
                 title="Longitudinal phase space", grid=grid, show_xtick_label=False, cmap=cmap, **filtered)

    if inverse_tau:
        arrow = r"$\Longrightarrow$"
        label = "Head " + arrow
        location = "upper right"
    else:
        arrow = r"$\Longleftarrow$"
        label = arrow + " Head"
        location = "upper left"

    if headtail:
        # Use previous legend's properties to set this AnchoredText instance to look
        # just like a legend (to match the style).
        frame = leg.get_frame()
        anchored_text = AnchoredText(
            label,
            loc=location,
            prop=dict(fontsize=leg.get_texts()[0].get_fontsize()))
        anchored_text.patch.set_boxstyle(frame.get_boxstyle())
        anchored_text.patch.set_alpha(frame.get_alpha())
        anchored_text.patch.set_edgecolor(leg.get_frame().get_edgecolor())
        ax_ps.add_artist(anchored_text)

    if filename is not None:
        plt.savefig(filename)


def show_phase_space(p_array, nparts_in_slice=5000, smooth_param=0.05, nbins_x=200, nbins_y=200,
                     interpolation="bilinear", inverse_tau=False,
                     show_moments=False, nfig=40, title=None, figsize=None, grid=True):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param nbins_x: number of bins for 2D hist. in horz. plane
    :param nbins_y: number of bins for 2D hist. in vertical plane
    :param interpolation: "bilinear", and acceptable values are 'none’, ‘nearest’, ‘bilinear’, ‘bicubic’, ‘spline16’,
                        ‘spline36’, ‘hanning’, ‘hamming’, ‘hermite’, ‘kaiser’, ‘quadric’, ‘catrom’, ‘gaussian’, ‘bessel’
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param show_moments: False, show moments (X_mean_slice and Y_mean_slice) in the density distribution
    :param nfig: number of the figure
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param grid: True, show grid
    :return:
    """
    p_array_copy = deepcopy(p_array)
    if inverse_tau:
        p_array_copy.tau()[:] *= -1
    slice_params = global_slice_analysis(p_array_copy, nparts_in_slice, smooth_param, 2, 2)

    fig = plt.figure(nfig, figsize=figsize)
    if title != None:
        fig.suptitle(title)
    ax_sp = plt.subplot(224)
    show_density(p_array_copy.tau() * 1e3, p_array_copy.py() * 1e3, ax=ax_sp, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, xlabel='s [mm]', ylabel=r'$p_y$ [mrad]', nfig=50,
                 title="", figsize=None, grid=grid, show_xtick_label=True)
    if show_moments:
        plt.plot(slice_params.s * 1e3, slice_params.myp * 1e3, "k", lw=2)
        plt.plot(slice_params.s * 1e3, (slice_params.myp + slice_params.sig_yp) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * 1e3, (slice_params.myp - slice_params.sig_yp) * 1e3, "w", lw=1)

    ax_ys = plt.subplot(222, sharex=ax_sp)

    show_density(p_array_copy.tau() * 1e3, p_array_copy.y() * 1e3, ax=ax_ys, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, ylabel='y [mm]', nfig=50,
                 title="Side view", figsize=None, grid=grid, show_xtick_label=False)
    if show_moments:
        plt.plot(slice_params.s * 1e3, slice_params.my * 1e3, "k", lw=2)
        plt.plot(slice_params.s * 1e3, (slice_params.my + slice_params.sig_y) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * 1e3, (slice_params.my - slice_params.sig_y) * 1e3, "w", lw=1)

    ax_xs = plt.subplot(221, sharex=ax_sp)

    show_density(p_array_copy.tau() * 1e3, p_array_copy.x() * 1e3, ax=ax_xs, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, ylabel='x [mm]',
                 title="Top view", grid=grid, show_xtick_label=False)
    if show_moments:
        plt.plot(slice_params.s * 1e3, slice_params.mx * 1e3, "k", lw=2)
        plt.plot(slice_params.s * 1e3, (slice_params.mx + slice_params.sig_x) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * 1e3, (slice_params.mx - slice_params.sig_x) * 1e3, "w", lw=1)

    ax_xp = plt.subplot(223, sharex=ax_sp)

    show_density(p_array_copy.tau() * 1e3, p_array_copy.px() * 1e3, ax=ax_xp, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, xlabel='s [mm]', ylabel=r'$p_x$ [mrad]',
                 title="", grid=grid, show_xtick_label=True)
    if show_moments:
        plt.plot(slice_params.s * 1e3, slice_params.mxp * 1e3, "k", lw=2)
        plt.plot(slice_params.s * 1e3, (slice_params.mxp + slice_params.sig_xp) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * 1e3, (slice_params.mxp - slice_params.sig_xp) * 1e3, "w", lw=1)

    return slice_params.s, slice_params.myp


def compare_beams(p_array_1, p_array_2, nparts_in_slice1=5000, nparts_in_slice2=5000, smoth_param=0.05,
                  inverse_tau=False, nfig=40, title=None, figsize=None, legend_beam1=None, legend_beam2=None):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array_1: ParticleArray
    :param p_array_2: ParticleArray
    :param nparts_in_slice1: number of particles per slice in p_array_1
    :param nparts_in_slice2: number of particles per slice in p_array_2
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param nfig: number of the figure
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param legend_beam1: None, legend for beam N1
    :param legend_beam2: None, legend for beam N1
    :return:
    """
    tau_label = r"$s\,[\mathrm{\mu{}m}]$"
    if legend_beam1 == None:
        legend_beam1 = "beam 1"
    if legend_beam2 == None:
        legend_beam2 = "beam 2"
    p_array_copy1 = deepcopy(p_array_1)
    p_array_copy2 = deepcopy(p_array_2)
    if inverse_tau:
        p_array_copy1.tau()[:] *= -1
        p_array_copy2.tau()[:] *= -1

    slice_params1 = global_slice_analysis(p_array_copy1, nparts_in_slice1, smoth_param, 2, 2)
    slice_params2 = global_slice_analysis(p_array_copy2, nparts_in_slice2, smoth_param, 2, 2)

    fig = plt.figure(nfig, figsize=figsize)
    if title != None:
        fig.suptitle(title)
    ax_sp = plt.subplot(325)
    plt.title("Energy Spread")
    plt.plot(slice_params1.s * 1e6, slice_params1.se * 1e-3, "r",
             label=rf"$\sigma_E$: {legend_beam1}")
    plt.plot(slice_params2.s * 1e6, slice_params2.se * 1e-3, "b",
             label=rf"$\sigma_E$: {legend_beam2}")
    plt.legend()
    plt.xlabel(tau_label)

    plt.ylabel(r"$\Delta{} E\,[\mathrm{keV}]$")
    plt.grid(True)

    ax_em = plt.subplot(323, sharex=ax_sp)
    plt.title("Horizontal Emittances")
    emitxn1_mm_mrad = np.round(slice_params1.emitxn * 1e6, 2)
    emitxn2_mm_mrad = np.round(slice_params2.emitxn * 1e6, 2)
    plt.plot(slice_params1.s * 1e6, slice_params1.ex, "r",
             label=rf"$\varepsilon_x = {emitxn1_mm_mrad}\,\mathrm{{mm\cdot{{}}mrad}}$: {legend_beam1}")
    plt.plot(slice_params2.s * 1e6, slice_params2.ex, "b",
             label=rf"$\varepsilon_x = {emitxn2_mm_mrad}\,\mathrm{{mm\cdot{{}}mrad}}$: {legend_beam2}")

    plt.legend()
    plt.setp(ax_em.get_xticklabels(), visible=False)
    plt.ylabel(r"$\varepsilon_x\,[\mathrm{mm\cdot{}mrad}]$")
    plt.grid(True)

    ax_c = plt.subplot(321, sharex=ax_sp)
    plt.title("Current")
    i1max = np.max(slice_params1.I)
    i2max = np.max(slice_params2.I)
    plt.plot(slice_params1.s * 1e6, slice_params1.I, "r",
             label=fr"$I_{{\mathrm{{max}}}}={i1max:.0f}\,\mathrm{{A}}$ {legend_beam1}")
    plt.plot(slice_params2.s * 1e6, slice_params2.I, "b",
             label=fr"$I_{{\mathrm{{max}}}}={i2max:.0f}\,\mathrm{{A}}$ {legend_beam2}")
    plt.legend()
    plt.setp(ax_c.get_xticklabels(), visible=False)
    plt.ylabel("$I\,[\mathrm{A}]$")
    plt.grid(True)

    # my_rainbow = deepcopy(plt.get_cmap('rainbow'))
    # my_rainbow.set_under('w')

    ax_mp = plt.subplot(326)
    plt.title("Energy Distribution")

    plt.plot(slice_params1.s * 1e6, slice_params1.mp * 1e2, "r", label=r"$\Delta{} E/E$: " + legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.mp * 1e2, "b", label=r"$\Delta{} E/E$: " + legend_beam2)
    plt.legend()
    plt.xlabel(tau_label)
    plt.ylabel(r'$\delta{}_E\,[\%]$')
    plt.grid(True)

    ax_em = plt.subplot(324, sharex=ax_mp)
    plt.title("Vertical Emittances")
    emityn1_mm_mrad = np.round(slice_params1.emityn * 1e6, 2)
    emityn2_mm_mrad = np.round(slice_params2.emityn * 1e6, 2)
    plt.plot(slice_params1.s * 1e6, slice_params1.ey, "r",
             label=r"$\varepsilon_y = $" + str(emityn1_mm_mrad) + r" $\mathrm{mm\cdot{}mrad}$: " + legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.ey, "b",
             label=r"$\varepsilon_y = $" + str(emityn2_mm_mrad) + r" $\mathrm{mm\cdot{}mrad}$: " + legend_beam2)
    plt.legend()
    plt.setp(ax_em.get_xticklabels(), visible=False)
    plt.ylabel(r"$\varepsilon_y\,[\mathrm{mm\cdot{}mrad}]$")
    plt.grid(True)

    ax_e = plt.subplot(322, sharex=ax_mp)
    plt.title("Slice Positions")
    plt.plot(slice_params1.s * 1e6, slice_params1.mx * 1e6, "r", label=r"$x_\mathrm{slice}$: " + legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.mx * 1e6, "b", label=r"$x_\mathrm{slice}$: " + legend_beam2)
    plt.plot(slice_params1.s * 1e6, slice_params1.my * 1e6, "r--", label=r"$y_\mathrm{slice}$: " + legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.my * 1e6, "b--", label=r"$y_\mathrm{slice}$: " + legend_beam2)
    plt.legend()
    plt.setp(ax_e.get_xticklabels(), visible=False)
    plt.ylabel(r"$x_{\mathrm{slice}},\,y_{\mathrm{slice}}\,[\mathrm{\mu{}m}]$")
    plt.grid(True)


def compare_beams_reduced(p_array_1, p_array_2, nparts_in_slice=5000, smoth_param=0.05,
                          inverse_tau=True, nfig=40, title=None, figsize=None, legend_beam1=None, legend_beam2=None):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array_1: ParticleArray
    :param p_array_2: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param nfig: number of the figure
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param legend_beam1: None, legend for beam N1
    :param legend_beam2: None, legend for beam N1
    :return:
    """
    if legend_beam1 == None:
        legend_beam1 = "beam1"
    if legend_beam2 == None:
        legend_beam2 = "beam2"
    p_array_copy1 = deepcopy(p_array_1)
    p_array_copy2 = deepcopy(p_array_2)
    if inverse_tau:
        p_array_copy1.tau()[:] *= -1
        p_array_copy2.tau()[:] *= -1

    slice_params1 = global_slice_analysis(p_array_copy1, nparts_in_slice, smoth_param, 2, 2)
    slice_params2 = global_slice_analysis(p_array_copy2, nparts_in_slice, smoth_param, 2, 2)

    fig = plt.figure(nfig, figsize=figsize)
    if title != None:
        fig.suptitle(title)

    ax_sp = plt.subplot(223)
    plt.title("Slice energy spread")
    plt.plot(slice_params1.s * 1e6, slice_params1.se * 1e-3, "r", label=legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.se * 1e-3, "b", label=legend_beam2)
    # plt.legend()
    plt.xlabel(r"s [$\mu m$]")
    plt.ylabel(r"$\sigma_E$ [keV]")

    ax_mp = plt.subplot(224)
    plt.title("Mean slice energy")

    plt.plot(slice_params1.s * 1e6, slice_params1.mp * 1e2, "r", label=legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.mp * 1e2, "b", label=legend_beam2)
    # plt.legend()
    plt.xlabel(r"s [$\mu m$]")
    plt.ylabel(r"$\delta_E$ [%]")

    ax_em = plt.subplot(222, sharex=ax_mp)
    plt.title("Horizontal slice emittance")
    emitxn1_mm_mrad = np.round(slice_params1.emitxn * 1e6, 2)
    emitxn2_mm_mrad = np.round(slice_params2.emitxn * 1e6, 2)
    plt.plot(slice_params1.s * 1e6, slice_params1.ex, "r",
             label=legend_beam1 + ": " + r"$\varepsilon_x^{proj} = $" + str(emitxn1_mm_mrad))
    plt.plot(slice_params2.s * 1e6, slice_params2.ex, "b",
             label=legend_beam2 + ": " + r"$\varepsilon_x^{proj} = $" + str(emitxn2_mm_mrad))
    plt.legend()
    plt.setp(ax_em.get_xticklabels(), visible=False)
    plt.ylabel(r"$\varepsilon_x$ [$\mu m$]")

    ax_c = plt.subplot(221, sharex=ax_sp)
    plt.title("Current")
    plt.plot(slice_params1.s * 1e6, slice_params1.I, "r", label=legend_beam1)
    plt.plot(slice_params2.s * 1e6, slice_params2.I, "b", label=legend_beam2)
    # plt.legend()
    plt.setp(ax_c.get_xticklabels(), visible=False)
    plt.ylabel("I [A]")


def show_e_beam_reduced(p_array, nparts_in_slice=5000, smooth_param=0.05, nbins_x=200, nbins_y=200,
                        interpolation="bilinear", inverse_tau=False,
                        show_moments=False, nfig=40, title=None, figsize=None, grid=True, filename=None):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param nbins_x: number of bins for 2D hist. in horz. plane
    :param nbins_y: number of bins for 2D hist. in vertical plane
    :param interpolation: "bilinear", and acceptable values are 'none’, ‘nearest’, ‘bilinear’, ‘bicubic’, ‘spline16’,
                        ‘spline36’, ‘hanning’, ‘hamming’, ‘hermite’, ‘kaiser’, ‘quadric’, ‘catrom’, ‘gaussian’, ‘bessel’
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param show_moments: False, show moments (X_mean_slice and Y_mean_slice) in the density distribution
    :param nfig: number of the figure
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param grid: True, show grid
    :param filename: None or str,  filename to save picture in the file
    :return:
    """
    p_array_copy = deepcopy(p_array)
    if inverse_tau:
        p_array_copy.tau()[:] *= -1
    slice_params = global_slice_analysis(p_array_copy, nparts_in_slice, smooth_param, 2, 2)

    fig = plt.figure(nfig, figsize=figsize)
    if title != None:
        fig.suptitle(title)

    # ax_sp = plt.subplot(225)
    # plt.title("Energy spread")
    # plt.plot(slice_params.s * 1e3, slice_params.se*1e-3, "b")
    # plt.legend()
    # plt.xlabel("s [mm]")
    # plt.ylabel("$\sigma_E$ [keV]")
    # plt.grid(grid)

    ax_sp = plt.subplot(223)
    plt.title("Emittances")
    emitxn_mm_mrad = np.round(slice_params.emitxn * 1e6, 2)
    emityn_mm_mrad = np.round(slice_params.emityn * 1e6, 2)
    plt.plot(slice_params.s * 1e3, slice_params.ex, "r", label=r"$\varepsilon_x^{proj} = $" + str(emitxn_mm_mrad))
    plt.plot(slice_params.s * 1e3, slice_params.ey, "b", label=r"$\varepsilon_y^{proj} = $" + str(emityn_mm_mrad))
    plt.legend()
    # plt.setp(ax_sp.get_xticklabels(), visible=False)
    plt.xlabel("s [mm]")
    plt.ylabel(r"$\varepsilon_{x,y}$ [$\mu m$]")
    plt.grid(grid)

    ax_c = plt.subplot(221, sharex=ax_sp)
    plt.title("Current")
    plt.plot(slice_params.s * 1e3, slice_params.I, "b", label=r"$I_{max}=$" + str(np.round(np.max(slice_params.I), 1)))
    # plt.legend()
    plt.setp(ax_c.get_xticklabels(), visible=False)
    plt.ylabel("I [A]")
    plt.grid(grid)

    ax_ys = plt.subplot(224, sharex=ax_sp)

    show_density(p_array_copy.tau() * 1e3, p_array_copy.y() * 1e3, ax=ax_ys, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, xlabel='s [mm]', ylabel='y [mm]', nfig=50,
                 title="Side view", figsize=None, grid=grid)
    if show_moments:
        plt.plot(slice_params.s * 1e3, slice_params.my * 1e3, "k", lw=2)
        plt.plot(slice_params.s * 1e3, (slice_params.my + slice_params.sig_y) * 1e3, "w", lw=1)
        plt.plot(slice_params.s * 1e3, (slice_params.my - slice_params.sig_y) * 1e3, "w", lw=1)

    # ax_xs = plt.subplot(324, sharex=ax_sp)
    #
    # show_density(p_array_copy.tau() * 1e3, p_array_copy.x() * 1e3, ax=ax_xs, nbins_x=nbins_x, nbins_y=nbins_y,
    #             interpolation=interpolation, ylabel='x [mm]',
    #             title="Top view", grid=grid, show_xtick_label=False)
    # if show_moments:
    #    plt.plot(slice_params.s * 1e3, slice_params.mx * 1e3, "k", lw=2)
    #    plt.plot(slice_params.s * 1e3, (slice_params.mx + slice_params.sig_x) * 1e3, "w", lw=1)
    #    plt.plot(slice_params.s * 1e3, (slice_params.mx - slice_params.sig_x) * 1e3, "w", lw=1)
    #
    ax_ps = plt.subplot(222, sharex=ax_sp)

    show_density(p_array_copy.tau() * 1e3, p_array_copy.p() * 1e2, ax=ax_ps, nbins_x=nbins_x, nbins_y=nbins_y,
                 interpolation=interpolation, ylabel=r'$\delta_E$ [%]',
                 title="Longitudinal phase space", grid=grid, show_xtick_label=False)

    if filename is not None:
        plt.savefig(filename)


def show_e_beam_slices(p_array, nparts_in_slice=5000, smooth_param=0.05, inverse_tau=False, figname=50,
                       title=None, figsize=None, grid=True,
                       filename=None, headtail=True,
                       filter_base=2, filter_iter=2, tau_units="mm", slice=0):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param show_moments: False, show moments (X_mean_slice and Y_mean_slice) in the density distribution
    :param figname: number of the figure or name
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param grid: True, show grid
    :param headtail: True, show direction of the bunch
    :param tau_units: [m], [mm] or [um]
    :param slice: None, str or float:
                                None - ignore
                                float - s position of the slice
                                "center" - center of mass of the bunch
                                "Imax" - current maximum
                                "Emax" - maximum energy
    :return:
    """

    if tau_units == "m":
        tau_factor = 1  # m
        tau_label = "s [mm]"
    elif tau_units == "um":
        tau_factor = 1e6  # um
        tau_label = r"s [$\mu$m]"
    else:
        tau_factor = 1e3  # mm
        tau_label = "s [mm]"

    p_array_copy = deepcopy(p_array)
    if inverse_tau:
        p_array_copy.tau()[:] *= -1
    slice_params = global_slice_analysis(p_array_copy, nparts_in_slice, smooth_param, filter_base, filter_iter)

    s = slice_params.s * tau_factor
    if isinstance(slice, numbers.Number):
        ind0 = np.argsort(np.abs(s - slice))[0]
    elif isinstance(slice, str):
        if slice == "center":
            ind0 = np.argsort(np.abs(s - np.mean(s)))[0]
        elif slice == "Imax":
            ind0 = np.argmax(slice_params.I)
        elif slice == "Emax":
            ind0 = np.argmax(slice_params.me)
    else:
        ind0 = np.argsort(np.abs(slice_params.s))[0]

    fig = plt.figure(figname, figsize=figsize)
    if title is not None:
        fig.suptitle(title)
    ax_sp = plt.subplot(325)
    plt.title("Energy spread")
    plt.plot(s, slice_params.se * 1e-3, "b")
    fm = "bo"
    label = r"$\sigma_E=$" + str(np.round(slice_params.se[ind0] * 1e-3, 2)) + " keV"
    if slice is None:
        fm = "b"
        label = None
    plt.plot(s[ind0], slice_params.se[ind0] * 1e-3, fm, label=label)
    if slice is not None:
        plt.legend()
    plt.xlabel(tau_label)
    plt.ylabel(r"$\sigma_E$ [keV]")
    #plt.ylim([200, 300])
    plt.grid(grid)

    ax_em = plt.subplot(323, sharex=ax_sp)

    plt.title("Emittances")

    plt.plot(s, slice_params.ex, "r")
    fmx = "ro"
    fmy = "bo"
    labelx = r"$\varepsilon_x^{0} = $" + str(np.round(slice_params.ex[ind0], 2)) + r" $\mu m \cdot rad$"
    labely = r"$\varepsilon_y^{0} = $" + str(np.round(slice_params.ey[ind0], 2)) + r" $\mu m \cdot rad$"
    if slice is None:
        fmx = "r"
        fmy = "b"
        labelx = r"$\varepsilon_x$"
        labely = r"$\varepsilon_y$"
    plt.plot(s[ind0], slice_params.ex[ind0], fmx, label=labelx)
    plt.plot(s, slice_params.ey, "b")
    plt.plot(s[ind0], slice_params.ey[ind0], fmy, label=labely)
    plt.legend()
    plt.setp(ax_em.get_xticklabels(), visible=False)
    plt.ylabel(r"$\varepsilon_{x,y}$ [$\mu m \cdot rad$]")
    plt.grid(grid)

    ax_c = plt.subplot(321, sharex=ax_sp)
    plt.title("Current")

    if inverse_tau:
        arrow = r"$\Longrightarrow$"
        label_arr = "head " + arrow
        location = "upper right"
    else:
        arrow = r"$\Longleftarrow$"
        label_arr = arrow + " head"
        location = "upper left"

    plt.plot(s, slice_params.I, "b")
    fm = "bo"
    label = "$I_0=$" + str(np.round(slice_params.I[ind0], 1)) + " A"
    if slice is None:
        fm = "b"
        label = None
    plt.plot(s[ind0], slice_params.I[ind0], fm, label=label)
    # label = r"$I_{max}=$" + str(np.round(np.max(slice_params.I), 1))
    if slice is not None:
        plt.legend()
    if headtail:
        leg = ax_c.legend([label_arr], handlelength=0, handletextpad=0, fancybox=True, loc=location)
        for item in leg.legend_handles:
            item.set_visible(False)
    plt.setp(ax_c.get_xticklabels(), visible=False)
    plt.ylabel("I [A]")
    plt.grid(grid)

    ax_alpha = plt.subplot(326, sharex=ax_sp)
    plt.title(r"$\alpha_{x,y}$")
    plt.plot(s, slice_params.alpha_x, "r")
    fmx = "ro"
    fmy = "bo"
    labelx = r"$\alpha_x=$" + str(np.round(slice_params.alpha_x[ind0], 2))
    labely = r"$\alpha_y=$" + str(np.round(slice_params.alpha_y[ind0], 2))
    if slice is None:
        fmx = "r"
        fmy = "b"
        labelx = r"$\alpha_x$"
        labely = r"$\alpha_y$"
    plt.plot(s[ind0], slice_params.alpha_x[ind0], fmx, label=labelx)
    plt.plot(s, slice_params.alpha_y, "b")
    plt.plot(s[ind0], slice_params.alpha_y[ind0], fmy, label=labely)
    plt.legend()
    plt.xlabel(tau_label)
    plt.ylabel(r"$\alpha_{x,y}$")
    plt.grid(grid)

    ax_b = plt.subplot(324, sharex=ax_sp)
    plt.title(r"$\beta_{x,y}$")
    plt.plot(s, slice_params.beta_x, "r")
    fmx = "ro"
    fmy = "bo"
    labelx = r"$\beta_x=$" + str(np.round(slice_params.beta_x[ind0], 2)) + " m"
    labely = r"$\beta_y=$" + str(np.round(slice_params.beta_y[ind0], 2)) + " m"
    if slice is None:
        fmx = "r"
        fmy = "b"
        labelx = r"$\beta_x$"
        labely = r"$\beta_y$"
    plt.plot(s[ind0], slice_params.beta_x[ind0], fmx, label=labelx)
    plt.plot(s, slice_params.beta_y, "b")
    plt.plot(s[ind0], slice_params.beta_y[ind0], fmy, label=labely)
    plt.setp(ax_b.get_xticklabels(), visible=False)
    plt.ylabel(r"$\beta_{x,y}$")
    plt.grid(grid)
    plt.legend()

    ax_me = plt.subplot(322, sharex=ax_sp)
    plt.title(r"Longitudinal phase space")
    MeV = 1e-6
    plt.plot(s, slice_params.me * MeV, "r")
    fm = "ro"
    label = r"$E_0=$" + str(np.round(slice_params.me[ind0] * MeV, 1)) + " MeV"
    if slice is None:
        fm = "r"
        label = None
    plt.plot(s[ind0], slice_params.me[ind0] * MeV, fm, label=label)
    # plt.plot(slice_params.s[ind0] * tau_factor, bx, "ro", label=r"$\beta_x=$"+str(np.round(bx, 2)))
    plt.setp(ax_me.get_xticklabels(), visible=False)
    plt.ylabel(r"$E$ [MeV]")
    #plt.ylim([17100, 17300])
    plt.grid(grid)
    if slice is not None:
        plt.legend()
    if filename is not None:
        plt.savefig(filename)


def beam_jointplot(p_array, show_plane="x", nparts_in_slice=5000, smooth_param=0.05, nbins_x=200, nbins_y=200,
                   interpolation="bilinear", inverse_tau=True, show_head=True,
                   show_moments=False, nfig=40, title=None, figsize=None, grid=True, filename=None):
    """
    Shows e-beam slice parameters (current, emittances, energy spread)
    and beam distributions (dE/(p0 c), X, Y) against long. coordinate (S)
    Note: beam head is on the left side

    :param p_array: ParticleArray
    :param nparts_in_slice: number of particles per slice
    :param smoth_param: 0.05, smoothing parameters to calculate the beam current: sigma = smoth_param * np.std(p_array.tau())
    :param nbins_x: number of bins for 2D hist. in horz. plane
    :param nbins_y: number of bins for 2D hist. in vertical plane
    :param interpolation: "bilinear", and acceptable values are 'none’, ‘nearest’, ‘bilinear’, ‘bicubic’, ‘spline16’,
                        ‘spline36’, ‘hanning’, ‘hamming’, ‘hermite’, ‘kaiser’, ‘quadric’, ‘catrom’, ‘gaussian’, ‘bessel’
    :param inverse_tau: False, inverse tau - head will be on the right side of figure
    :param show_moments: False, show moments (X_mean_slice and Y_mean_slice) in the density distribution
    :param nfig: number of the figure
    :param title: None or string - title of the figure
    :param figsize: None or e.g. (8, 6)
    :param grid: True, show grid
    :param filename: None or str,  filename to save picture in the file
    :return:
    """
    p_array_copy = deepcopy(p_array)
    if inverse_tau:
        p_array_copy.tau()[:] *= -1
    slice_params = global_slice_analysis(p_array_copy, nparts_in_slice, smooth_param, 2, 2)

    fig = plt.figure(nfig, figsize=figsize)
    if title is not None:
        fig.suptitle(title)

    ax_top = plt.subplot(211)
    plt.title("Current")

    if inverse_tau:
        arrow = r"$\Longrightarrow$"
        label = "head " + arrow
        location = "upper right"
    else:
        arrow = r"$\Longleftarrow$"
        label = arrow + " head"
        location = "upper left"

    plt.plot(slice_params.s * 1e3, slice_params.I, "b")
    # label = r"$I_{max}=$" + str(np.round(np.max(slice_params.I), 1))
    if show_head:
        leg = ax_top.legend([label], handlelength=0, handletextpad=0, fancybox=True, loc=location)
        for item in leg.legend_handles:
            item.set_visible(False)
    plt.setp(ax_top.get_xticklabels(), visible=False)
    plt.ylabel("I [A]")
    plt.grid(grid)
    if show_plane == "y":
        ax_down = plt.subplot(212, sharex=ax_top)

        show_density(p_array_copy.tau() * 1e3, p_array_copy.y() * 1e3, ax=ax_down, nbins_x=nbins_x, nbins_y=nbins_y,
                     interpolation=interpolation, xlabel='s [mm]', ylabel='y [mm]', nfig=50,
                     title="Side view", figsize=None, grid=grid, show_xtick_label=True)
        if show_moments:
            plt.plot(slice_params.s * 1e3, slice_params.my * 1e3, "k", lw=2)
            plt.plot(slice_params.s * 1e3, (slice_params.my + slice_params.sig_y) * 1e3, "w", lw=1)
            plt.plot(slice_params.s * 1e3, (slice_params.my - slice_params.sig_y) * 1e3, "w", lw=1)
    elif show_plane == "x":
        ax_down = plt.subplot(212, sharex=ax_top)

        show_density(p_array_copy.tau() * 1e3, p_array_copy.x() * 1e3, ax=ax_down, nbins_x=nbins_x, nbins_y=nbins_y,
                     interpolation=interpolation, ylabel='x [mm]', xlabel='s [mm]',
                     title="Top view", grid=grid, show_xtick_label=True)
        if show_moments:
            plt.plot(slice_params.s * 1e3, slice_params.mx * 1e3, "k", lw=2)
            plt.plot(slice_params.s * 1e3, (slice_params.mx + slice_params.sig_x) * 1e3, "w", lw=1)
            plt.plot(slice_params.s * 1e3, (slice_params.mx - slice_params.sig_x) * 1e3, "w", lw=1)
    else:
        ax_down = plt.subplot(212, sharex=ax_top)

        show_density(p_array_copy.tau() * 1e3, p_array_copy.p() * 1e2, ax=ax_down, nbins_x=nbins_x, nbins_y=nbins_y,
                     interpolation=interpolation, ylabel=r'$\delta_E$ [%]', xlabel='s [mm]',
                     title="Longitudinal phase space", grid=grid, show_xtick_label=True)

    if filename is not None:
        plt.savefig(filename)

    return ax_top, ax_down


class Save3DBeamDensity(PhysProc):
    def __init__(self):
        PhysProc.__init__(self)
        self.energy = None
        self.napply = 0

    def apply_3d(self, p_array, dz):
        # _logger.debug(" SaveBeam applied, dz =", dz)

        my_rainbow = copy.deepcopy(plt.get_cmap('rainbow'))
        my_rainbow.set_under('w')

        y = p_array.x()[::50] * 1e6  # 10*np.random.normal(mu, sigma, 5000)
        z = p_array.y()[::50] * 1e6  # 10*np.random.normal(mu, sigma, 5000)
        x = p_array.tau()[::50] * 1e6  # 10*np.random.normal(mu, sigma, 5000)

        xyz = np.vstack([x, y, z])
        density = stats.gaussian_kde(xyz)(xyz)

        idx = density.argsort()
        x, y, z, density = x[idx], y[idx], z[idx], density[idx]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=density, cmap=my_rainbow)
        ax.set_xlim(20, 70)
        ax.set_ylim(-1500, 1500)
        ax.set_zlim(-1500, 1500)

        dig = str(self.napply)
        name = "0" * (4 - len(dig)) + dig

        plt.savefig(name)
        self.napply += 1
        plt.clf()

    def apply(self, p_array, dz):
        nbins_x = 400
        nbins_y = 400
        interpolation = "bilinear"
        fig = plt.figure(figsize=None)
        ax_xs = plt.subplot(211)

        show_density(p_array.tau() * 1e3, p_array.x() * 1e3, ax=ax_xs, nbins_x=nbins_x, nbins_y=nbins_y,
                     interpolation=interpolation, ylabel='x [mm]',
                     title="Top view", grid=True, show_xtick_label=False, limits=[[0.025, 0.075], [-2, 2]])
        ax_ys = plt.subplot(212, sharex=ax_xs)

        show_density(p_array.tau() * 1e3, p_array.y() * 1e3, ax=ax_ys, nbins_x=nbins_x, nbins_y=nbins_y,
                     interpolation=interpolation, xlabel="s, [mm]", ylabel='y [mm]', nfig=50,
                     title="Side view", figsize=None, grid=True, show_xtick_label=True,
                     limits=[[0.025, 0.075], [-2, 2]])

        dig = str(self.napply)
        name = "0" * (4 - len(dig)) + dig

        plt.savefig(name)
        self.napply += 1
        plt.clf()
