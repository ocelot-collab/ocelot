"""Artist and coordinate checks for the optics and bunch diagnostics."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.backend_bases import MouseEvent
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from ocelot.common.globals import speed_of_light
from ocelot.cpbd.beam import ParticleArray, Twiss
from ocelot.cpbd.elements import Cavity, Drift, Hcor, Monitor, Multipole, Quadrupole, SBend, Sextupole, Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.gui.accelerator import plot_elems, plot_opt_func, show_density, show_e_beam


@pytest.fixture
def figure():
    fig = plt.figure()
    yield fig
    plt.close(fig)


def test_optics_refreshes_supplied_figure_and_keeps_shared_positions(figure):
    q = Quadrupole(l=0.2, k1=1)
    lat = MagneticLattice([q, Drift(l=0.6), q])
    samples = [Twiss(s=s, E=0.13, beta_x=2, beta_y=3) for s in (12, 13)]
    inset = figure.add_axes([0.01, 0.8, 0.03, 0.1])
    other, other_ax = plt.subplots()
    try:
        with mpl.rc_context({"axes.grid": False, "grid.color": "green"}):
            assert plot_opt_func(lat, samples, plt_figure=figure, title="Initial") is None
            plot_opt_func(lat, samples, top_plot=["Dx", "E"], plt_figure=figure)
            assert mpl.rcParams["axes.grid"] is False
            assert mpl.rcParams["grid.color"] == "green"
        assert other_ax.get_xlim() == (0, 1)
        assert figure.axes[0] is inset
        top, beta, magnets = figure.axes[1:]
        assert all(ax.get_xlim() == (12, 13) for ax in (top, beta, magnets))
        top.set_xlim(12.1, 12.5)
        assert beta.get_xlim() == magnets.get_xlim() == top.get_xlim()
        assert [mpl.colors.to_hex(line.get_color()) for line in beta.lines] == ["#0000ff", "#ff0000"]
        assert [line.get_label().split()[-1] for line in top.lines] == ["[m]", "[GeV]"]
        assert magnets.get_xlabel() == "s [m]"
        for patch, bounds in zip(magnets.patches, [(12, 12.2), (12.2, 12.8), (12.8, 13)]):
            x = patch.get_xy()[:, 0]
            np.testing.assert_allclose([x.min(), x.max()], bounds)
        assert all(text.get_text() != "Initial" for text in figure.texts)
    finally:
        plt.close(other)


@pytest.mark.parametrize("with_undulator", [False, True])
def test_multipole_scale_is_independent_of_undulator(figure, with_undulator):
    ax = figure.add_subplot()
    elements = [Drift(l=1), Multipole(kn=[0, 0.2]), Multipole(kn=[0, 0.4]), Hcor()]
    if with_undulator:
        elements.append(Undulator(lperiod=0.1, nperiods=2, Ky=1))
    lat = MagneticLattice(elements)
    plot_elems(figure, ax, lat, excld_legend=[Hcor])
    _, weak, strong, corrector = ax.patches[:4]
    assert np.ptp(weak.get_xy()[:, 1]) * 2 == pytest.approx(np.ptp(strong.get_xy()[:, 1]))
    if with_undulator:
        assert np.ptp(ax.patches[-1].get_xy()[:, 1]) > 0
    assert np.ptp(corrector.get_xy()[:, 1]) == 0
    assert ax.texts[1].xy == (1, 0)
    assert "cor" not in ax.get_legend_handles_labels()[1]


def test_element_key_and_hover_in_embedded_figure():
    fig = Figure(figsize=(6.4, 4.8))
    FigureCanvasAgg(fig)
    lat = MagneticLattice([Quadrupole(l=1, k1=1), Cavity(l=1), SBend(l=1), Monitor(), Hcor(),
                           Sextupole(l=1), Undulator(lperiod=0.1, nperiods=2, Kx=1)])
    samples = [Twiss(s=0), Twiss(s=lat.totalLen)]
    plot_opt_func(lat, samples, plt_figure=fig)
    fig.canvas.draw()
    magnets = fig.axes[-1]
    key = magnets.get_legend().get_window_extent()
    assert magnets.bbox.x0 <= key.x0 < key.x1 <= magnets.bbox.x1
    assert key.y0 > max(patch.get_window_extent().y1 for patch in magnets.patches)
    x, y = magnets.transData.transform((0.5, 0.5))
    event = MouseEvent("motion_notify_event", fig.canvas, x, y)
    fig.canvas.callbacks.process(event.name, event)
    assert magnets.texts[0].get_visible()


def test_drift_only_optics_has_no_element_key(figure):
    lat = MagneticLattice([Drift(l=1)])
    plot_opt_func(lat, [Twiss(s=0), Twiss(s=1)], top_plot=["alpha_x", "alpha_y"], plt_figure=figure)
    top, _, magnets = figure.axes
    assert "[m]" not in top.get_ylabel()
    assert magnets.get_legend() is None
    figure.canvas.draw()
    assert magnets.xaxis.label.get_window_extent().y0 >= figure.bbox.y0


def test_beam_units_orientation_and_redraw(figure, tmp_path):
    rng = np.random.default_rng(8)
    beam = ParticleArray(1000)
    beam.E = 0.13
    beam.rparticles[:] = rng.normal(size=(6, beam.n)) * np.array([1e-4, 2e-5, 1e-4, 2e-5, 1e-4, 1e-3])[:, None]
    beam.q_array[:] = 1e-12
    initial = beam.rparticles.copy()
    assert show_e_beam(beam, nparts_in_slice=50, figname=figure.number, title="Bunch") is None
    phase = figure.axes[-1]
    extent = phase.images[0].get_extent()
    np.testing.assert_allclose(extent[:2], [beam.tau().min() * 1e3, beam.tau().max() * 1e3])
    assert phase.images[0].get_cmap().name == "rainbow"
    assert len(phase.artists) == 1
    saved = tmp_path / "beam.png"
    show_e_beam(beam, nparts_in_slice=50, nfig=figure.number, inverse_tau=True, tau_units="fs",
                p_units="MeV", density_mode="scatter", headtail=False, show_moments=True, filename=saved)
    spread, emit, current, side, top, phase = figure.axes
    assert saved.stat().st_size > 0
    assert len(phase.artists) == 0
    assert phase.get_shared_x_axes().joined(phase, current)
    assert "fs" in spread.get_xlabel() and "fs" in side.get_xlabel()
    assert r"\Delta E" in phase.get_ylabel()
    assert len(emit.lines) == 2 and len(side.lines) == len(top.lines) == 3
    np.testing.assert_allclose(phase.collections[0].get_offsets()[:, 0], -beam.tau() * 1e15 / speed_of_light)
    np.testing.assert_allclose(phase.collections[0].get_offsets()[:, 1], beam.p() * beam.p0c * 1e3)
    np.testing.assert_array_equal(beam.rparticles, initial)


@pytest.mark.parametrize("mode", ["hist", "scatter"])
def test_density_integer_coordinates_and_empty_window(figure, mode):
    ax = figure.add_subplot()
    xy = np.array([0, 1, 2])
    limits = [[0, 1], [0, 1]]
    artist = show_density(xy, xy, ax=ax, nbins_x=2, nbins_y=2, limits=limits, density_mode=mode)
    if mode == "hist":
        assert artist.get_array().sum() == pytest.approx(2 / 3)
    else:
        np.testing.assert_allclose(artist.get_offsets(), [[0, 0], [1, 1]])
        np.testing.assert_allclose(artist.get_array(), 1 / 3)
    assert ax.get_xlim() == ax.get_ylim() == (0, 1)
    show_density(xy, xy, ax=ax, limits=[[5, 6], [5, 6]], density_mode=mode)
    figure.canvas.draw()
