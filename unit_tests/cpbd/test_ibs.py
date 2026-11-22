"""IBS limits, coordinate conversions, and diagnostic-window regressions."""

import numpy as np
import pytest
from scipy.special import elliprd

from ocelot.common.globals import m_e_GeV, q_e, ro_e, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.physics_proc import IBS


@pytest.fixture
def beam():
    rng = np.random.default_rng(31)
    particles = ParticleArray(n=20000)
    particles.E = 0.1
    particles.q_array[:] = 0.25e-9 / particles.n
    particles.rparticles[:] = rng.normal(size=(6, particles.n)) * np.array(
        [250e-6, 2e-5, 250e-6, 2e-5, 1e-3, 1e-5])[:, None]
    return particles


def fixed_parameters(monkeypatch, energy, xi, method="Nagaitsev"):
    """Supply known round-beam moments; sample kicks with exactly unit variance."""
    particles = ParticleArray(n=2)
    particles.E = energy
    particles.q_array[:] = 0.25e-9 / particles.n
    gamma = energy / m_e_GeV
    beta = np.sqrt(1 - gamma**-2)
    theta = 2e-5
    sigma = 250e-6
    sigma_delta = np.sqrt(xi) * gamma * theta
    proc = IBS(method=method, update_Clog=False)
    proc.sigma_xy = sigma
    proc.sigma_z = 1e-3
    proc.emit_n = beta * gamma * sigma * theta
    proc.sigma_dgamma0 = beta**2 * gamma * sigma_delta
    monkeypatch.setattr(proc, "get_beam_params", lambda p: None)
    monkeypatch.setattr(np.random, "randn", lambda n: np.array([-1., 1.]))
    return particles, proc, beta, gamma, theta


@pytest.mark.parametrize("energy", [0.001, 0.1])
@pytest.mark.parametrize("xi", [0., 1e-20, 0.01, 0.1, 0.5])
def test_nagaitsev_against_paper_rate(monkeypatch, energy, xi):
    particles, proc, beta, gamma, theta = fixed_parameters(monkeypatch, energy, xi)
    # Exact F, Nagaitsev Eqs. (8), (48), independent of the approximate G in IBS.
    exact_f = (1. if xi == 0 else
               4 / (3 * np.pi) * (elliprd(xi, 1., 1.) - xi * elliprd(1., 1., xi)))
    # Eq. (46): ring rate of relative MOMENTUM variance per unit time.
    # Its bunch length is spatial; Ocelot's longitudinal coordinate is c*delta_t.
    sigma_s = beta * proc.sigma_z
    ring_rate = (np.sum(particles.q_array) / q_e * ro_e**2 * speed_of_light * proc.Clog /
                 (8 * beta**3 * gamma**3 * sigma_s * proc.sigma_xy**2 * theta) * exact_f)
    dz = 0.7
    expected = 2 * ring_rate * dz / (beta * speed_of_light) * (beta**2 * gamma)**2
    proc.apply(particles, dz)
    actual = np.var(particles.p()) * (particles.p0c / m_e_GeV)**2
    # The quoted G approximation error is 20%; its cold limit is much tighter.
    np.testing.assert_allclose(actual, expected, rtol=2e-5 if xi < 1e-15 else 0.2, atol=0)


@pytest.mark.parametrize("energy", [0.001, 0.1])
@pytest.mark.parametrize("xi", [0., 1e-20])
def test_nagaitsev_cold_limit_matches_huang(monkeypatch, energy, xi):
    particles, proc, _, _, _ = fixed_parameters(monkeypatch, energy, xi)
    proc.apply(particles, 1.)
    nagaitsev_variance = np.var(particles.p())
    particles.p()[:] = 0
    proc.method = "Huang"
    proc.apply(particles, 1.)
    np.testing.assert_allclose(nagaitsev_variance, np.var(particles.p()), rtol=2e-5, atol=0)
    assert nagaitsev_variance > 0


def test_huang_published_examples(monkeypatch):
    particles, proc, beta, _, _ = fixed_parameters(monkeypatch, 100., 0., method="Huang")
    particles.q_array[:] = 6.2e9 * q_e / particles.n
    proc.sigma_xy = 100e-6
    proc.emit_n = 1e-6
    for sigma_z, initial, length, expected in [(830e-6, 0.006, 40, 0.00913528),
                                               (190e-6, 0.04, 360, 0.05886992)]:
        particles.p()[:] = 0
        proc.sigma_z = sigma_z / beta
        proc.apply(particles, length)
        added_variance = np.var(particles.p()) * (particles.p0c / m_e_GeV)**2
        np.testing.assert_allclose(np.sqrt(initial**2 + added_variance), expected, rtol=1e-6)


@pytest.mark.parametrize("xi", [2., 10.])
def test_cooling_rejected_without_modifying_particles(monkeypatch, xi):
    particles, proc, _, _, _ = fixed_parameters(monkeypatch, 0.1, xi)
    before = particles.rparticles.copy()
    with pytest.raises(ValueError, match="Longitudinal cooling"):
        proc.apply(particles, 1.)
    np.testing.assert_array_equal(particles.rparticles, before)


def test_equal_temperatures_give_no_net_heating(monkeypatch):
    particles, proc, _, _, _ = fixed_parameters(monkeypatch, 0.1, 1.)
    proc.apply(particles, 1.)
    np.testing.assert_allclose(particles.p(), 0., atol=1e-12)


@pytest.mark.parametrize("spread", [0., 1e-14])
def test_cold_particle_array_is_heated(beam, spread, monkeypatch):
    beam.p()[:] *= spread / 1e-5
    before = beam.p().copy()
    rng = np.random.default_rng(8)
    monkeypatch.setattr(np.random, "randn", rng.standard_normal)
    with np.errstate(all="raise"):
        IBS(method="Nagaitsev", update_Clog=False).apply(beam, 1.)
    assert np.all(np.isfinite(beam.rparticles))
    assert np.std(beam.p() - before) > 1e-8


@pytest.mark.parametrize("bounds,slice_name", [(None, None), ([-0.5, 0.5], "Imax"),
                                               ([-0.5, 0.5], None)])
def test_chirp_fit_preserves_particles_and_intrinsic_spread(beam, bounds, slice_name):
    proc = IBS(bounds=bounds, slice=slice_name)
    proc.get_beam_params(beam)
    original_spread = proc.sigma_dgamma0
    beam.p()[:] += 0.02 * beam.tau() / np.std(beam.tau()) + 0.001
    before = beam.rparticles.copy()
    proc.get_beam_params(beam)
    np.testing.assert_allclose(proc.sigma_dgamma0, original_spread, rtol=1e-12)
    np.testing.assert_array_equal(beam.rparticles, before)


def test_spread_is_measured_in_selected_window(beam):
    sigma_z = np.std(beam.tau())
    outside = np.abs(beam.tau() - np.mean(beam.tau())) > 0.5 * sigma_z
    beam.p()[outside] *= 100
    proc = IBS(slice=None)
    proc.get_beam_params(beam)
    # Huge tail spread must not enter the central-window estimate.
    assert proc.sigma_dgamma0 * m_e_GeV / beam.p0c < 1.1e-5


@pytest.mark.parametrize("method", ["Huang", "Nagaitsev"])
def test_zero_step_preserves_beam_and_rng(beam, method, monkeypatch):
    beam.p()[:] = 0
    before = beam.rparticles.copy()
    def unexpected_draw(*args):
        pytest.fail("Zero IBS step must not draw random numbers")
    monkeypatch.setattr(np.random, "randn", unexpected_draw)
    IBS(method=method).apply(beam, 0.)
    np.testing.assert_array_equal(beam.rparticles, before)


@pytest.mark.parametrize("update", [True, False])
def test_effective_coulomb_log(beam, update):
    proc = IBS(Clog=8, update_Clog=update)
    proc.apply(beam, 1.)
    assert proc.Clog == 8
    assert proc.Clog_used == (proc.estimate_Clog() if update else 8)
    if update:
        assert proc.Clog_used > 13


@pytest.mark.parametrize("dz", [-1., np.nan, np.inf])
def test_invalid_step_rejected(beam, dz):
    before = beam.rparticles.copy()
    with pytest.raises(ValueError, match="path length"):
        IBS().apply(beam, dz)
    np.testing.assert_array_equal(beam.rparticles, before)


@pytest.mark.parametrize("clog", [-1., np.nan, np.inf])
def test_invalid_coulomb_log_rejected(beam, clog):
    before = beam.rparticles.copy()
    with pytest.raises(ValueError, match="Coulomb logarithm"):
        IBS(Clog=clog, update_Clog=False).apply(beam, 1.)
    np.testing.assert_array_equal(beam.rparticles, before)


@pytest.mark.parametrize("bounds", [[1, -1], [0], [0, np.inf], [100, 101]])
def test_invalid_analysis_window_rejected(beam, bounds):
    before = beam.rparticles.copy()
    with pytest.raises(ValueError, match="bounds|analysis window"):
        IBS(bounds=bounds).apply(beam, 1.)
    np.testing.assert_array_equal(beam.rparticles, before)


@pytest.mark.parametrize("problem,match", [("energy", "reference energy"),
                                           ("length", "bunch length"),
                                           ("emittance", "emittances"),
                                           ("coordinates", "finite coordinates")])
def test_invalid_beam_rejected(beam, problem, match):
    if problem == "energy":
        beam.E = m_e_GeV
    elif problem == "length":
        beam.tau()[:] = 0
    elif problem == "emittance":
        beam.px()[:] = 0
    else:
        beam.x()[0] = np.nan
    before = beam.rparticles.copy()
    with pytest.raises(ValueError, match=match):
        IBS().apply(beam, 1.)
    np.testing.assert_array_equal(beam.rparticles, before)


def test_unknown_method_rejected(beam):
    with pytest.raises(ValueError, match="Invalid method"):
        IBS(method="unknown").apply(beam, 1.)


def test_zero_charge_is_no_op(beam):
    beam.q_array[:] = 0
    before = beam.rparticles.copy()
    IBS().apply(beam, 1.)
    np.testing.assert_array_equal(beam.rparticles, before)
