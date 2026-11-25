"""Spontaneous radiation moments and recoil in undulators and uniform bends."""

from inspect import Parameter, signature

import numpy as np
import pytest

from ocelot.common.globals import lambda_C_r, m_e_GeV, ro_e
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.physics_proc import SpontanRadEffects


def test_constructor_accepts_named_options_and_rejects_typos():
    options = signature(SpontanRadEffects).parameters
    defaults = dict(energy_loss=True, quant_diff=True, filling_coeff=1.0, radius=np.inf, step=1)
    for name, value in defaults.items():
        assert options[name].kind == Parameter.KEYWORD_ONLY
        assert options[name].default == value
    assert SpontanRadEffects().step == 1
    process = SpontanRadEffects(3.7, 0.04, "helical", energy_loss=False, quant_diff=True,
                                filling_coeff=0.8, step=2)
    assert (process.K, process.lperiod, process.type) == (3.7, 0.04, "helical")
    assert (process.energy_loss, process.quant_diff, process.filling_coeff, process.step) == (False, True, 0.8, 2)
    with pytest.raises(TypeError, match="quant_dif"):
        SpontanRadEffects(quant_dif=False)


@pytest.fixture
def cold_beam(monkeypatch):
    particles = ParticleArray(n=2)
    particles.E = 10.0  # GeV
    # Two standardized kicks isolate the mean and variance without sampling error.
    monkeypatch.setattr(np.random, "randn", lambda n: np.linspace(-1., 1., n))
    return particles


def assert_energy_kick_moments(particles, mean_loss, added_variance):
    relative_energy_change = particles.p() * particles.p0c / particles.E
    np.testing.assert_allclose(-relative_energy_change.mean(), mean_loss, atol=0, rtol=1e-8)
    np.testing.assert_allclose(relative_energy_change.var(), added_variance, atol=0, rtol=1e-8)


@pytest.mark.parametrize("geometry", ["planar", "helical"])
@pytest.mark.parametrize("K", [0.2, -0.2, 3.5, -3.5])
@pytest.mark.parametrize("filling", [1., 0.25])
def test_undulator_energy_kick_moments(cold_beam, geometry, K, filling):
    period = 0.04  # m
    length = 2.0  # m
    process = SpontanRadEffects(type=geometry, K=K, lperiod=period, filling_coeff=filling)
    process.apply(cold_beam, length)

    gamma = cold_beam.E / m_e_GeV
    ku = 2 * np.pi / period
    active_length = filling * length
    peak_curvature = abs(K) * ku / gamma
    mean_curvature_squared = peak_curvature**2 * (0.5 if geometry == "planar" else 1.)
    mean_loss = (2 / 3) * ro_e * gamma**3 * mean_curvature_squared * active_length

    # Saldin, Schneidmiller and Yurkov, NIM A 381 (1996) 545, peak K convention.
    k = abs(K)
    if geometry == "planar":
        fit = 1.20 * k + 1 / (1 + 1.33 * k + 0.40 * k**2)
        coefficient = 7 / 15
    else:
        fit = 1.42 * k + 1 / (1 + 1.50 * k + 0.95 * k**2)
        coefficient = 14 / 15
    added_variance = coefficient * lambda_C_r * ro_e * gamma**2 * ku**3 * k**2 * fit * active_length
    assert_energy_kick_moments(cold_beam, mean_loss, added_variance)


@pytest.mark.parametrize("radius", [100., -100.])
@pytest.mark.parametrize("filling", [1., 0.25])
def test_dipole_energy_kick_moments(cold_beam, radius, filling):
    length = 1.0  # m
    process = SpontanRadEffects(type="dipole", radius=radius, filling_coeff=filling)
    process.apply(cold_beam, length)

    # Yampolsky and Carlsten, arXiv:1207.4833, Eqs. (17)–(20).
    gamma = cold_beam.E / m_e_GeV
    active_length = filling * length
    mean_loss = (2 / 3) * ro_e * gamma**3 * active_length / radius**2
    added_variance = (55 / (24 * np.sqrt(3))) * lambda_C_r * ro_e * gamma**5 * active_length / abs(radius)**3
    assert_energy_kick_moments(cold_beam, mean_loss, added_variance)
    assert (process.K, process.lperiod) == (0., 0.)


@pytest.fixture
def chirped_beam():
    beam = ParticleArray(n=3)
    beam.E = 0.15  # GeV
    beam.tau()[:] = [-0.001, 0., 0.001]
    beam.p()[:] = [-0.2, 0.1, 0.4]
    beam.px()[:] = [0.004, -0.002, 0.003]
    beam.py()[:] = [-0.002, 0.003, 0.004]
    return beam


@pytest.mark.parametrize("device", ["planar", "dipole"])
def test_mean_loss_preserves_direction_in_a_chirped_beam(chirped_beam, device):
    beam = chirped_beam
    length = 0.8
    if device == "dipole":
        radius = 0.6  # design radius in meters at the reference momentum
        radiation = SpontanRadEffects(type=device, radius=radius, quant_diff=False)
        curvature_squared = (beam.p0c / (radius * beam.momenta)) ** 2
    else:
        period = 0.03
        K = 4.
        radiation = SpontanRadEffects(type=device, lperiod=period, K=K, quant_diff=False)
        curvature_squared = 0.5 * (2 * np.pi * K / (period * beam.gamma)) ** 2
    # Integrate the classical power using each particle's mean squared curvature.
    energy_loss = (2 / 3) * m_e_GeV * ro_e * beam.gamma**4 * curvature_squared * length
    p_before = beam.p().copy()
    p0c = beam.p0c
    slopes_before = np.array([beam.px(), beam.py()]) / beam.pz

    radiation.apply(beam, length)

    assert beam.E == 0.15
    np.testing.assert_allclose(p_before - beam.p(), energy_loss / p0c, rtol=1e-7, atol=0)
    slopes_after = np.array([beam.px(), beam.py()]) / beam.pz
    np.testing.assert_allclose(slopes_after, slopes_before, rtol=5e-15, atol=0)


def test_off_reference_dipole_diffusion_uses_fixed_field(chirped_beam, monkeypatch):
    beam = chirped_beam
    length = 0.8
    radius = 0.6
    particle_radii = radius * beam.momenta / beam.p0c
    variance_gamma = (55 / (24 * np.sqrt(3)) * ro_e * lambda_C_r * beam.gamma**7
                      * length / particle_radii**3)
    draws = np.array([-1.5, 0.5, 2.])
    monkeypatch.setattr(np.random, "randn", lambda n: draws.copy())
    expected_kicks = draws * m_e_GeV * np.sqrt(variance_gamma) / beam.p0c
    p_before = beam.p().copy()
    slopes_before = np.array([beam.px(), beam.py()]) / beam.pz

    SpontanRadEffects(type="dipole", radius=radius, energy_loss=False).apply(beam, length)

    assert beam.E == 0.15
    np.testing.assert_allclose(beam.p() - p_before, expected_kicks, rtol=1e-8, atol=0)
    slopes_after = np.array([beam.px(), beam.py()]) / beam.pz
    np.testing.assert_allclose(slopes_after, slopes_before, rtol=5e-15, atol=0)


@pytest.mark.parametrize("kwargs", [
    {}, {"type": "dipole", "radius": np.inf},
    {"K": 3.7, "lperiod": 0.04, "filling_coeff": 0.},
])
def test_no_active_field_leaves_beam_untouched(chirped_beam, kwargs):
    before = chirped_beam.rparticles.copy()
    SpontanRadEffects(**kwargs).apply(chirped_beam, 1.)
    np.testing.assert_array_equal(chirped_beam.rparticles, before)


def test_excessive_energy_loss_leaves_beam_untouched(chirped_beam):
    before = chirped_beam.rparticles.copy()
    process = SpontanRadEffects(type="dipole", radius=1.e-6, quant_diff=False)
    with pytest.raises(ValueError, match="step length"):
        process.apply(chirped_beam, 1.)
    np.testing.assert_array_equal(chirped_beam.rparticles, before)
