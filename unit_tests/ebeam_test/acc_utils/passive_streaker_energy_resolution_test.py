"""Focused tests for passive-streaker resolution utilities."""

from types import SimpleNamespace

import numpy as np

from ocelot.utils.acc_utils import (
    convolve_beam,
    passive_streaker_energy_effects,
    passive_streaker_energy_resolution,
    passive_streaker_resolutions,
    single_plane_dipole_wake,
    single_plate_quadrupole_wake,
)
from ocelot.common.globals import Z0, speed_of_light


def test_passive_streaker_energy_effects_matches_appendix_a():
    w_l = np.array([1.0, 2.0])
    w_ld = np.array([3.0, 4.0])
    w_lq = np.array([5.0, 6.0])
    sigma_x = 0.2
    sigma_y = 0.1
    length = 2.0

    mean, spread = passive_streaker_energy_effects(
        w_l, w_ld, w_lq, sigma_x, sigma_y, length
    )

    expected_mean = length * (
        w_l + 0.5 * (sigma_x ** 2 - sigma_y ** 2) * w_lq
    )
    expected_spread = length * np.sqrt(
        (w_ld * sigma_y) ** 2
        + 0.5 * w_lq ** 2 * (sigma_x ** 4 + sigma_y ** 4)
    )
    np.testing.assert_allclose(mean, expected_mean)
    np.testing.assert_allclose(spread, expected_spread)


def test_round_beam_has_no_quadrupole_mean_energy_correction():
    w_l = np.array([1.0, 2.0])
    w_lq = np.array([100.0, 200.0])
    mean, _ = passive_streaker_energy_effects(
        w_l, np.zeros(2), w_lq, sigma_x=0.1, sigma_y=0.1, length=3.0
    )
    np.testing.assert_allclose(mean, 3.0 * w_l)


def test_passive_streaker_energy_resolution_uses_outgoing_chirp():
    spectrometer = np.array([3.0, 4.0])
    time_resolution = np.array([2.0, 2.0])
    incoming_chirp = np.array([1.0, -1.0])
    induced_spread = np.array([4.0, 0.0])
    ps_chirp = np.array([0.5, 0.25])

    result = passive_streaker_energy_resolution(
        spectrometer,
        time_resolution,
        incoming_chirp,
        induced_spread,
        ps_chirp,
    )
    expected = np.sqrt(
        spectrometer ** 2
        + induced_spread ** 2
        + (time_resolution * (incoming_chirp + ps_chirp)) ** 2
    )
    np.testing.assert_allclose(result, expected)


def test_passive_streaker_energy_resolution_optional_terms_default_to_zero():
    result = passive_streaker_energy_resolution(3.0, 2.0, 1.0)
    np.testing.assert_allclose(result, np.sqrt(13.0))


def test_convolve_beam_does_not_modify_current():
    current = np.array([[2.0, 1.0], [3.0, 2.0], [4.0, 3.0]])
    original = current.copy()
    convolve_beam(current, lambda s: 1.0)
    np.testing.assert_array_equal(current, original)


def test_single_plate_wakes_use_qin_corrected_characteristic_scales():
    """Check Qin et al., PRAB 26, 064402 (2023), Eqs. (60)-(66)."""

    p = 0.5e-3
    t = 0.25e-3
    b = 500e-6
    length = 5.0
    s = np.linspace(0.0, 50e-6, 20)
    alpha = 1.0 - 0.465 * np.sqrt(t / p) - 0.070 * t / p
    prefactor = Z0 * speed_of_light / (4.0 * np.pi)

    s_monopole = 8.0 * b ** 2 * t / (9.0 * np.pi * alpha ** 2 * p ** 2)
    u_monopole = np.sqrt(s / s_monopole)
    expected_monopole = (
        length
        * 2.0
        / b ** 3
        * s_monopole
        * (1.0 - (1.0 + u_monopole) * np.exp(-u_monopole))
        * prefactor
    )

    s_quadrupole = b ** 2 * t / (2.0 * np.pi * alpha ** 2 * p ** 2)
    u_quadrupole = np.sqrt(s / s_quadrupole)
    expected_quadrupole = (
        length
        * 3.0
        / b ** 4
        * s_quadrupole
        * (1.0 - (1.0 + u_quadrupole) * np.exp(-u_quadrupole))
        * prefactor
    )

    monopole_wake = single_plane_dipole_wake(p=p, t=t, b=b, l=length)
    quadrupole_wake = single_plate_quadrupole_wake(p=p, t=t, b=b, l=length)
    np.testing.assert_allclose(monopole_wake(s), expected_monopole)
    np.testing.assert_allclose(quadrupole_wake(s), expected_quadrupole)


def test_legacy_energy_resolution_is_nonnegative_for_negative_dispersion():
    s = np.array([0.0, 1.0, 2.0])
    dipole_kick = np.column_stack((s, np.array([0.0, 1.0, 2.0])))
    quadrupole_kick = np.column_stack((s, np.zeros(3)))
    r_matrix = np.eye(6)
    r_matrix[0, 5] = -0.5
    r_matrix[2, 3] = 1.0
    twiss = SimpleNamespace(
        beta_x=1.0,
        alpha_x=0.0,
        gamma_x=1.0,
        beta_y=1.0,
        alpha_y=0.0,
        gamma_y=1.0,
    )

    _, energy_resolution, _, _ = passive_streaker_resolutions(
        dipole_kick,
        quadrupole_kick,
        r_matrix,
        twiss,
        kick="vert",
        energy=1.0,
    )
    assert np.all(energy_resolution[:, 1] >= 0.0)
