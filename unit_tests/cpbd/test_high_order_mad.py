"""Independent checks of the MAD sector-body map, including issue #326."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.high_order import _mad_focusing, _mad_mixed_integrals, t_nnn_mad


def integrated_body_tensor(length, h, k1, k2, energy):
    """Integrate the Taylor-map equations of the cubic sector Hamiltonian.

    This uses the local equations of motion, not the closed-form J integrals.
    Coordinates are Ocelot's (x, px, y, py, tau, dE/(p0*c)); tau has the
    opposite sign to MAD's time coordinate. See MAD8 Physical Methods Manual,
    sections 2.1 and 5.3.2.
    """
    igamma2 = (m_e_GeV / energy) ** 2 if energy else 0.
    beta = np.sqrt(1 - igamma2)
    linear = np.zeros((6, 6))
    linear[0, 1] = linear[2, 3] = 1
    linear[1, 0] = -(k1 + h * h)
    linear[1, 5] = linear[4, 0] = h / beta
    linear[3, 2] = k1
    linear[4, 5] = -igamma2 / beta ** 2

    quadratic = np.zeros((6, 6, 6))

    def coefficient(i, j, k, value):
        quadratic[i, j, k] = quadratic[i, k, j] = (
            value if j == k else value / 2
        )

    coefficient(0, 0, 1, h)
    coefficient(0, 1, 5, -1 / beta)
    coefficient(1, 0, 0, -(k2 + 2 * h * k1) / 2)
    coefficient(1, 1, 1, -h / 2)
    coefficient(1, 2, 2, (k2 + h * k1) / 2)
    coefficient(1, 3, 3, -h / 2)
    coefficient(1, 5, 5, -h * igamma2 / (2 * beta ** 2))
    coefficient(2, 0, 3, h)
    coefficient(2, 3, 5, -1 / beta)
    coefficient(3, 0, 2, k2 + h * k1)
    coefficient(4, 1, 1, 1 / (2 * beta))
    coefficient(4, 3, 3, 1 / (2 * beta))
    coefficient(4, 0, 5, -h * igamma2 / beta ** 2)
    coefficient(4, 5, 5, 1.5 * igamma2 / beta ** 3)

    def rhs(_, state):
        r = state[:36].reshape(6, 6)
        t = state[36:].reshape(6, 6, 6)
        dr = linear @ r
        dt = np.einsum('ia,ajk->ijk', linear, t)
        dt += np.einsum('iab,aj,bk->ijk', quadratic, r, r)
        return np.concatenate((dr.ravel(), dt.ravel()))

    solution = solve_ivp(
        rhs, (0, length), np.concatenate((np.eye(6).ravel(), np.zeros(216))),
        method='DOP853', rtol=2e-13, atol=2e-15,
    )
    assert solution.success
    symmetric = solution.y[36:, -1].reshape(6, 6, 6)
    # Ocelot stores each mixed monomial once, in the upper triangle.
    result = np.zeros_like(symmetric)
    for j in range(6):
        result[:, j, j] = symmetric[:, j, j]
        for k in range(j + 1, 6):
            result[:, j, k] = symmetric[:, j, k] + symmetric[:, k, j]
    return result


@pytest.mark.parametrize('energy', [0., 0.002])
@pytest.mark.parametrize('length,h,k1,k2', [
    (0., 1., -0.2, 0.7),
    (1.3, 0., 0., 0.),
    (0.7, 0., 0.8, 0.),
    (1.3, 0., -0.8, 0.),
    (1.3, 0., 0., 0.7),
    (1.3, 0.2, 0., 0.),
    (0.7, 0.3, 0.4, 0.7),
    (1.3, 0.3, -0.4, 0.7),
    (1.3, 0.3, 0., 0.7),
    (1., 1., -0.2, 0.7),  # kx2 = 4*ky2 != 0
    (1., 1., -0.2 + 1e-12, 0.7),
    (1., 1., -0.2 - 1e-12, 0.7),
    (1., 1., -1., 0.7),  # kx2 = 0, ky2 != 0
    (1., 1., -1. + 1e-10, 0.7),
    (1., 1., -1. - 1e-10, 0.7),
    (1., 1e-5, -2e-11, 0.7),  # weak, resonant focusing
    (1e-5, 1., -0.2, 0.7),
    (-0.7, 0.3, 0.4, 0.7),
], ids=[
    'zero-length', 'drift', 'focusing-quad', 'defocusing-quad', 'sextupole',
    'dipole', 'combined-focusing', 'combined-defocusing', 'bend-sextupole',
    'resonance', 'above-resonance', 'below-resonance', 'zero-kx2',
    'small-positive-kx2', 'small-negative-kx2', 'weak-resonance',
    'short-resonance', 'negative-length',
])
def test_mad_tensor_against_integrated_hamiltonian(length, h, k1, k2, energy):
    expected = integrated_body_tensor(length, h, k1, k2, energy)
    actual = t_nnn_mad(length, h, k1, k2, energy)
    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-12)


@pytest.mark.parametrize('a', [-0.100001, -0.1, -0.099999, 0.099999, 0.1, 0.100001])
def test_mad_tensor_across_series_boundary(a):
    # A nonzero bend keeps the J3 contribution observable as kx2 changes sign.
    h, k1, k2 = 0.7, a - 0.7 ** 2, 0.8
    expected = integrated_body_tensor(1., h, k1, k2, 0.)
    np.testing.assert_allclose(
        t_nnn_mad(1., h, k1, k2), expected, rtol=2e-10, atol=2e-12,
    )


@pytest.mark.parametrize('length', [0.03, 1., 2.7])
@pytest.mark.parametrize('a,b', [
    (0., 0.), (1e-12, 1e-12), (0.099999, 0.1), (0.1, 0.1),
    (0.100001, 0.1), (0.8, 0.8), (0.8 + 1e-12, 0.8),
    (0.8 - 1e-12, 0.8), (-0.8, 1.2), (0.8, -0.2),
])
def test_mixed_integrals_against_green_function_ode(length, a, b):
    # Work in normalized length so that small Jf values cannot pass solely
    # because an absolute tolerance masks a wrong power of the element length.
    def rhs(s, y):
        return [y[1], np.cos(np.sqrt(b + 0.j) * s).real - a * y[0],
                y[0], y[2], y[3]]

    solution = solve_ivp(
        rhs, (0, 1), np.zeros(5), method='DOP853', rtol=2e-13, atol=2e-15,
    )
    assert solution.success
    expected = solution.y[[0, 2, 3, 4], -1]
    kx2, ky2 = a / length ** 2, b / (4 * length ** 2)
    actual = _mad_mixed_integrals(length, kx2, ky2, *_mad_focusing(length, kx2))
    normalized = np.asarray(actual) / length ** np.arange(2, 6)
    np.testing.assert_allclose(normalized, expected, rtol=2e-11, atol=2e-14)
