import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.elements import Bend, Drift, Matrix, Quadrupole
from ocelot.cpbd.high_order import fringe_ent, fringe_ext
from ocelot.cpbd.r_matrix import bend_edge_matrix, linear_magnet_matrix, uni_matrix
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transformation import TMTypes
from ocelot.cpbd.transformations.transfer_map import TransferMap


def test_matrix_exposes_explicit_full_element_r_b_t_and_delta_e():
    matrix = Matrix(
        l=2.0,
        delta_e=0.05,
        eid="MX1",
        R11=1.0,
        R12=3.0,
        R22=1.0,
        R33=1.0,
        R34=4.0,
        R44=1.0,
        R55=1.0,
        R66=1.0,
        B1=0.2,
        T111=0.7,
        tm=SecondTM,
    )
    energy = 1.7

    r = matrix.R(energy)[0]
    b = matrix.B(energy)[0]
    t = matrix.T(energy)[0]

    assert r[0, 0] == 1.0
    assert r[0, 1] == 3.0
    assert r[2, 3] == 4.0
    assert b[0, 0] == 0.2
    assert t[0, 0, 0] == 0.7
    assert matrix.tms[0].get_delta_e() == 0.05


def test_matrix_partial_slice_requires_first_order_only():
    matrix = Matrix(l=2.0, delta_e=0.05, eid="MX2", tm=SecondTM)

    with pytest.raises(RuntimeError, match="Matrix supports partial slicing only for first_order_only=True"):
        matrix.get_section_tms(start_l=0.4, delta_l=0.5)


def test_matrix_first_order_only_slice_uses_drift_approximation_and_scaled_delta_e():
    matrix = Matrix(
        l=2.0,
        delta_e=0.05,
        eid="MX3",
        R11=1.0,
        R12=3.0,
        R22=1.0,
        R33=1.0,
        R34=4.0,
        R44=1.0,
        R55=1.0,
        R66=1.0,
        B1=0.2,
        T111=0.7,
        tm=SecondTM,
    )
    energy = 1.7

    section = matrix.get_section_tms(start_l=0.4, delta_l=0.5, first_order_only=True)
    params = section[0].get_params(energy)
    expected_r = uni_matrix(0.5, 0.0, hx=0.0, energy=energy)

    assert section[0].tm_type == TMTypes.MAIN
    assert isinstance(section[0], TransferMap)
    np.testing.assert_allclose(params.R, expected_r)
    np.testing.assert_allclose(params.B, np.zeros((6, 1)))
    assert params.R[0, 1] == 0.5
    assert params.R[2, 3] == 0.5
    assert section[0].get_delta_e() == pytest.approx(0.05 * 0.5 / 2.0)


def test_linear_magnet_matrix_handles_zero_horizontal_focusing_limit():
    energy = 1.7
    z = 0.8
    h = 0.2
    k1 = -(h * h)
    gamma = energy / m_e_GeV
    igamma2 = 1.0 / (gamma * gamma)
    beta = np.sqrt(1.0 - igamma2)
    ky = np.sqrt(-k1)
    cy = np.cos(ky * z)
    sy = np.sin(ky * z) / ky
    dx = h * z * z / 2.0
    r56 = h * h * z ** 3 / (6.0 * beta * beta) - z * igamma2 / (beta * beta)

    expected_r = np.array([
        [1.0, z, 0.0, 0.0, 0.0, dx / beta],
        [0.0, 1.0, 0.0, 0.0, 0.0, z * h / beta],
        [0.0, 0.0, cy, sy, 0.0, 0.0],
        [0.0, 0.0, k1 * sy, cy, 0.0, 0.0],
        [z * h / beta, dx / beta, 0.0, 0.0, 1.0, r56],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ])

    np.testing.assert_allclose(linear_magnet_matrix(z, k1, hx=h, energy=energy), expected_r)
    np.testing.assert_allclose(uni_matrix(z, k1, hx=h, energy=energy), expected_r)


def test_bend_edge_matrix_matches_first_order_fringe_matrix():
    params = dict(h=0.13, k1=0.02, e=0.11, h_pole=0.03, gap=0.04, fint=0.7)
    expected_r = bend_edge_matrix(params["h"], params["e"], gap=params["gap"], fint=params["fint"])

    r_ent, _ = fringe_ent(**params)
    r_ext, _ = fringe_ext(**params)

    np.testing.assert_allclose(r_ent, expected_r)
    np.testing.assert_allclose(r_ext, expected_r)


def test_linear_r_matches_existing_first_order_r_api():
    energy = 1.7
    drift = Drift(l=0.9, eid="D1")
    quad = Quadrupole(l=0.35, k1=1.4, tilt=0.17, eid="Q1")
    bend = Bend(l=1.2, angle=0.18, k1=0.04, e1=0.03, e2=0.05, gap=0.02, fint=0.6, tilt=0.11, eid="B1")

    np.testing.assert_allclose(drift.linear_r(energy=energy), drift.R(energy)[0])
    np.testing.assert_allclose(quad.linear_r(energy=energy), quad.R(energy)[0])
    np.testing.assert_allclose(bend.linear_r(energy=energy), bend.R(energy)[1])

    block_expected = bend.R(energy)
    block_actual = bend.linear_r_blocks(energy=energy)
    assert len(block_actual) == 3
    for actual, expected in zip(block_actual, block_expected):
        np.testing.assert_allclose(actual, expected)

    full_expected = block_expected[0]
    for block in block_expected[1:]:
        full_expected = block @ full_expected
    np.testing.assert_allclose(bend.linear_r_full(energy=energy), full_expected)


def test_linear_r_supports_explicit_parameter_overrides_without_mutation():
    energy = 0.0
    quad = Quadrupole(l=0.3, k1=1.2, tilt=0.1, eid="Q2")
    bend = Bend(l=1.0, angle=0.12, e1=0.02, e2=0.04, gap=0.03, fint=0.5, tilt=0.09, eid="B2")

    quad_override = quad.linear_r(energy=energy, l=0.45, k1=1.7, tilt=0.2)
    quad_expected = Quadrupole(l=0.45, k1=1.7, tilt=0.2).R(energy)[0]
    np.testing.assert_allclose(quad_override, quad_expected)
    assert quad.l == pytest.approx(0.3)
    assert quad.k1 == pytest.approx(1.2)
    assert quad.tilt == pytest.approx(0.1)

    bend_override = bend.linear_r_full(
        energy=energy,
        l=1.3,
        angle=0.16,
        e1=0.025,
        e2=0.045,
        gap=0.05,
        fint=0.7,
        fintx=0.8,
        tilt=0.12,
    )
    bend_expected_elem = Bend(
        l=1.3,
        angle=0.16,
        e1=0.025,
        e2=0.045,
        gap=0.05,
        fint=0.7,
        fintx=0.8,
        tilt=0.12,
    )
    bend_expected = bend_expected_elem.R(energy)[0]
    for block in bend_expected_elem.R(energy)[1:]:
        bend_expected = block @ bend_expected
    np.testing.assert_allclose(bend_override, bend_expected)
    assert bend.l == pytest.approx(1.0)
    assert bend.angle == pytest.approx(0.12)
    assert bend.fint == pytest.approx(0.5)
