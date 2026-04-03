import numpy as np
import pytest

from ocelot.cpbd.elements import Matrix
from ocelot.cpbd.r_matrix import uni_matrix
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
