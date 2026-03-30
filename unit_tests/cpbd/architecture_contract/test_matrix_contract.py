import numpy as np

from ocelot.cpbd.elements import Matrix
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transformation import TMTypes


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


def test_matrix_slice_uses_drift_like_first_order_map_instead_of_slicing_stored_r():
    matrix = Matrix(
        l=2.0,
        delta_e=0.05,
        eid="MX2",
        R11=1.0,
        R12=3.0,
        R22=1.0,
        R33=1.0,
        R34=4.0,
        R44=1.0,
        R55=1.0,
        R66=1.0,
    )
    energy = 1.7

    section = matrix.get_section_tms(start_l=0.4, delta_l=0.5)
    params = section[0].get_params(energy)
    expected_r = uni_matrix(0.5, 0.0, hx=0.0)

    assert section[0].tm_type == TMTypes.MAIN
    np.testing.assert_allclose(params.R, expected_r)
    assert params.R[0, 1] == 0.5
    assert params.R[2, 3] == 0.5

