import numpy as np
import pytest

from ocelot.cpbd.elements import RBend
from ocelot.cpbd.high_order import fringe_ent, fringe_ext
from ocelot.cpbd.r_matrix import bend_edge_matrix
from ocelot.cpbd.transformations.transformation import TMTypes


def test_rbend_angle_update_reapplies_default_half_angle_shift():
    rbend = RBend(l=1.2, angle=0.2)

    assert rbend.e1 == pytest.approx(0.1)
    assert rbend.e2 == pytest.approx(0.1)

    rbend.angle = 0.4

    assert rbend.e1 == pytest.approx(0.2)
    assert rbend.e2 == pytest.approx(0.2)


def test_rbend_angle_update_preserves_custom_edge_offsets():
    rbend = RBend(l=1.2, angle=0.2, e1=0.01, e2=0.02)

    assert rbend.e1 == pytest.approx(0.11)
    assert rbend.e2 == pytest.approx(0.12)

    rbend.angle = 0.4

    assert rbend.e1 == pytest.approx(0.21)
    assert rbend.e2 == pytest.approx(0.22)


def test_rbend_edge_first_order_matrix_is_shared_with_high_order_fringe_code():
    rbend = RBend(
        l=1.2,
        angle=0.24,
        k1=0.3,
        e1=0.01,
        e2=0.02,
        gap=0.03,
        h_pole1=0.04,
        h_pole2=0.05,
        fint=0.6,
        fintx=0.7,
    )
    h = rbend.angle / rbend.l

    entrance_r = rbend.get_tm(TMTypes.ENTRANCE, first_order_only=True).get_params(0.0).R
    exit_r = rbend.get_tm(TMTypes.EXIT, first_order_only=True).get_params(0.0).R
    high_order_entrance_r, _ = fringe_ent(
        h=h,
        k1=rbend.k1,
        e=rbend.e1,
        h_pole=rbend.h_pole1,
        gap=rbend.gap,
        fint=rbend.fint,
    )
    high_order_exit_r, _ = fringe_ext(
        h=h,
        k1=rbend.k1,
        e=rbend.e2,
        h_pole=rbend.h_pole2,
        gap=rbend.gap,
        fint=rbend.fintx,
    )

    np.testing.assert_allclose(entrance_r, bend_edge_matrix(h=h, edge=rbend.e1, gap=rbend.gap, fint=rbend.fint))
    np.testing.assert_allclose(exit_r, bend_edge_matrix(h=h, edge=rbend.e2, gap=rbend.gap, fint=rbend.fintx))
    np.testing.assert_allclose(high_order_entrance_r, entrance_r)
    np.testing.assert_allclose(high_order_exit_r, exit_r)
