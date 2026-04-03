import pytest

from ocelot.cpbd.elements import RBend


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
