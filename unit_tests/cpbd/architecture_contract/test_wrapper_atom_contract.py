import numpy as np

from ocelot.cpbd.elements import Cavity, Quadrupole
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


def test_quadrupole_wrapper_forwards_attribute_reads_and_writes_to_atom():
    quad = Quadrupole(l=0.4, k1=1.2, eid="Q1")

    assert quad.k1 == quad.element.k1

    quad.k1 = -0.7

    assert quad.k1 == -0.7
    assert quad.element.k1 == -0.7


def test_quadrupole_parameter_change_rebuilds_first_order_and_active_maps():
    quad = Quadrupole(l=0.4, k1=1.2, tm=SecondTM, eid="Q2")
    energy = 1.3

    first_order_before = quad.first_order_tms
    active_before = quad.tms
    r_before = quad.R(energy)[0].copy()
    t_before = quad.T(energy)[0].copy()

    quad.k1 = -0.3

    first_order_after = quad.first_order_tms
    active_after = quad.tms

    assert first_order_after is not first_order_before
    assert active_after is not active_before
    assert all(isinstance(tm, TransferMap) for tm in first_order_after)
    assert all(isinstance(tm, SecondTM) for tm in active_after)
    assert not np.allclose(r_before, quad.R(energy)[0])
    assert not np.allclose(t_before, quad.T(energy)[0])


def test_cavity_wrapper_passes_plotting_kwargs_to_atom():
    cavity = Cavity(width=0.2, height=0.3, color="red", eid="C_STYLE")

    assert cavity.width == 0.2
    assert cavity.height == 0.3
    assert cavity.color == "red"
