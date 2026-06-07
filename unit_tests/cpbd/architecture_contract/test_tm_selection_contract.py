import pytest

from ocelot.cpbd.elements import Cavity, Drift, Multipole, TWCavity, XYQuadrupole
from ocelot.cpbd.elements.cavity import CavityTM
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM


def test_drift_can_switch_to_second_order_tracking():
    drift = Drift(l=0.35, tm=TransferMap)

    drift.set_tm(SecondTM)

    assert all(isinstance(tm, SecondTM) for tm in drift.tms)


def test_cavity_rejects_unsupported_explicit_tm_request():
    with pytest.raises(RuntimeError, match="Cavity does not declare support for TransferMap"):
        Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, tm=TransferMap, eid="C1")


def test_cavity_keeps_first_order_maps_for_optics():
    cavity = Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, eid="C2")

    assert all(isinstance(tm, TransferMap) for tm in cavity.first_order_tms)
    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)


def test_multipole_keeps_multipole_tm_as_active_method():
    multipole = Multipole(kn=[0.2, -0.5], eid="M1")

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)
    assert all(isinstance(tm, TransferMap) for tm in multipole.first_order_tms)

    with pytest.raises(RuntimeError, match="Multipole does not declare support for TransferMap"):
        multipole.set_tm(TransferMap)

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)


def test_twcavity_rejects_unsupported_explicit_tm_request():
    with pytest.raises(RuntimeError, match="TWCavity does not declare support for TransferMap"):
        TWCavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, tm=TransferMap, eid="TW1")

    cavity = TWCavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, eid="TW2")
    assert all(isinstance(tm, TWCavityTM) for tm in cavity.tms)


def test_xyquadrupole_rejects_unsupported_explicit_tm_request():
    quad = XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-2e-3, k1=0.4, eid="XY1")

    with pytest.raises(RuntimeError, match="XYQuadrupole does not declare support for SecondTM"):
        quad.set_tm(SecondTM)

    assert all(isinstance(tm, TransferMap) for tm in quad.tms)
