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
    assert drift.tm_policy == "generic"


def test_cavity_normalizes_unsupported_tm_to_cavity_tm():
    with pytest.warns(UserWarning, match="pins active tracking to CavityTM"):
        cavity = Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, tm=TransferMap, eid="C1")

    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)
    assert cavity.tm_policy == "pinned"


def test_cavity_keeps_first_order_maps_for_optics():
    cavity = Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, eid="C2")

    assert all(isinstance(tm, TransferMap) for tm in cavity.first_order_tms)
    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)


def test_multipole_keeps_multipole_tm_as_active_method():
    with pytest.warns(UserWarning, match="pins active tracking to MultipoleTM"):
        multipole = Multipole(kn=[0.2, -0.5], tm=TransferMap, eid="M1")

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)
    assert all(isinstance(tm, TransferMap) for tm in multipole.first_order_tms)
    assert multipole.tm_policy == "pinned"

    with pytest.warns(UserWarning, match="pins active tracking to MultipoleTM"):
        multipole.set_tm(TransferMap)

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)


def test_twcavity_uses_shared_pinned_tm_policy():
    with pytest.warns(UserWarning, match="pins active tracking to TWCavityTM"):
        cavity = TWCavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, tm=TransferMap, eid="TW1")

    assert cavity.tm_policy == "pinned"
    assert all(isinstance(tm, TWCavityTM) for tm in cavity.tms)


def test_xyquadrupole_uses_shared_pinned_tm_policy():
    quad = XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-2e-3, k1=0.4, eid="XY1")

    with pytest.warns(UserWarning, match="pins active tracking to TransferMap"):
        quad.set_tm(SecondTM)

    assert quad.tm_policy == "pinned"
    assert all(isinstance(tm, TransferMap) for tm in quad.tms)
