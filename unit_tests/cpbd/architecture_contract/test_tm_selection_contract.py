from ocelot.cpbd.elements import Cavity, Drift, Multipole
from ocelot.cpbd.elements.cavity import CavityTM
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


def test_drift_can_switch_to_second_order_tracking():
    drift = Drift(l=0.35, tm=TransferMap)

    drift.set_tm(SecondTM)

    assert all(isinstance(tm, SecondTM) for tm in drift.tms)


def test_cavity_normalizes_unsupported_tm_to_cavity_tm():
    cavity = Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, tm=TransferMap, eid="C1")

    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)


def test_cavity_keeps_first_order_maps_for_optics():
    cavity = Cavity(l=0.346, v=0.0025, freq=3.9e9, phi=180.0, eid="C2")

    assert all(isinstance(tm, TransferMap) for tm in cavity.first_order_tms)
    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)


def test_multipole_keeps_multipole_tm_as_active_method():
    multipole = Multipole(kn=[0.2, -0.5], tm=TransferMap, eid="M1")

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)
    assert all(isinstance(tm, TransferMap) for tm in multipole.first_order_tms)

    multipole.set_tm(TransferMap)

    assert all(isinstance(tm, MultipoleTM) for tm in multipole.tms)
