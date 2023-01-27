from ocelot.cpbd.elements.cavity import Cavity, CavityTM
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.second_order import SecondTM


def test_try_to_set_first_order_for_cavity():
    cavity = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.1.I1', vx_up=0.1, vy_up=0.002, vxx_up=0.1, vxy_up=0.2,
                    vx_down=0.1, vy_down=0.2, vxx_down=0.1, vxy_down=0.2, tm=TransferMap)
    for tm in cavity.tms:
        assert CavityTM == tm.__class__


def test_try_to_set_second_order_for_frift():
    drift = Drift(l=0.346, tm=TransferMap)
    drift.set_tm(SecondTM)
    for tm in drift.tms:
        assert SecondTM == tm.__class__