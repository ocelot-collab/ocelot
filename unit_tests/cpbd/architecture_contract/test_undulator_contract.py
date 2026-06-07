from ocelot.cpbd.elements import Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.undulator_test import UndulatorTestTM


def test_undulator_can_switch_to_runge_kutta_and_keep_first_order_maps():
    und = Undulator(lperiod=0.03, nperiods=10, Kx=1.2, eid="U1")

    und.set_tm(RungeKuttaTM, npoints=25)

    assert all(isinstance(tm, RungeKuttaTM) for tm in und.tms)
    assert all(isinstance(tm, TransferMap) for tm in und.first_order_tms)
    assert und.tms[0].npoints == 25


def test_undulator_can_construct_undulator_test_tm():
    und = Undulator(lperiod=0.03, nperiods=10, Kx=1.2, eid="U2")

    und.set_tm(UndulatorTestTM, ndiv=7)

    assert all(isinstance(tm, UndulatorTestTM) for tm in und.tms)
    assert all(isinstance(tm, TransferMap) for tm in und.first_order_tms)
    assert und.tms[0].ndiv == 7


def test_magnetic_lattice_can_select_runge_kutta_for_undulator_family():
    und = Undulator(lperiod=0.03, nperiods=10, Kx=1.2, eid="U3")

    MagneticLattice((und,), method={Undulator: RungeKuttaTM})

    assert all(isinstance(tm, RungeKuttaTM) for tm in und.tms)
    assert all(isinstance(tm, TransferMap) for tm in und.first_order_tms)


def test_undulator_declares_second_order_tracking():
    und = Undulator(lperiod=0.03, nperiods=10, Kx=1.2, eid="U4")

    und.set_tm(SecondTM)

    assert all(isinstance(tm, SecondTM) for tm in und.tms)
    assert all(isinstance(tm, TransferMap) for tm in und.first_order_tms)


def test_global_second_order_request_keeps_undulator_on_second_tm():
    und = Undulator(lperiod=0.03, nperiods=10, Kx=1.2, eid="U5")

    MagneticLattice((und,), method={"global": SecondTM})

    assert all(isinstance(tm, SecondTM) for tm in und.tms)
