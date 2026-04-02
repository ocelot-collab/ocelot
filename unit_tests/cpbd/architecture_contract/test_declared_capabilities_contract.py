import pytest

from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.cavity import Cavity, CavityTM
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.hcor import Hcor
from ocelot.cpbd.elements.matrix import Matrix
from ocelot.cpbd.elements.marker import Marker
from ocelot.cpbd.elements.monitor import Monitor
from ocelot.cpbd.elements.multipole import Multipole
from ocelot.cpbd.elements.octupole import Octupole
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.rbend import RBend
from ocelot.cpbd.elements.sbend import SBend
from ocelot.cpbd.elements.sextupole import Sextupole
from ocelot.cpbd.elements.tdcavity import TDCavity
from ocelot.cpbd.elements.twcavity import TWCavity
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.elements.vcor import Vcor
from ocelot.cpbd.elements.xyquadruple import XYQuadrupole
from ocelot.cpbd.elements.quadrupole import Quadrupole
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM
from ocelot.cpbd.transformations.undulator_test import UndulatorTestTM


DECLARED_TM_FACTORIES = [
    ("Marker", lambda: Marker(eid="MK_META"), (TransferMap, SecondTM)),
    ("Monitor", lambda: Monitor(l=0.1, eid="MON_META"), (TransferMap, SecondTM)),
    ("Drift", lambda: Drift(l=0.2, eid="D_META"), (TransferMap, SecondTM, KickTM, RungeKuttaTM, RungeKuttaTrTM)),
    ("Hcor", lambda: Hcor(l=0.2, angle=1e-4, eid="HC_META"), (TransferMap, SecondTM)),
    ("Vcor", lambda: Vcor(l=0.2, angle=1e-4, eid="VC_META"), (TransferMap, SecondTM)),
    ("Quadrupole", lambda: Quadrupole(l=0.4, k1=1.2, k2=0.3, eid="Q_META"), (TransferMap, SecondTM, KickTM)),
    ("Sextupole", lambda: Sextupole(l=0.2, k2=2.5, eid="SX_META"), (TransferMap, SecondTM, KickTM)),
    ("Octupole", lambda: Octupole(l=0.2, k3=3.5, eid="OC_META"), (TransferMap, KickTM)),
    ("Bend", lambda: Bend(l=1.1, angle=0.15, e1=0.03, e2=0.02, eid="B_META"), (TransferMap, SecondTM, RungeKuttaTM, RungeKuttaTrTM)),
    ("SBend", lambda: SBend(l=1.1, angle=0.15, e1=0.03, e2=0.02, eid="SB_META"), (TransferMap, SecondTM, RungeKuttaTM, RungeKuttaTrTM)),
    ("RBend", lambda: RBend(l=1.1, angle=0.15, e1=0.03, e2=0.02, eid="RB_META"), (TransferMap, SecondTM, RungeKuttaTM, RungeKuttaTrTM)),
    ("Matrix", lambda: Matrix(l=0.2, delta_e=0.001, eid="M_META"), (TransferMap, SecondTM)),
    ("Cavity", lambda: Cavity(l=0.3, v=0.002, freq=1.3e9, phi=0.0, eid="C_META"), (CavityTM,)),
    ("TDCavity", lambda: TDCavity(l=0.3, v=0.002, freq=2.8e9, phi=0.0, eid="TDS_META"), (TransferMap, SecondTM)),
    ("TWCavity", lambda: TWCavity(l=0.3, v=0.002, freq=1.3e9, phi=0.0, eid="TWC_META"), (TWCavityTM,)),
    ("Multipole", lambda: Multipole(kn=[0.1, 0.2], eid="MP_META"), (MultipoleTM,)),
    ("XYQuadrupole", lambda: XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-1e-3, k1=0.2, eid="XY_META"), (TransferMap,)),
    ("Undulator", lambda: Undulator(lperiod=0.03, nperiods=5, Kx=1.0, eid="U_META"), (TransferMap, RungeKuttaTM, RungeKuttaTrTM, UndulatorTestTM)),
]


class BrokenRungeKuttaAtom(Element):
    def __init__(self):
        super().__init__(eid="BROKEN")
        self.l = 0.2


class BrokenRungeKuttaElement(OpticElement):
    default_tm = TransferMap
    supported_tms = {TransferMap, RungeKuttaTM}

    def __init__(self):
        super().__init__(BrokenRungeKuttaAtom(), tm=TransferMap, default_tm=TransferMap)


def test_representative_families_expose_current_default_tm_on_instances():
    drift = Drift(l=0.2, eid="D_META")
    matrix = Matrix(l=0.2, eid="M_META")
    cavity = Cavity(l=0.3, v=0.002, freq=1.3e9, phi=0.0, eid="C_META")
    tdcavity = TDCavity(l=0.3, v=0.002, freq=2.8e9, phi=0.0, eid="TDS_META")
    multipole = Multipole(kn=[0.1, 0.2], eid="MP_META")
    xyquad = XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-1e-3, k1=0.2, eid="XY_META")
    undulator = Undulator(lperiod=0.03, nperiods=5, Kx=1.0, eid="U_META")

    assert drift.default_tm is TransferMap
    assert matrix.default_tm is TransferMap
    assert cavity.default_tm is CavityTM
    assert tdcavity.default_tm is TransferMap
    assert multipole.default_tm is MultipoleTM
    assert xyquad.default_tm is TransferMap
    assert undulator.default_tm is TransferMap


def test_representative_families_expose_tm_policy_on_instances():
    drift = Drift(l=0.2, eid="D_POLICY")
    cavity = Cavity(l=0.3, v=0.002, freq=1.3e9, phi=0.0, eid="C_POLICY")
    twcavity = TWCavity(l=0.3, v=0.002, freq=1.3e9, phi=0.0, eid="TWC_POLICY")
    multipole = Multipole(kn=[0.1, 0.2], eid="MP_POLICY")
    xyquad = XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-1e-3, k1=0.2, eid="XY_POLICY")

    assert drift.tm_policy == "generic"
    assert cavity.tm_policy == "pinned"
    assert twcavity.tm_policy == "pinned"
    assert multipole.tm_policy == "pinned"
    assert xyquad.tm_policy == "pinned"


@pytest.mark.parametrize(
    ("family_name", "factory", "declared_tms"),
    DECLARED_TM_FACTORIES,
    ids=[name for name, _, _ in DECLARED_TM_FACTORIES],
)
def test_declared_supported_tms_are_buildable(family_name, factory, declared_tms):
    element = factory()

    assert element.supported_tms == set(declared_tms)

    for tm_class in declared_tms:
        element.set_tm(tm_class)
        assert all(isinstance(tm, tm_class) for tm in element.tms), family_name


def test_undeclared_tm_request_warns_and_falls_back_to_default():
    quad = Quadrupole(l=0.3, k1=0.4, eid="Q_WARN")

    with pytest.warns(UserWarning, match="Quadrupole does not declare support for RungeKuttaTM"):
        quad.set_tm(RungeKuttaTM)

    assert all(isinstance(tm, TransferMap) for tm in quad.tms)


def test_constructor_request_for_undeclared_tm_warns_and_uses_default():
    with pytest.warns(UserWarning, match="Cavity pins active tracking to CavityTM"):
        cavity = Cavity(l=0.7, freq=1.3e9, phi=20.0, v=0.005, tm=TransferMap, eid="C_WARN")

    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)


def test_xyquadrupole_stays_first_order_only_with_explicit_fallback_warning():
    quad = XYQuadrupole(l=0.3, x_offs=1e-3, y_offs=-2e-3, k1=0.4, tm=TransferMap, eid="XY1")

    with pytest.warns(UserWarning, match="XYQuadrupole pins active tracking to TransferMap"):
        quad.set_tm(SecondTM)

    assert all(isinstance(tm, TransferMap) for tm in quad.tms)


def test_declared_but_unbuildable_tm_raises_contract_error():
    element = BrokenRungeKuttaElement()

    with pytest.raises(RuntimeError, match="declares support for RungeKuttaTM"):
        element.set_tm(RungeKuttaTM)


def test_octupole_declared_kick_tm_preserves_k3_in_kick_params():
    octupole = Octupole(l=0.25, k3=4.2, tm=KickTM, eid="OC_KICK")

    params = octupole.tms[0].get_params()

    assert params.k3 == pytest.approx(4.2)
