from ocelot.cpbd.elements import Cavity, Drift, Quadrupole
from ocelot.cpbd.elements.cavity import CavityTM
from ocelot.cpbd.latticeIO import LatticeIO
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


def test_magnetic_lattice_preserves_family_specific_tm_behavior():
    drift = Drift(l=0.3, eid="D2")
    cavity = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3.9e9, eid="C6")

    lat = MagneticLattice((drift, cavity), method={"global": SecondTM})

    assert all(isinstance(tm, SecondTM) for tm in drift.tms)
    assert all(isinstance(tm, CavityTM) for tm in cavity.tms)
    assert all(isinstance(tm, TransferMap) for tm in cavity.first_order_tms)


def test_latticeio_element_definition_reads_state_from_wrapper_element():
    quad = Quadrupole(l=0.4, k1=1.2, eid="Q3")
    quad.name = "q3"
    quad.k1 = -0.6

    line = LatticeIO.element_def_string(quad)

    assert line.startswith("q3 = Quadrupole(")
    assert "l=0.4" in line
    assert "k1=-0.6" in line
    assert "eid='Q3'" in line
