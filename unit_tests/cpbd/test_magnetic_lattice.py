import pytest


from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.elements import Drift, SBend


def test_magnetic_lattice_find_indices():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)
    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    indices = mlat.find_indices(SBend)
    assert indices == [1, 3]


def test_magnetic_lattice_find_indices_by_predicate():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)

    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    def pred(element):
        return element.l == 0.2 and isinstance(element, SBend)

    indices = mlat.find_indices_by_predicate(pred)
    assert indices == [1, 3]
