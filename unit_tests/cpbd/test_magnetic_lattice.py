import pytest


from ocelot.cpbd.magnetic_lattice import MagneticLattice, insert_markers_by_predicate
from ocelot.cpbd.elements import Drift, SBend, Marker


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


def test_magnetic_lattice_insert_markers_by_predicate():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)

    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    def pred(element):
        return element.l == 0.2 and isinstance(element, SBend)

    before_suffix = "before-suffix"
    after_suffix = "after-suffix"    
    insert_markers_by_predicate(sequence, pred,
                                before_suffix=before_suffix,
                                after_suffix=after_suffix)
    isb = sequence.index(sb)
    before = sequence[isb - 1]
    after = sequence[isb + 1]    



    assert isinstance(before, Marker)
    assert isinstance(after, Marker)    
    assert before.id == f"{sb.id}{before_suffix}"
    assert after.id == f"{sb.id}{after_suffix}"    

    
