"""Test parameters description file"""

import pytest

from ocelot import *

"""Lattice elements defenition"""

Q1 = Quadrupole(l=0.4, k1=-1.3, eid="Q1")
Q2 = Quadrupole(l=0.8, k1=1.4, eid="Q2")
Q3 = Quadrupole(l=0.4, k1=-1.7, eid="Q3")
Q4 = Quadrupole(l=0.5, k1=1.4, eid="Q4")

B = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., eid= "B")

SF = Sextupole(l=0.01, k2=1.5, eid="SF") #random value
SD = Sextupole(l=0.01, k2=-1.5, eid="SD") #random value

D1 = Drift(l=2., eid= "D1")
D2 = Drift(l=0.6, eid= "D2")
D3 = Drift(l=0.3, eid= "D3")
D4 = Drift(l=0.7, eid= "D4")
D5 = Drift(l=0.9, eid= "D5")
D6 = Drift(l=0.2, eid= "D6")


"""pytest fixtures defenition"""

@pytest.fixture(scope='module')
def cell():
    return (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)


@pytest.fixture(scope='module')
def method():
    return {'global': TransferMap}
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
