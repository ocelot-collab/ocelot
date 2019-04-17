"""Test parameters description file"""

import pytest

from ocelot import *

"""Lattice elements defenition"""

# elements
d0 = Drift(l=0.1)
sol1 = Solenoid(l=0.2, k=3)
d1 = Drift(l=0.1)
m_sol = Marker()
q1 = Quadrupole(l=0.2, k1=1)
d2 = Drift(l=0.1)
q2 = Quadrupole(l=0.2, k1=-1)
d3 = Drift(l=0.1)
b1 = Bend(l=0.2, k1=1, angle=0.01)
d4 = Drift(l=0.2)
b2 = Bend(l=0.2, k1=-1, angle=-0.01)
d5 = Drift(l=0.2)
start = Marker()
end = Marker()


"""pytest fixtures defenition"""

@pytest.fixture(scope='module')
def cell():
    return (start, d0, sol1, d1, m_sol, q1, d2, q2, d3, b1, d4, b2, d5, end)


@pytest.fixture(scope='module')
def method():
    return MethodTM()
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
