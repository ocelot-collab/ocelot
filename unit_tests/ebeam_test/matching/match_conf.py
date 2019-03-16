"""Test parameters description file"""

import pytest

from ocelot import *

"""Lattice elements defenition"""

# elements
d0 = Drift(l=0.1)
sol1 = Solenoid(l=0.2, k=4.5)
d1 = Drift(l=0.1)
start = Marker()
end = Marker()



"""pytest fixtures defenition"""

@pytest.fixture(scope='module')
def cell():
    return (start, d0, sol1, d1, end)


@pytest.fixture(scope='module')
def method():
    return MethodTM()
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
