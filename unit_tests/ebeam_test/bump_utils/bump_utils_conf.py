"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *
from ocelot.utils.bump_utils import *

"""Lattice elements definition"""

d = Drift(l=0.35)
d1 = Drift(l=0.6)
qf = Quadrupole(l=0.2, k1=4)
qd = Quadrupole(l=0.2, k1=-4)

c1 = Vcor(l=0.1)
c2 = Vcor(l=0.1)
c3 = Vcor(l=0.1)
c4 = Vcor(l=0.1)

m = Marker()
"""pytest fixtures definition"""

@pytest.fixture(scope='module')
def cell():
    return (d, qf, c1, d1, qd, c2, d1, m, qf, c3, d1, qd, c4, d1, qf, d)


@pytest.fixture(scope='module')
def method():
    
    m = MethodTM()
    m.global_method = TransferMap
    
    return m
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)

