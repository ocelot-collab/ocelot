"""Test parameters description"""

import pytest
import numpy as np
import copy

from ocelot import *

"""lattice elements descripteion"""

Q1 = Quadrupole(l=0.3, k1=5.0)
Q2 = Quadrupole(l=0.3, k1=-5.0)

D = Drift(l=0.5)


"""pytest fixtures descripteion"""

@pytest.fixture(scope='module')
def cell():
    cell = (D, Q1, D, Q2, D, Q1, D)
    return [copy.deepcopy(cell), copy.deepcopy(cell)]


@pytest.fixture(scope='module')
def method():

    mmm1 = MethodTM()
    mmm1.global_method = TransferMap

    mmm2 = MethodTM()
    mmm2.global_method = SecondTM

    return [mmm1, mmm2]
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):

    result = []
    for i in range(2):
        result.append(MagneticLattice(cell[i], method=method[i]))
        
    return result
