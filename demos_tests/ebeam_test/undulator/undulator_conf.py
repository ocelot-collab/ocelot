"""Test parameters description"""

import pytest
import numpy as np
import copy

from ocelot import *

"""lattice elements description"""

Q1 = Quadrupole(l=0.3, k1=3.0)
Q2 = Quadrupole(l=0.3, k1=-3.0)

D = Drift(l=0.5)

und = Undulator(Kx=2, nperiods=100, lperiod=0.01, eid="und")
und.ax = 0.1


"""pytest fixtures descripteion"""

@pytest.fixture(scope='module')
def cell():
    cell = (D, Q1, D, und, D, Q2, D)
    return [copy.deepcopy(cell), copy.deepcopy(cell)]


@pytest.fixture(scope='module')
def method():

    mmm1 = MethodTM()
    mmm1.params[Undulator] = UndulatorTestTM

    mmm2 = MethodTM()

    return [mmm1, mmm2]


@pytest.fixture(scope='module')
def lattice(cell, method):
    
    result = []
    for i in range(2):
        result.append(MagneticLattice(cell[i], method=method[i]))
        
    return result
