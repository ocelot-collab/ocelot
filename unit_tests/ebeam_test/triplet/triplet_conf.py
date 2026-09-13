"""Test parameters description"""

from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.quadrupole import Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


import pytest
import numpy as np
import copy


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

    mmm1 = {'global': TransferMap}

    mmm2 = {'global': SecondTM}

    return [mmm1, mmm2]
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):

    result = []
    for i in range(2):
        result.append(MagneticLattice(cell[i], method=method[i]))
        
    return result
