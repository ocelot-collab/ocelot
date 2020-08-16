"""Test parameters description"""

import pytest
import numpy as np
import copy

from ocelot import *

"""lattice elements description"""

Q1 = Quadrupole(l=0.2, k1=5.0)
Q2 = Quadrupole(l=0.2, k1=-5.0)

Db = Drift(l=2.0)
Dc = Drift(l=3/2.0)

angle = 45.0 * pi / 180.0
phi = Q1.l * np.sqrt(Q1.k1)
Lc = 2.0 * Dc.l + Q2.l
ro = (1.0 / np.sqrt(Q1.k1) * (Lc * np.sqrt(Q1.k1) * np.cos(phi) + 2.0 * np.sin(phi)) / (Lc * np.sqrt(Q1.k1) * np.sin(phi) - 2.0 *np.cos(phi)) - Db.l) / np.tan(angle/2.0)

B1 = SBend(l=ro*angle, angle=-angle)
B2 = SBend(l=ro*angle, angle=angle)


"""pytest fixtures descripteion"""

@pytest.fixture(scope='module')
def cell():
    cell = (B1, Db, Q1, Dc, Q2, Dc, Q1, Db, B2)
    return [copy.deepcopy(cell), copy.deepcopy(cell)]


@pytest.fixture(scope='module')
def method():

    mmm1 = MethodTM()

    mmm2 = MethodTM()
    mmm2.global_method = SecondTM

    return [mmm1, mmm2]

    
@pytest.fixture(scope='module')
def lattice(cell, method):

    result = []
    for i in range(2):
        result.append(MagneticLattice(cell[i], method=method[i]))
        
    return result
