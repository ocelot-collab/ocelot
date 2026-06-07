"""Test parameters description file"""

import pytest
import sase1
from ocelot import *


@pytest.fixture(scope='module')
def cell():
    return sase1.cell


@pytest.fixture(scope='module')
def method():
    return {'global': TransferMap}
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
