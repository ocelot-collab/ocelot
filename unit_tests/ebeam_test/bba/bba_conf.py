"""Test parameters description file"""

from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.transfer_map import TransferMap


import pytest
import sase1


@pytest.fixture(scope='module')
def cell():
    return sase1.cell


@pytest.fixture(scope='module')
def method():
    return {'global': TransferMap}
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
