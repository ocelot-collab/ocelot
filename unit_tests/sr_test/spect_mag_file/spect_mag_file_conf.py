"""Test parameters description file"""

import pytest
import os

from ocelot import *
from ocelot.rad import *


"""Lattice elements definition"""

mag_file = os.path.dirname(os.path.abspath(__file__)) + '/mag_file.txt'
und = Undulator(field_file=mag_file, eid="und")


"""pytest fixtures definition"""
    
@pytest.fixture(scope='module')
def lattice():
    return MagneticLattice((und))


@pytest.fixture(scope='function')
def beam():

    b = Beam()
    b.E = 1.25
    b.I = 0.1
    b.beta_x = 12.84
    b.beta_y = 6.11
    b.Dx = 0.526

    return b


@pytest.fixture(scope='function')
def screen():

    s = Screen()
    s.z = 50.0
    s.size_x = 0.0
    s.size_y = 0.0
    s.nx = 1
    s.ny = 1
    s.start_energy = 0.0001
    s.end_energy = 0.1
    s.num_energy = 101

    return s
